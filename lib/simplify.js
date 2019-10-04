
var Util = require('./util');

/**
 * The direction of the point relative to the line segment as determined by the cross-product
 * @param {number[]} segmentStart 
 * @param {number[]} segmentEnd 
 * @param {number[]} point 
 * @returns {number} right (>0), left (<0), or straight (0)
 */
var pointDirection = function(segmentStart, segmentEnd, point) {
    var segment = {
        x: segmentEnd[0] - segmentStart[0],
        y: segmentEnd[1] - segmentStart[1]
    };
    var p = {
        x: point[0] - segmentStart[0],
        y: point[1] - segmentStart[1]
    }
    return segment.x * p.y - segment.y * p.x;
}

var distanceSquared = function(p1, p2) {
    var dx = p2[0] - p1[0];
    var dy = p2[1] - p1[1];
    return dx * dx + dy * dy;
}

var Simplify = {
    /**
     * Whether subpoly is entirely contained within poly
     * @param {number[][]} poly A convex poly
     * @param {number[][]} subpoly 
     * @returns {boolean}
     */
    regionContains: function(poly, subpoly) {
        for (var i = 0; i < subpoly.length; ++i) {
            var subpolyPoint = subpoly[i];
            var firstSubpolyPointDirection = pointDirection(poly[poly.length - 1], poly[0], subpolyPoint);
            for (var j = 1; j < poly.length; ++j) {
                var subpolyPointDirection = pointDirection(poly[j - 1], poly[j], subpolyPoint);
                if (firstSubpolyPointDirection === 0) {
                    firstSubpolyPointDirection = subpolyPointDirection;
                } else if (subpolyPointDirection !== 0 && (subpolyPointDirection > 0) !== (firstSubpolyPointDirection > 0)) {
                    return false;
                }
            }
        }
        return true;
    },

    /**
     * Carve a poly completely contained within another convex poly. Modifies the larger region.
     * @param {number[][]} region a convex poly
     * @param {number[][]} subregion
     */
    carve: function(region, subregion) {
        var targetPoint = region[region.length - 1];
        var closestPoint, closestPointDistance = Infinity;
        for (var i = 0; i < subregion.length; ++i) {
            var distance = distanceSquared(targetPoint, subregion[i]);
            if (distance < closestPointDistance) {
                closestPointDistance = distance;
                closestPoint = i;
            }
        }

        var regionWinding = pointDirection(region[0], region[1], region[2]) > 0;
        var subregionWinding = pointDirection(subregion[0], subregion[1], subregion[2]) > 0
        var backward = regionWinding === subregionWinding;
        for (var i = 0; i < subregion.length; ++i) {
            region.push(subregion[Util.normalizeIndex(subregion, closestPoint + (backward ? -i : i))]);
        }
        region.push(subregion[closestPoint]);
        region.push(targetPoint);
    }
};

module.exports = Simplify;
