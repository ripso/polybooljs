
var Util = require('./util');

var TAU = Math.PI * 2;
var HALF_PI = Math.PI * 0.5;

var makePoly = function(poly, index1, index2) {
    var result = [];
    for (var i = index1; i !== index2 + 1; ++i) {
        i = Util.normalizeIndex(poly, i);
        result.push(poly[i]);
    }
    return result;
}
var splitPoly = function(poly, index1, index2) {
    return [ makePoly(poly, index1, index2), makePoly(poly, index2, index1) ];
}

var getPolyTurns = function(poly) {
    // Check for too few points
    if (poly.length < 3)
        return 0;

    var oldPoly = poly[poly.length - 2];
    var newPoly = poly[poly.length - 1];
    var oldX = oldPoly[0];
    var oldY = oldPoly[1];
    var newX = newPoly[0];
    var newY = newPoly[1];
    var newDirection = Math.atan2(newY - oldY, newX - oldX);
    var angleSum = 0;

    // Check each point, the side ending there, its angle, and accumulate angles
    for (var ndx = 0; ndx < poly.length; ++ndx) {
        var newPoint = poly[ndx];

        // Update point coordinates and side direction, check side length
        oldX = newX;
        oldY = newY;
        if (oldX === newPoint[0] && oldY === newPoint[1])
            continue;
        var oldDirection = newDirection;
        newX = newPoint[0];
        newY = newPoint[1];
        newDirection = Math.atan2(newY - oldY, newX - oldX);

        // Calculate and check the normalized direction change angle
        var angle = newDirection - oldDirection;
        while (angle <= -Math.PI)
            angle += TAU;
        while (angle > Math.PI)
            angle -= TAU;

        // Accumulate the direction change angle
        angleSum += angle;
    }

    // Check that the total number of full turns is plus or minus 1
    return Math.round(angleSum / TAU);
}

var normalizeAngle = function(angle) {
    while (angle <= -Math.PI)
        angle += TAU;
    while (angle > Math.PI)
        angle -= TAU;
    return angle;
}

var Convex = {
    /**
     * Check if a polygon is convex
     * @param {number[][]} poly - An array of points defining a polygon
     * @returns {boolean} True if the polygon is convex
     */
    isConvex: function(poly) {
        // From https://stackoverflow.com/questions/471962/how-do-i-efficiently-determine-if-a-polygon-is-convex-non-convex-or-complex/45372025#45372025

        // Check for too few points
        if (poly.length < 3)
            return false;

        var oldPoly = poly[poly.length - 2];
        var newPoly = poly[poly.length - 1];
        var oldX = oldPoly[0];
        var oldY = oldPoly[1];
        var newX = newPoly[0];
        var newY = newPoly[1];
        var newDirection = Math.atan2(newY - oldY, newX - oldX);
        var angleSum = 0;
        var orientation;

        // Check each point, the side ending there, its angle, and accumulate angles
        for (var ndx = 0; ndx < poly.length; ++ndx) {
            var newPoint = poly[ndx];

            // Update point coordinates and side direction, check side length
            oldX = newX;
            oldY = newY;
            if (oldX === newPoint[0] && oldY === newPoint[1])
                continue;
            var oldDirection = newDirection; 
            newX = newPoint[0];
            newY = newPoint[1];
            newDirection = Math.atan2(newY - oldY, newX - oldX);
            
            // Calculate and check the normalized direction change angle
            var angle = normalizeAngle(newDirection - oldDirection);
            if (typeof orientation === 'undefined') {
                // First time through the loop, initialize orientation
                if (angle === 0)
                    continue;
                orientation = angle > 0 ? 1 : -1;
            } else if (orientation * angle < 0) { // Check orientation is stable
                return false;
            }
            // Accumulate the direction change angle
            angleSum += angle;
        }

        // Check that the total number of full turns is plus or minus 1
        return Math.abs(Math.round(angleSum / TAU)) === 1;
    },

    /**
     * Split a single polygon into multiple convex polygons. Does not handle complex polygons (those that intersect themselves).
     * @param {number[][]} poly - An array of points defining a polygon
     * @returns {number[][][]} An array of convex polygons
     */
    makeConvex: function(poly) {
        // Check for too few points
        if (poly.length < 3)
            return [];
        else if (poly.length === 3)
            return [poly];

        var orientation = getPolyTurns(poly);
        if (Math.abs(orientation) !== 1)
            throw 'Convex.makeConvex does not support complex polygons. Poly given was: ' + JSON.stringify(poly);

        var oldPoly = poly[poly.length - 2];
        var newPoly = poly[poly.length - 1];
        var oldX = oldPoly[0];
        var oldY = oldPoly[1];
        var newX = newPoly[0];
        var newY = newPoly[1];
        var newDirection = Math.atan2(newY - oldY, newX - oldX);
        var concavePointIndex;
        var concavePointDirection;
        var foundConvexPoint = false;
        var firstConvexPointIndex = 0;

        // Check each point, the side ending there, its angle, and accumulate angles
        for (var ndx = 0; ndx - firstConvexPointIndex < poly.length; ++ndx) {
            var newPoint = poly[Util.normalizeIndex(poly, ndx)];

            // Update point coordinates and side direction, check side length
            oldX = newX;
            oldY = newY;
            if (oldX === newPoint[0] && oldY === newPoint[1])
                continue;
            var oldDirection = newDirection;
            newX = newPoint[0];
            newY = newPoint[1];
            newDirection = Math.atan2(newY - oldY, newX - oldX);

            // Calculate and check the normalized direction change angle
            var angle = normalizeAngle(newDirection - oldDirection);

            // Check orientation is stable
            if (orientation * angle >= 0 && (typeof concavePointDirection === 'undefined' || Math.abs(normalizeAngle(oldDirection - concavePointDirection)) < HALF_PI)) {
                // Convex point
                if (!foundConvexPoint) {
                  foundConvexPoint = true;
                  firstConvexPointIndex = ndx;
                }
            } else if (foundConvexPoint) { // Don't start looking until we're on a convex surface
                // Concave point
                if (typeof concavePointIndex === 'undefined') {
                    if (foundConvexPoint) {
                        concavePointIndex = ndx - 1;
                        concavePointDirection = newDirection;
                    }
                } else {
                    var nextConcavePointIndex = ndx - 1;
                    if (nextConcavePointIndex === concavePointIndex + 1) {
                        concavePointIndex = nextConcavePointIndex;
                    } else {
                        var index1 = Util.normalizeIndex(poly, concavePointIndex);
                        var index2 = Util.normalizeIndex(poly, ndx - 1);
                        var split = splitPoly(poly, index1, index2);
                        return [split[0]].concat(this.makeConvex(split[1]));
                    }
                }
            }
        }

        if (typeof concavePointIndex !== 'undefined') {
            // Break up poly
            var index1 = Util.normalizeIndex(poly, concavePointIndex);
            var index2 = Util.normalizeIndex(poly, concavePointIndex + 2);
            var split = splitPoly(poly, index1, index2);
            return [split[0]].concat(this.makeConvex(split[1]));
        }

        return [poly];
    }
};

/*
if (window) {
    window.debugPolyToHtmlPolygon = function(poly, maxWidth, maxHeight) {
        var minX = Infinity;
        var maxX = -Infinity;
        var minY = Infinity;
        var maxY = -Infinity;

        poly.forEach(function(point) {
            minX = Math.min(minX, point[0]);
            maxX = Math.max(maxX, point[0]);
            minY = Math.min(minY, point[1]);
            maxY = Math.max(maxY, point[1]);
        });

        var scaleX = maxWidth / (maxX - minX);
        var scaleY = maxHeight / (maxY - minY);
        var scale = Math.min(scaleX, scaleY);

        var result = '<polygon fill="lime" stroke="purple" points="';
        result += poly.map(function(point) { return ((point[0] - minX) * scale) + ',' + ((point[1] - minY) * scale); }).join(' ');
        result += '"/>';
        return result;
    }
}
*/

module.exports = Convex;
