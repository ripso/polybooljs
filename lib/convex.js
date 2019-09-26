

var TAU = Math.PI * 2;

var normalizeIndex = function(array, index) {
    while (index < 0)
        index += array.length;
    while (index >= array.length)
        index -= array.length;
    return index;
}

var makePoly = function(poly, index1, index2) {
    var result = [];
    for (var i = index1; i !== index2 + 1; ++i) {
        i = normalizeIndex(poly, i);
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
            var angle = newDirection - oldDirection;
            if (angle <= -Math.PI)
                angle += TAU;
            else if (angle > Math.PI)
                angle -= TAU;
            if (typeof (orientation) === 'undefined') {
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
        var convexPointIndex;

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
            if (angle <= -Math.PI)
                angle += TAU;
            else if (angle > Math.PI)
                angle -= TAU;

            if (orientation * angle < 0) { // Check orientation is stable
                // Convex point
                if (typeof(convexPointIndex) === 'undefined') {
                    convexPointIndex = normalizeIndex(poly, ndx - 1);
                } else {
                    var split = splitPoly(poly, convexPointIndex, normalizeIndex(poly, ndx - 1));
                    return [split[0]].concat(this.makeConvex(split[1]));
                }
            }
        }

        if (typeof(convexPointIndex) !== 'undefined') {
            // Break up poly
            return splitPoly(poly, convexPointIndex, normalizeIndex(poly, Math.floor(convexPointIndex + poly.length * 0.5)));
        }

        return [poly];
    }
};

module.exports = Convex;
