
var Convex = require('../lib/convex');

describe('Convex', () => {
    describe('#isConvex', () => {
        var convexPolys = [
            [[0,0],[1,0],[1,1],[0,1]],
            [[0,0],[0,1],[.5,.5]],
            [[0,0],[1,0],[2,2],[3,5],[2,8],[0,6]],
            [[0,0],[1,0],[2,2],[3,5],[2,8],[0,6],[0,2]]
        ];

        var concavePolys = [
            [[0,0],[1,0],[1,1],[0.6,.5]],
            [[0,0],[1,0],[2,2],[3,5],[2,8],[0,6]],
            [[0,0],[1,0],[2,2],[3,5],[2,8],[0,6],[.1,2]]
        ];

        it('should properly identify convex polygons', () => {
            expect(convexPolys.every((poly) => Convex.isConvex(poly))).toBe(true);
        });

        it('should properly identify concave polygons', () => {
            expect(concavePolys.every((poly) => Convex.isConvex(poly))).toBe(false);
        });
    });
});
