
var Convex = require('../lib/convex');

describe('Convex', () => {
    var convexPolys = [
        [[0, 0], [1, 0], [1, 1], [0, 1]],
        [[0, 0], [0, 1], [.5, .5]],
        [[0, 0], [1, 0], [2, 2], [3, 5], [2, 8], [0, 6]],
        [[0, 0], [1, 0], [2, 2], [3, 5], [2, 8], [0, 6], [0, 2]]
    ];

    var concavePolys = [
        [[0, 0], [1, 0], [1, 1], [0.6, .5]],
        [[0, 0], [1, 0], [2, 2], [3, 5], [2, 8], [0, 6]],
        [[0, 0], [1, 0], [2, 2], [3, 5], [2, 8], [0, 6], [.1, 2]],
        [[400, 200], [400, 174], [396, 177], [383, 174], [381, 174], [372, 174], [368, 170], [356, 170], [353, 166], [350, 166], [350,200]]
    ];

    describe('#isConvex', () => {
        it('should properly identify convex polygons', () => {
            expect(convexPolys.every((poly) => Convex.isConvex(poly))).toBe(true);
        });

        it('should properly identify concave polygons', () => {
            expect(concavePolys.every((poly) => Convex.isConvex(poly))).toBe(false);
        });
    });

    describe('#makeConvex', () => {
        it('should break down concave polygons into convex polygons', () => {
            expect(concavePolys.every((poly) => Convex.makeConvex(poly).every((poly) => Convex.isConvex(poly)))).toBe(true);
        })
    });
});
