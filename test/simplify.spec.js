
var Convex = require('../lib/convex');
var Simplify = require('../lib/simplify');

describe('Simplify', () => {
    var polys = [
        [[0,0],[10,0],[10,10],[0,10]],
        [[0,0],[0,10],[10,10],[10,0]],
        [[4,2],[8,8],[2,5]],
        [[4,-2],[8,8],[2,5]]
    ];
    var insideTestPoints = [
        [1,1],
        [1,7],
        [3,7],
        [6,9],
        [5,7]
    ];
    var outsideTestPoints = [
        [-1,4],
        [11,7],
        [4,5],
    ];
    var sometimesOutsideTestPoints = [
        [4,0],
        [4.1,2]
    ];

    var regionsContain = function(regions, point) {
        return regions.some((region) => Simplify.regionContains(region, [point]));
    }

    describe('#regionContains', () => {
        it('should return true if a poly is contained within another poly', () => {
            expect(Simplify.regionContains(polys[0], polys[2])).toBe(true);
        });

        it('shouldn\'t care about vertex winding order', () => {
            expect(Simplify.regionContains(polys[1], polys[2])).toBe(true);
        });

        it('should return false if a poly is not contained within another poly', () => {
            expect(Simplify.regionContains(polys[0], polys[3])).toBe(false);
        });
    });

    describe('#carve', () => {
        it('should properly difference polys completely contained within other polys', () => {
            var poly = [...polys[0]];
            Simplify.carve(poly, polys[2]);
            var convexPoly = Convex.makeConvex(poly);

            expect(insideTestPoints.every((point) => regionsContain(convexPoly, point))).toBe(true);
            expect(outsideTestPoints.some((point) => regionsContain(convexPoly, point))).toBe(false);
        });
    });
});
