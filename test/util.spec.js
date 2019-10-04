
var Util = require('../lib/util');

describe('Util', () => {
    describe('#normalizeIndex', () => {
        var array = [5,6,7,8,9];

        it('should return the value passed in for indices in range', () => {
            expect(Util.normalizeIndex(array, 3)).toBe(3);
        });

        it('should wrap indices < 0', () => {
            expect(Util.normalizeIndex(array, -1)).toBe(4);
        });

        it('should wap indices < -length', () => {
            expect(Util.normalizeIndex(array, -7)).toBe(3);
        });

        it('should wrap indices > length', () => {
            expect(Util.normalizeIndex(array, 6)).toBe(1);
        });

        it('should wrap indices > 2*length', () => {
            expect(Util.normalizeIndex(array, 12)).toBe(2);
        });
    });
});
