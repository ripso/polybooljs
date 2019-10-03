
var Util = {
    normalizeIndex: function (array, index) {
        while (index < 0)
            index += array.length;
        while (index >= array.length)
            index -= array.length;
        return index;
    }
};

module.exports = Util;
