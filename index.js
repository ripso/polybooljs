/*
 * @copyright 2016 Sean Connelly (@voidqk), http://syntheti.cc
 * @license MIT
 * @preserve Project Home: https://github.com/voidqk/polybooljs
 */

var BuildLog = require('./lib/build-log');
var Epsilon = require('./lib/epsilon');
var Intersecter = require('./lib/intersecter');
var SegmentChainer = require('./lib/segment-chainer');
var SegmentSelector = require('./lib/segment-selector');
var GeoJSON = require('./lib/geojson');
var Convex = require('./lib/convex');
var Simplify = require('./lib/simplify');
var earcut = require('earcut');

var buildLog = false;
var epsilon = Epsilon();

var PolyBool;
PolyBool = {
	// getter/setter for buildLog
	buildLog: function(bl){
		if (bl === true)
			buildLog = BuildLog();
		else if (bl === false)
			buildLog = false;
		return buildLog === false ? false : buildLog.list;
	},
	// getter/setter for epsilon
	epsilon: function(v){
		return epsilon.epsilon(v);
	},

	// core API
	segments: function(poly){
		var i = Intersecter(true, epsilon, buildLog);
		poly.regions.forEach(i.addRegion);
		return {
			segments: i.calculate(poly.inverted),
			inverted: poly.inverted
		};
	},
	combine: function(segments1, segments2){
		var i3 = Intersecter(false, epsilon, buildLog);
		return {
			combined: i3.calculate(
				segments1.segments, segments1.inverted,
				segments2.segments, segments2.inverted
			),
			inverted1: segments1.inverted,
			inverted2: segments2.inverted
		};
	},
	selectUnion: function(combined){
		return {
			segments: SegmentSelector.union(combined.combined, buildLog),
			inverted: combined.inverted1 || combined.inverted2
		}
	},
	selectIntersect: function(combined){
		return {
			segments: SegmentSelector.intersect(combined.combined, buildLog),
			inverted: combined.inverted1 && combined.inverted2
		}
	},
	selectDifference: function(combined){
		return {
			segments: SegmentSelector.difference(combined.combined, buildLog),
			inverted: combined.inverted1 && !combined.inverted2
		}
	},
	selectDifferenceRev: function(combined){
		return {
			segments: SegmentSelector.differenceRev(combined.combined, buildLog),
			inverted: !combined.inverted1 && combined.inverted2
		}
	},
	selectXor: function(combined){
		return {
			segments: SegmentSelector.xor(combined.combined, buildLog),
			inverted: combined.inverted1 !== combined.inverted2
		}
	},
	polygon: function(segments){
		return {
			regions: SegmentChainer(segments.segments, epsilon, buildLog),
			inverted: segments.inverted
		};
	},

	// GeoJSON converters
	polygonFromGeoJSON: function(geojson){
		return GeoJSON.toPolygon(PolyBool, geojson);
	},
	polygonToGeoJSON: function(poly){
		return GeoJSON.fromPolygon(PolyBool, epsilon, poly);
	},

	// helper functions for common operations
	union: function(poly1, poly2){
		return operate(poly1, poly2, PolyBool.selectUnion);
	},
	intersect: function(poly1, poly2){
		return operate(poly1, poly2, PolyBool.selectIntersect);
	},
	difference: function(poly1, poly2){
		return operate(poly1, poly2, PolyBool.selectDifference);
	},
	differenceRev: function(poly1, poly2){
		return operate(poly1, poly2, PolyBool.selectDifferenceRev);
	},
	xor: function(poly1, poly2){
		return operate(poly1, poly2, PolyBool.selectXor);
	},

	// Convex
	isConvex: function(poly) {
		return poly.regions.every(function(poly) { return Convex.isConvex(poly); });
	},
	makeConvex: function(poly) {
		var regions = [];
		poly.regions.forEach(function (polyRegion) {
			var transformed = [];
			polyRegion.forEach(function (point) {
				transformed.push(point[0]);
				transformed.push(point[1]);
			});
			var indices = earcut(transformed);
			for (var i = 0; i < indices.length; i += 3)
				regions.push([polyRegion[indices[i]], polyRegion[indices[i+1]], polyRegion[indices[i+2]]]);
		});
		return {
			regions: regions,
			inverted: false
		}
	},

	// Simplify
	/**
	 * Removes polygons defining negative space (polygons completely contained within other polygons)
	 * @param {*} poly The PolyBool poly structure to subtract from
	 * @param {number[][]} region A polygon to subtract from the (presumably) larger poly
	 * @returns {*} The PolyBool poly structure resulting from the operation
	 */
	differenceAndSimplify: function(poly, region) {
		// Make regions convex
		poly = this.makeConvex(poly);

		// Check for region entirely contained within poly regions
		for (var i = 0; i < poly.regions.length; ++i) {
			var polyRegion = poly.regions[i];
			if (Simplify.regionContains(polyRegion, region)) {
				// Divide manually
				Simplify.carve(polyRegion, region);
				return poly;
			}
		};

		// Apply difference per region
		var result = { regions: [], inverted: false };
		for (var i = 0; i < poly.regions.length; ++i) {
			var polyRegion = poly.regions[i];
			Array.prototype.push.apply(
				result.regions,
				this.polygon(this.selectDifference(this.combine(this.regionSegments(polyRegion), this.regionSegments(region)))).regions
			);
		}
		return result;
	},
	
	/**
	 * Faster segment creation for a single non-inverted region
	 * @param {number[][]} region The region
	 * @returns {*} The segments for the region
	 */
	regionSegments: function(region) {
		var i = Intersecter(true, epsilon, buildLog);
		i.addRegion(region);
		return {
			segments: i.calculate(false),
			inverted: false
		};
	}
};

function operate(poly1, poly2, selector){
	var seg1 = PolyBool.segments(poly1);
	var seg2 = PolyBool.segments(poly2);
	var comb = PolyBool.combine(seg1, seg2);
	var seg3 = selector(comb);
	return PolyBool.polygon(seg3);
}

if (typeof window === 'object')
	window.PolyBool = PolyBool;

module.exports = PolyBool;
