/**
 * @license Quaternion.js v0.0.0 22/02/2016
 *
 * Copyright (c) 2016, Robert Eisele (robert@xarg.org)
 * Dual licensed under the MIT or GPL Version 2 licenses.
 **/
(function(root) {

  'use strict';

  var EPSILON = 1e-16;

  Math.hypot = Math.hypot || function() {

    var sum = 0;
    for (var i = 0; i < arguments.length; i++) {
      sum += arguments[i] * arguments[i];
    }
    return Math.sqrt(sum);
  };

  /*
   * Default is the multiplicative one element
   *
   * Proof: TODO
   */
  var P = {
    'w': 1,
    'x': 0,
    'y': 0,
    'z': 0
  };

  function parse(dest, w, x, y, z) {

    // Most common internal use case with 4 params
    if (z !== undefined) {
      dest['w'] = w;
      dest['x'] = x;
      dest['y'] = y;
      dest['z'] = z;
      return;
    }

    if (typeof w === 'object') {

      // Check for quats, for example when an object gets cloned
      if ('w' in w && 'x' in w && 'y' in w && 'z' in w) {
        dest['w'] = w['w'];
        dest['x'] = w['x'];
        dest['y'] = w['y'];
        dest['z'] = w['z'];
        return;
      }

      // Check for complex numbers
      if ('re' in w && 'im' in w) {
        dest['w'] = w['re'];
        dest['x'] = w['im'];
        dest['y'] = 0;
        dest['z'] = 0;
        return;
      }

      // Check for arrays
      if (w.length === 4) {
        if (w[0].length === 4) {
          if (checkFor4x4Matrix(w)) {
            dest['w'] = w[0][0];
            dest['x'] = -w[0][1];
            dest['y'] = -w[0][2];
            dest['z'] = -w[0][3];
          } else {
            throw "cant convert matrix to quaternion";
          }
        } else {
          dest['w'] = w[0];
          dest['x'] = w[1];
          dest['y'] = w[2];
          dest['z'] = w[3];
        }
        return;
      }

      // if a modulus, argument and unit vecor as passed in
      if (typeof w === 'object' && 'mod' in w && 'arg' in w && 'unit' in w &&
        ('x' in w.unit && 'y' in w.unit && 'z' in w.unit && ( !('w' in w.unit) || w.unit.w === 0))) {
        /*
         * quaternion can be represented in modulus argument form like:
         *
         * q = |q|(sin(t) + n * cos(t))
         *
         * where:
         *  q = a + bi + cj + dk = a + v
         *  n = v/|v|
         *  t = arccos(a/|q|)
         *
         * here the mod is passed in as a, the arg as b and the unit vector as c
         */

        dest['w'] = w.mod * Math.cos(w.arg);
        var sinarg = Math.sin(w.arg);
        dest['x'] = sinarg * w.mod * w.unit['x'];
        dest['y'] = sinarg * w.mod * w.unit['y'];
        dest['z'] = sinarg * w.mod * w.unit['z'];
        return
      }

        throw 'Invalid object';
    }

    // Parse string values
    if (typeof w === 'string') {

      var tokens = w.match(/\d+\.?\d*e[+-]?\d+|\d+\.?\d*|\.\d+|./g);
      var plus = 1;
      var minus = 0;

      var iMap = {'i': 'x', 'j': 'y', 'k': 'z'};

      if (tokens === null) {
        throw 'Parse error';
      }

      // Reset the current state
      dest['w'] =
              dest['x'] =
              dest['y'] =
              dest['z'] = 0;

      for (var i = 0; i < tokens.length; i++) {

        var c = tokens[i];
        var d = tokens[i + 1];

        if (c === ' ' || c === '\t' || c === '\n') {
          /* void */
        } else if (c === '+') {
          plus++;
        } else if (c === '-') {
          minus++;
        } else {

          if (plus + minus === 0) {
            throw 'Parse error' + c;
          }
          var g = iMap[c];

          // Is the current token an imaginary sign?
          if (g !== undefined) {

            // Is the following token a number?
            if (d !== ' ' && !isNaN(d)) {
              c = d;
              i++;
            } else {
              c = '1';
            }

          } else {

            if (isNaN(c)) {
              throw 'Parser error';
            }

            g = iMap[d];

            if (g !== undefined) {
              i++;
            }
          }

          dest[g || 'w'] += parseFloat((minus % 2 ? '-' : '') + c);
          plus = minus = 0;
        }
      }

      // Still something on the stack
      if (plus + minus > 0) {
        throw 'Parser error';
      }
      return;
    }

    // If no single variable was given AND it was the constructor, set it to the identity
    if (w === undefined && dest !== P) {
      dest['w'] = 1;
      dest['x'] = 0;
      dest['y'] = 0;
      dest['z'] = 0;
    } else {

      dest['w'] = w || 0;

      if (x && x.length === 3) {
        dest['x'] = x[0];
        dest['y'] = x[1];
        dest['z'] = x[2];
      } else {
        dest['x'] = x || 0;
        dest['y'] = y || 0;
        dest['z'] = z || 0;
      }
    }
  }

  function numToStr(n, char, prev) {

    var ret = '';

    if (n !== 0) {

      if (prev !== '') {
        ret += n < 0 ? ' - ' : ' + ';
      } else if (n < 0) {
        ret += '-';
      }

      n = Math.abs(n);

      if (1 !== n || char === '') {
        ret += n;
      }
      ret += char;
    }
    return ret;
  }

  function checkFor4x4Matrix(a) {
    /*
     * for a real 4x4 matrix to be converted to a quaternion (a + bi + cj + dk) it must take the form:
     *
     *      a -b -c -d
     *      b  a -d  c
     * M =  c  d  a -b
     *      d -c  b  a
     *
     * easiest check is to check that array[i][j] === -array[j][i] unless i==j
     * and that the lead diagonal elements are all equal
     */

    var a;
    if (!Array.isArray(a)) {
      return false
    }

    if (!(Array.isArray(a)) || (a.length !== 4)) {
      return false;
    }

    for (var i = 0; i < 4; i++) {
      if (!(Array.isArray(a[i])) || a[i].length !== 4) {
        return false;
      } else {
        for (var j = 0; j < 4; j++) {
          if (typeof a[i][j] !== 'number') {
            return false;
          }
        }
      }
    }
    return true;
  }

  /**
   * Quaternion constructor
   *
   * @constructor
   * @param {number|Object|string} w real
   * @param {number=} x imag
   * @param {number=} y imag
   * @param {number=} z imag
   * @returns {Quaternion}
   */
  function Quaternion(w, x, y, z) {

    if (!(this instanceof Quaternion)) {
      return new Quaternion(w, x, y, z);
    }

    parse(this, w, x, y, z);
  }

  Quaternion.prototype = {
    'w': 1,
    'x': 0,
    'y': 0,
    'z': 0,
    /**
     * Adds two quaternions Q1 and Q2
     *
     * @param {number|Object|string} w real
     * @param {number=} x imag
     * @param {number=} y imag
     * @param {number=} z imag
     * @returns {Quaternion}
     */
    'add': function(w, x, y, z) {

      parse(P, w, x, y, z);

      // Q1 + Q2 := [w1, v1] + [w2, v2] = [w1 + w2, v1 + v2]

      return new Quaternion(
              this['w'] + P['w'],
              this['x'] + P['x'],
              this['y'] + P['y'],
              this['z'] + P['z']);
    },
    /**
     * Subtracts a quaternions Q2 from Q1
     *
     * @param {number|Object|string} w real
     * @param {number=} x imag
     * @param {number=} y imag
     * @param {number=} z imag
     * @returns {Quaternion}
     */
    'sub': function(w, x, y, z) {

      parse(P, w, x, y, z);

      // Q1 - Q2 := Q1 + (-Q2)
      //          = [w1, v1] - [w2, v2] = [w1 - w2, v1 - v2]

      return new Quaternion(
              this['w'] - P['w'],
              this['x'] - P['x'],
              this['y'] - P['y'],
              this['z'] - P['z']);
    },
    /**
     * Calculates the additive inverse, or simply it negates the quaternion
     *
     * @returns {Quaternion}
     */
    'neg': function() {

      // -Q := [-w, -v]

      return new Quaternion(-this['w'], -this['x'], -this['y'], -this['z']);
    },
    /**
     * Calculates the length/modulus or the norm of a quaternion
     *
     * @returns {number}
     */
    'norm': function() {

      // |Q| := sqrt(|Q|^2)

      // The unit quaternion has |Q| = 1

      var w = this['w'];
      var x = this['x'];
      var y = this['y'];
      var z = this['z'];

      return Math.hypot(w, x, y, z);
    },
    /**
     * Calculates the squared length/modulus or the norm of a quaternion
     *
     * @returns {number}
     */
    'normSq': function() {

      // |Q|^2 := [w, v] * [w, -v]
      //        = [w^2 + dot(v, v), -w * v + w * v + cross(v, -v)]
      //        = [w^2 + |v|^2, 0]
      //        = [w^2 + dot(v, v), 0]
      //        = dot(Q, Q)
      //        = Q * Q'

      var w = this['w'];
      var x = this['x'];
      var y = this['y'];
      var z = this['z'];

      return w * w + x * x + y * y + z * z;
    },
    /**
     * Normalizes the quaternion to have |Q| = 1 as long as the norm is not zero
     * Alternative names are the signum, unit or versor
     *
     * @returns {Quaternion}
     */
    'normalize': function() {

      // Q* := Q / |Q|

      var w = this['w'];
      var x = this['x'];
      var y = this['y'];
      var z = this['z'];

      var norm = Math.hypot(w, x, y, z);

      if (norm === 0) {
        return Quaternion['ZERO']; // TODO: Is the result zero or one when the norm is zero? -> limes
      }

      norm = 1 / norm;

      return new Quaternion(w * norm, x * norm, y * norm, z * norm);
    },
    /**
     * Calculates the Hamilton product of two quaternions
     * Leaving out the imaginary part results in just scaling the quat
     *
     * @param {number|Object|string} w real
     * @param {number=} x imag
     * @param {number=} y imag
     * @param {number=} z imag
     * @returns {Quaternion}
     */
    'mul': function(w, x, y, z) {

      parse(P, w, x, y, z);

      // Q1 * Q2 = [w1 * w2 - dot(v1, v2), w1 * v2 + w2 * v1 + cross(v1, v2)]

      // Not commutative because cross(v1, v2) != cross(v2, v1)!

      var w1 = this['w'];
      var x1 = this['x'];
      var y1 = this['y'];
      var z1 = this['z'];

      var w2 = P['w'];
      var x2 = P['x'];
      var y2 = P['y'];
      var z2 = P['z'];

      return new Quaternion(
              w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
              w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
              w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2,
              w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2);
    },
    /**
     * Scales a quaternion by a scalar, faster than using multiplication
     *
     * @param {number} s scaling factor
     * @returns {Quaternion}
     */
    'scale': function(s) {

      return new Quaternion(
              this['w'] * s,
              this['x'] * s,
              this['y'] * s,
              this['z'] * s);
    },
    /**
     * Calculates the dot product of two quaternions
     *
     * @param {number|Object|string} w real
     * @param {number=} x imag
     * @param {number=} y imag
     * @param {number=} z imag
     * @returns {number}
     */
    'dot': function(w, x, y, z) {

      parse(P, w, x, y, z);

      // dot(Q1, Q2) := w1 * w2 + dot(v1, v2)

      return this['w'] * P['w'] + this['x'] * P['x'] + this['y'] * P['y'] + this['z'] * P['z'];
    },
    /**
     * Calculates the inverse of a quat such that
     * Q^-1 * Q = 1 and Q * Q^-1 = 1
     *
     * @returns {Quaternion}
     */
    'inverse': function() {

      // Q^-1 := Q' / |Q|^2
      //       = [w / (w^2 + |v|^2), -v / (w^2 + |v|^2)]

      // Proof:
      // Q * Q^-1 = [w, v] * [w / (w^2 + |v|^2), -v / (w^2 + |v|^2)]
      //          = [1, 0]
      // Q^-1 * Q = [w / (w^2 + |v|^2), -v / (w^2 + |v|^2)] * [w, v]
      //          = [1, 0].

      var w = this['w'];
      var x = this['x'];
      var y = this['y'];
      var z = this['z'];

      var normSq = w * w + x * x + y * y + z * z;

      if (normSq === 0) {
        return Quaternion['ZERO']; // TODO: Is the result zero or one when the norm is zero?
      }

      normSq = 1 / normSq;

      return new Quaternion(w * normSq, -x * normSq, -y * normSq, -z * normSq);
    },
    /**
     * Multiplies a quaternion with the inverse of a second quaternion
     *
     * @param {number|Object|string} w real
     * @param {number=} x imag
     * @param {number=} y imag
     * @param {number=} z imag
     * @returns {Quaternion}
     */
    'div': function(w, x, y, z) {

      parse(P, w, x, y, z);

      // Q1 / Q2 := Q1 * Q2^-1

      var w1 = this['w'];
      var x1 = this['x'];
      var y1 = this['y'];
      var z1 = this['z'];

      var w2 = P['w'];
      var x2 = P['x'];
      var y2 = P['y'];
      var z2 = P['z'];

      var normSq = w2 * w2 + x2 * x2 + y2 * y2 + z2 * z2;

      if (normSq === 0) {
        return Quaternion['ZERO']; // TODO: Is the result zero or one when the norm is zero?
      }

      normSq = 1 / normSq;

      return new Quaternion(
              (w1 * w2 + x1 * x2 + y1 * y2 + z1 * z2) * normSq,
              (x1 * w2 - w1 * x2 + y1 * z2 - z1 * y2) * normSq,
              (y1 * w2 - w1 * y2 + z1 * x2 - x1 * z2) * normSq,
              (z1 * w2 - w1 * z2 + x1 * y2 - y1 * x2) * normSq);
    },
    /**
     * Calculates the conjugate of a quaternion
     *
     * @returns {Quaternion}
     */
    'conjugate': function() {

      // Q' = [s, -v]

      // If the quaternion is normalized,
      // the conjugate is the inverse of the quaternion - but faster
      // Q' * Q = Q * Q' = 1

      // Additionally, the conjugate of a unit quaternion is a rotation with the same
      // angle but the opposite axis.

      // Moreover the following property holds:
      // (Q1 * Q2)' = Q2' * Q1'

      return new Quaternion(this['w'], -this['x'], -this['y'], -this['z']);
    },
    /**
     * Checks if two quats are the same
     *
     * @param {number|Object|string} w real
     * @param {number=} x imag
     * @param {number=} y imag
     * @param {number=} z imag
     * @returns {boolean}
     */
    'equals': function(w, x, y, z) {

      parse(P, w, x, y, z);

      // maybe check for NaN's here?
      return Math.abs(P['w'] - this['w']) < EPSILON && Math.abs(P['x'] - this['x']) < EPSILON && Math.abs(P['y'] - this['y']) < EPSILON && Math.abs(P['z'] - this['z']) < EPSILON;
    },
    /**
     * Checks if all parts of a quaternion are finite
     *
     * @returns {boolean}
     */
    'isFinite': function() {

      return isFinite(this['w']) && isFinite(this['x']) && isFinite(this['y']) && isFinite(this['z']);
    },
    /**
     * Checks if any of the parts of the quaternion is not a number
     *
     * @returns {boolean}
     */
    'isNaN': function() {

      return isNaN(this['w']) || isNaN(this['x']) || isNaN(this['y']) || isNaN(this['z']);
    },
    /**
     * Gets the Quaternion as a well formatted string
     *
     * @returns {string}
     */
    'toString': function() {

      var w = this['w'];
      var x = this['x'];
      var y = this['y'];
      var z = this['z'];
      var ret = '';

      if (isNaN(w) || isNaN(x) || isNaN(y) || isNaN(z)) {
        return 'NaN';
      }

      // Alternative design?
      // '(%f, [%f %f %f])'

      ret = numToStr(w, '', ret);
      ret += numToStr(x, 'i', ret);
      ret += numToStr(y, 'j', ret);
      ret += numToStr(z, 'k', ret);

      if ('' === ret)
        return '0';

      return ret;
    },
    /**
     * Returns the real part of the quaternion
     *
     * @returns {number}
     */
    'real': function() {

      return this['w'];
    },
    /**
     * Returns the imaginary part of the quaternion as a new quaternion with real part zero
     *
     * @returns {Quaternion}
     */
    'imag': function() {

      return new Quaternion(
              0,
              this['x'],
              this['y'],
              this['z']);
    },
    /**
     * Gets the actual quaternion as an array
     *
     * @returns {Array}
     */
    'toArray': function() {

      return [this['w'], this['x'], this['y'], this['z']];
    },
    /**
     * Calculates the 3x3 rotation matrix for the current quat
     *
     * @see https://en.wikipedia.org/wiki/Rotation_matrix#Quaternion
     * @returns {Array}
     */
    'toMatrix': function() {

      var w = this['w'];
      var x = this['x'];
      var y = this['y'];
      var z = this['z'];

      var n = w * w + x * x + y * y + z * z;
      var s = n === 0 ? 0 : 2 / n;
      var wx = s * w * x, wy = s * w * y, wz = s * w * z;
      var xx = s * x * x, xy = s * x * y, xz = s * x * z;
      var yy = s * y * y, yz = s * y * z, zz = s * z * z;

      return [
        1 - (yy + zz), xy - wz, xz + wy,
        xy + wz, 1 - (xx + zz), yz - wx,
        xz - wy, yz + wx, 1 - (xx + yy)];
    },
    /**
     * Calculates the 4x4 rotation matrix for the current quat
     *
     * @returns {Array}
     */
    'toMatrix4': function() {

      var w = this['w'];
      var x = this['x'];
      var y = this['y'];
      var z = this['z'];

      var n = w * w + x * x + y * y + z * z;
      var s = n === 0 ? 0 : 2 / n;
      var wx = s * w * x, wy = s * w * y, wz = s * w * z;
      var xx = s * x * x, xy = s * x * y, xz = s * x * z;
      var yy = s * y * y, yz = s * y * z, zz = s * z * z;

      return [
        1 - (yy + zz), xy - wz, xz + wy, 0,
        xy + wz, 1 - (xx + zz), yz - wx, 0,
        xz - wy, yz + wx, 1 - (xx + yy), 0,
        0, 0, 0, 1];
    },
    /**
     * Clones the actual object
     *
     * @returns {Quaternion}
     */
    'clone': function() {

      return new Quaternion(this);
    },
    /**
     * Rotates a vector according to the current quaternion
     *
     * @param {Array} v The vector to be rotated
     * @returns {Array}
     */
    'rotateVector': function(v) {

      // [0, v'] = Q * [0, v] * Q'

      // Q
      var w1 = this['w'];
      var x1 = this['x'];
      var y1 = this['y'];
      var z1 = this['z'];

      // [0, v]
      var w2 = 0;
      var x2 = v[0];
      var y2 = v[1];
      var z2 = v[2];

      // Q * [0, v]
      var w3 = /*w1 * w2*/ -x1 * x2 - y1 * y2 - z1 * z2;
      var x3 = w1 * x2 + /*x1 * w2 +*/ y1 * z2 - z1 * y2;
      var y3 = w1 * y2 + /*y1 * w2 +*/ z1 * x2 - x1 * z2;
      var z3 = w1 * z2 + /*z1 * w2 +*/ x1 * y2 - y1 * x2;

      var w4 = w3 * w1 + x3 * x1 + y3 * y1 + z3 * z1;
      var x4 = x3 * w1 - w3 * x1 - y3 * z1 + z3 * y1;
      var y4 = y3 * w1 - w3 * y1 - z3 * x1 + x3 * z1;
      var z4 = z3 * w1 - w3 * z1 - x3 * y1 + y3 * x1;

      return [x4, y4, z4];
    },
    /**
     * Replaces the quaternion by a rotation given by axis and angle
     *
     * @see http://www.euclideanspace.com/maths/geometry/rotations/conversions/angleToQuaternion/index.htm
     * @param {Array} axis The axis around which to rotate
     * @param {number} angle The angle in radians
     * @returns {Quaternion}
     */
    'setFromAxisAngle': function(axis, angle) {

      // Q = [cos(angle / 2), v * sin(angle / 2)]

      var halfAngle = angle * 0.5;

      var a = axis[0];
      var b = axis[1];
      var c = axis[2];

      var sin = Math.sin(halfAngle);
      var cos = Math.cos(halfAngle);

      var sin_norm = sin / Math.hypot(a, b, c);

      this['w'] = cos;
      this['x'] = a * sin_norm;
      this['y'] = b * sin_norm;
      this['z'] = c * sin_norm;

      return this;
    },
    /**
     * Replaces the quaternion by a rotation given by Euler angles
     *
     * @param {number} alpha
     * @param {number} beta
     * @param {number} gamma
     * @returns {Quaternion}
     */
    'setFromEuler': function(alpha, beta, gamma, order) {

      var _x = beta;
      var _y = gamma;
      var _z = alpha;

      var cX = Math.cos(_x * 0.5);
      var cY = Math.cos(_y * 0.5);
      var cZ = Math.cos(_z * 0.5);

      var sX = Math.sin(_x * 0.5);
      var sY = Math.sin(_y * 0.5);
      var sZ = Math.sin(_z * 0.5);

      //
      // ZXY quaternion construction.
      //

      this['w'] = cX * cY * cZ - sX * sY * sZ;
      this['x'] = sX * cY * cZ - cX * sY * sZ;
      this['y'] = cX * sY * cZ + sX * cY * sZ;
      this['z'] = cX * cY * sZ + sX * sY * cZ;

      return this;
    },

    'ln' : function () {
      if (this.x === 0 && this.y === 0 && this.z === 0) {
        return new Quaternion(Math.log(x),0,0,0);
      } else {

        var n = Math.hypot(this['w'], this['x'], this['y'], this['z']);
        var vn = Math.pow(this['x'] * this['x']
          +this['y'] * this['y'] + this['z'] * this['z'], -0.5);
        var acos = Math.acos(this['w'] / n);

        var w = Math.log(n);
        var x = acos * this['x'] * vn;
        var y = acos * this['y'] * vn;
        var z = acos * this['z'] * vn;

        return new Quaternion(w, x, y, z);
      }
    },

    'exp' : function () {
      if (this.x === 0 && this.y === 0 && this.z === 0) {
        return new Quaternion(Math.exp(this['w']),0, 0, 0);
      } else {
        var eToW = Math.exp(this['w']);
        var vn = Math.hypot(this['x'], this['y'], this['z']);
        var rvn = 1/vn
        var sinvn = Math.sin(vn);

        var w = eToW * Math.cos(vn);
        var x = eToW * sinvn * this['x'] * rvn;
        var y = eToW * sinvn * this['y'] * rvn;
        var z = eToW * sinvn * this['z'] * rvn;

        return new Quaternion(w,x,y,z);
      }
    },

    'pow' : function (w,x,y,z) {
      if (w === 0 && x === 0 && y === 0 && z === 0) {
        return new Quaternion(1,0,0,0);
      }
      else if (Number.isInteger(w) && x === 0 && y === 0 && z === 0) {
        // if the imput quaternion is a integer and the rest of the quaternion components are 0
        var q = this;
        if (w < 0 ) {
          q = q.inverse()
          w = -w;
        }

        var w1 = q['w'];
        var x1 = q['x'];
        var y1 = q['y'];
        var z1 = q['z'];

        for (var i = 1; i < w; i++) {
          var w2 = q['w'];
          var x2 = q['x'];
          var y2 = q['y'];
          var z2 = q['z'];

          q.w = w2 * w1 - x2 * x1 - y2 * y1 - z2 * z1;
          q.x = w2 * x1 + x2 * w1 + y2 * z1 - z2 * y1;
          q.y = w2 * y1 + y2 * w1 + z2 * x1 - x2 * z1;
          q.z = w2 * z1 + z2 * w1 + x2 * y1 - y2 * x1;
        }
        return new Quaternion(q.w, q.x, q.y, q.z);
      } else if (x === 0 && y === 0 && z === 0) {// if the exponent is real but non integer
        /*
         * all quaternions can be writen as:
         * q = w +  xi + yj + zk = a + v
         * q = |q|*(cos(arg) + n*sin(arg))
         *
         * where:
         * n = v/|v|
         * w = |q|*cos(arg)
         *
         * following from this similar to D'Moivres Theorem:
         *
         * q^p = (|q|^p)*(cos(arg*p) + n*sin(arg*p))
         *
         * where p is any real number
         *
         * the modulus is raised the power, the argument is multiplied by the power
         * and the unit vector n remains the same
         */
         var norm = Math.hypot(this['w'], this['x'], this['y'], this['z']);
         var arg = Math.acos(this['w']/norm);
         var qtop = Math.pow(norm,w);
         var vn = Math.pow(this['x'] * this['x'] + this['y'] * this['y'] + this['z'] * this['z'], -0.5);
         var sinparg = Math.sin(arg * w)

         var w1 = qtop * Math.cos(arg * w);
         var x1 = qtop * this['x'] * vn * sinparg;
         var y1 = qtop * this['y'] * vn * sinparg;
         var z1 = qtop * this['z'] * vn * sinparg;

         return new Quaternion(w1, x1, y1, z1);
      } else { // if the exponent is non real
      /*
       * formula for evaluating a Quaternion to a Quaternion power
       *
       * q^p === e^ln(q^p) === e^[ln(q)*p]
       */
        // finds log base e of this
        var n = Math.hypot(this['w'], this['x'], this['y'], this['z']);
        var vn = Math.pow(this['x'] * this['x']
          +this['y'] * this['y'] + this['z'] * this['z'], -0.5);
        var acos = Math.acos(this['w'] / n);

        var w1 = Math.log(n);
        var x1 = acos * this['x'] * vn;
        var y1 = acos * this['y'] * vn;
        var z1 = acos * this['z'] * vn;

        // multiplies log base e of this by the exponent
        var w2 = w1 * w - x1 * x - y1 * y - z1 * z;
        var x2 = w1 * x + x1 * w + y1 * z - z1 * y;
        var y2 = w1 * y + y1 * w + z1 * x - x1 * z;
        var z2 = w1 * z + z1 * w + x1 * y - y1 * x;

        // raises all of the above to power e
        var eToW = Math.exp(w2);
        var vn = Math.hypot(x2, y2, z2);
        var rvn = 1/vn
        var sinvn = Math.sin(vn);

        var w3 = eToW * Math.cos(vn);
        var x3 = eToW * sinvn * x2 * rvn;
        var y3 = eToW * sinvn * y2 * rvn;
        var z3 = eToW * sinvn * z2 * rvn;

        return new Quaternion(w3, x3, y3, z3);
      }
    },

    'getArgument' : function () {
      return Math.acos(this.x * Math.pow(this['w'] * this['w'] + this['x'] * this['x']
          +this['y'] * this['y'] + this['z'] * this['z'],-0.5));
    },

    'toModArgUnitForm' : function () {
      var vn = Math.pow(this['x'] * this['x'] + this['y'] * this['y'] +
       this['z'] * this['z'], -0.5);
      var n = Math.hypot(this['w'], this['x'], this['y'], this['z']);
      var arg = Math.acos(this['w']/n);
      var unit = new Quaternion(0, this.x * vn, this.y * vn, this.z * vn);
      return {mod: n, arg: arg, unit: unit};
    },

    'round': function (places) {
      if (typeof places === 'undefined' || places < 0 || Math.round(places) !== places){
        places = 1;
      } else {
        places = Math.pow(10,places);
      }

      return new Quaternion(
        Math.round(this.w*places)/places, Math.round(this.x*places)/places,
        Math.round(this.y*places)/places, Math.round(this.z*places)/places
      );
    },

    'ceil': function (places) {
      if (typeof places === 'undefined' || places < 0 || Math.round(places) !== places){
        places = 1;
      } else {
        places = Math.pow(10,places);
      }

      return new Quaternion(
        Math.ceil(this.w*places)/places, Math.ceil(this.x*places)/places,
        Math.ceil(this.y*places)/places, Math.ceil(this.z*places)/places
      );
    },

    'floor': function(places) {
      if (typeof places === 'undefined' || places < 0 || Math.round(places) !== places){
        places = 1;
      } else {
        places = Math.pow(10,places);
      }

      return new Quaternion(
        Math.floor(this.w*places)/places, Math.floor(this.x*places)/places,
        Math.floor(this.y*places)/places, Math.floor(this.z*places)/places
      );
    },

    'fix': function () {
      return new Quaternion(
        (this.w > 0) ? Math.floor(this.w): Math.ceil(this.w),
        (this.x > 0) ? Math.floor(this.x): Math.ceil(this.x),
        (this.y > 0) ? Math.floor(this.y): Math.ceil(this.y),
        (this.z > 0) ? Math.floor(this.z): Math.ceil(this.z));

    },

    'dotDivide': function (a,b,c,d) {
      return new Quaternion(this.w/a, this.x/b, this.y/c, this.z/d);
    },

    'dotPow': function (a) {
      if ( (a % 1 !== 0 && a < 0) &&
        (this.w < 0 || this.x < 0 || this.y < 0 || this.z < 0)) {
        throw "invalid power";
      }
      return new Quaternion(Math.pow(this.w,a), Math.pow(this.x,a),
        Math.pow(this.y,a), Math.pow(this.z,a));
    },

    'toMatrix4x4' : function () {
      return [
        [this.w, -this.x, -this.y, -this.z],
        [this.x, this.w, -this.z, this.y],
        [this.y, this.z, this.w, -this.x],
        [this.z, -this.y, this.x, this.w]
      ];
    }


  };

  Quaternion['ZERO'] = new Quaternion(0, 0, 0, 0); // This is the additive identity Quaternion
  Quaternion['ONE'] = new Quaternion(1, 0, 0, 0); // This is the multiplicative identity Quaternion
  Quaternion['I'] = new Quaternion(0, 1, 0, 0);
  Quaternion['J'] = new Quaternion(0, 0, 1, 0);
  Quaternion['K'] = new Quaternion(0, 0, 0, 1);

  if (typeof define === 'function' && define['amd']) {
    define([], function() {
      return Quaternion;
    });
  } else if (typeof exports === 'object') {
    module['exports'] = Quaternion;
  } else {
    root['Quaternion'] = Quaternion;
  }

})(this);
