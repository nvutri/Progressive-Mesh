#ifndef _POINT_H_
#define _POINT_H_

#include <cmath>

/*!
* The Point Class
*/
class Point {
public:
    /*!
    *Constructor
    @param x x coordinate
    @param y y coordinate
    @param z z coordinate
    */
    Point(double x, double y, double z) {
        v[0] = x;
        v[1] = y;
        v[2] = z;
    }

    Point() {
        v[0] = v[1] = v[2] = 0;
    }

    /*! Default Constructor */
    ~Point() {
        ;
    }

    double &operator[](int i) {
        return v[i];
    }

    /*! Accessing the i-th coordinator */
    double norm() const {
        return sqrt(fabs(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]));
    }

    /*! Square root distance to the origin */
    double norm2() const {
        return fabs(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    }

    /*! Square distance to the origin */

    Point &operator+=(Point &p) {
        v[0] += p[0];
        v[1] += p[1];
        v[2] += p[2];
        return *this;
    }

    /*! Adding two point vectors (3-dimensional vectors) */
    Point &operator-=(Point &p) {
        v[0] -= p[0];
        v[1] -= p[1];
        v[2] -= p[2];
        return *this;
    }

    /*! Subtraction between two point vectors (3-dimensional vectors) */
    Point &operator*=(double s) {
        v[0] *= s;
        v[1] *= s;
        v[2] *= s;
        return *this;
    }

    /*! Scaling a (3-dimensional vectors) */
    Point &operator/=(double s) {
        v[0] /= s;
        v[1] /= s;
        v[2] /= s;
        return *this;
    }

    /*! Scale division (3-dimensional vectors) */

    double    operator*(Point &p) const {
        return v[0] * p[0] + v[1] * p[1] + v[2] * p[2];
    }

    Point    operator+(Point &p) const {
        return Point(v[0] + p[0], v[1] + p[1], v[2] + p[2]);
    }

    Point    operator-(Point &p) const {
        return Point(v[0] - p[0], v[1] - p[1], v[2] - p[2]);
    }

    Point    operator*(double s) const {
        return Point(v[0] * s, v[1] * s, v[2] * s);
    }

    Point    operator/(double s) const {
        return Point(v[0] / s, v[1] / s, v[2] / s);
    }

    Point    operator-() const {
        return Point(-v[0], -v[1], -v[2]);
    }

    Point    operator^(Point &p2) const {
        return Point(v[1] * p2[2] - v[2] * p2[1],
                v[2] * p2[0] - v[0] * p2[2],
                v[0] * p2[1] - v[1] * p2[0]);
    }


public:
    double v[3];
};

#endif