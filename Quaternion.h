#include <cmath>

#ifndef QUATERNION_H
#define QUATERNION_H


struct Quaternion {
  double x, y, z, w;

  Quaternion() : x(0.0), y(0.0), z(0.0), w(0.0) {}
  Quaternion(double w, double x, double y, double z) : x(x), y(y), z(z), w(w) {}

  double squaredNorm() const {
    return (x*x + y*y + z*z + w*w);
  }

  void normalize() {
    double n = std::sqrt(squaredNorm());
    if (n != 0.0) {
      x = x/n;
      y = y/n;
      z = z/n;
      w = w/n;
    }
  }

  Quaternion conjugate() const {
    Quaternion q;
    q.w = w;
    q.x = -x;
    q.y = -y;
    q.z = -z;
    return q;
  }

  Quaternion operator* (const Quaternion& b) const {
    Quaternion q;
    q.w = w * b.w - x * b.x - y * b.y - z * b.z;
    q.x = w * b.x + x * b.w + y * b.z - z * b.y;
    q.y = w * b.y + y * b.w + z * b.x - x * b.z;
    q.z = w * b.z + z * b.w + x * b.y - y * b.x;
    return q;
  }
};


#endif
