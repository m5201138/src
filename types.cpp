#include <cmath>
#include <iostream>

#include "types.h"


Vec3 Vec3::cross(const Vec3& rhs) const {
  Vec3 tmp;
  tmp.x = y * rhs.z - z * rhs.y;
  tmp.y = z * rhs.x - x * rhs.z;
  tmp.z = x * rhs.y - y * rhs.x;
  return tmp;
}

double Vec3::norm() const {
  return std::sqrt(x*x + y*y + z*z);
}

double Vec3::dot(const Vec3& rhs) const {
  return x*rhs.x + y*rhs.y + z*rhs.z;
}

Vec3 Vec3::normalized() const {
  double n = norm();
  Vec3 tmp; 
  if (std::fabs(n) > 1e-10) {
    tmp.x = x / n;
    tmp.y = y / n;
    tmp.z = z / n;
  }
  return tmp;
}

Vec3 Vec3::Zero() {
  Vec3 tmp;
  tmp.x = 0.0;
  tmp.y = 0.0;
  tmp.z = 0.0;
  return tmp;
}

Vec3& Vec3::operator+=(const Vec3& rhs) {
  x += rhs.x;
  y += rhs.y;
  z += rhs.z;
  return *this;
}

Vec3& Vec3::operator/=(double scale) {
  if (scale==0.0) {
    std::cerr << "scale is 0 in Vec3::operator/=" << std::endl;
  }
  x /= scale;
  y /= scale;
  z /= scale;
  return *this;
}

void Vec3::normalize() {
  double n = norm();
  if (std::fabs(n) > 1e-10) {
    x = x / n;
    y = y / n;
    z = z / n;
  }
}

Vec3& Vec3::operator-=(const Vec3& rhs) {
  x -= rhs.x;
  y -= rhs.y;
  z -= rhs.z;
  return *this;
}

Vec3 operator*(double scale, const Vec3& rhs) {
  Vec3 tmp;
  tmp.x = rhs.x * scale;
  tmp.y = rhs.y * scale;
  tmp.z = rhs.z * scale;
  return tmp;
}

Vec3 operator-(const Vec3& lhs, const Vec3& rhs) {
  Vec3 tmp;
  tmp.x = lhs.x - rhs.x;
  tmp.y = lhs.y - rhs.y;
  tmp.z = lhs.z - rhs.z;
  return tmp;
}

