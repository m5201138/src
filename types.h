#ifndef TYPES_H
#define TYPES_H

struct Vec3 {
  double x, y, z;

  Vec3() : x(0.0), y(0.0), z(0.0) {}
  Vec3(double x, double y, double z) : x(x), y(y), z(z) {}
  Vec3(const Vec3& other) : x(other.x), y(other.y), z(other.z) {}
  
  Vec3& operator=(const Vec3& other) {
    x = other.x;
    y = other.y;
    z = other.z;
    return *this;
  }

  Vec3 cross(const Vec3& rhs) const;
  double norm() const;
  double dot(const Vec3& rhs) const;
  Vec3 normalized() const;
  static Vec3 Zero();

  Vec3& operator+=(const Vec3& rhs);
  Vec3& operator/=(double scale);
  Vec3& operator-=(const Vec3& rhs);

  void normalize();

};

Vec3 operator-(const Vec3& lhs, const Vec3& rhs);
Vec3 operator*(double scale, const Vec3& rhs);

#endif // TYPES_H
