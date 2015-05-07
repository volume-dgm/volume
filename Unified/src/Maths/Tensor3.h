#pragma once
#include "Vector3.h"
#include "Matrix3x3.h"

template<typename T>
class Tensor3
{
public:
  Tensor3()
  {
    xx = xy = xz = yy = yz = zz = 0;
  }
  Tensor3(T xx, T xy, T xz, T yy, T yz, T zz)
  {
    this->xx = xx;
    this->xy = xy;
    this->xz = xz;
    this->yy = yy;
    this->yz = yz;
    this->zz = zz;
  }
  Tensor3(T diag)
  {
    xx = yy = zz = diag;
    xy = xz = yz = 0;
  }
  Vector3<T> operator ()(const Vector3<T> &v) const
  {
    return Vector3<T>(
      xx * v.x + xy * v.y + xz * v.z,
      xy * v.x + yy * v.y + yz * v.z,
      xz * v.x + yz * v.y + zz * v.z);
  }
    
  T GetTrace()
  {
    return xx + yy + zz;
  }

  T& operator[](const int i)
  {
    return *(&(xx) + i);
  }

  const T& operator[](const int i) const
  {
    return *(&(xx)+i);
  }

  Tensor3& operator +=(const T &v)
  {
    xx += v;
    yy += v;
    zz += v;
    return *this;
  }

  Tensor3& operator +=(const Tensor3& t)
  {
    xx += t.xx;
    xy += t.xy;
    xz += t.xz;
    yy += t.yy;
    yz += t.yz;
    zz += t.zz;
    return *this;
  }

  T xx, xy, xz, yy, yz, zz;
};

template <typename T>
Tensor3<T> TensorProduct(const Vector3<T> &v)
{
  return Tensor3<T>(
    v.x * v.x, //xx
    v.x * v.y, //xy
    v.x * v.z, //xz
    v.y * v.y, //yy
    v.y * v.z, //yz
    v.z * v.z  //zz
  );
}

template <typename T>
Tensor3<T> SymmetricTensorProduct(const Vector3<T> &v0, const Vector3<T> &v1) {
  return Tensor3<T>(
    T(2.0) * v0.x * v1.x, //xx
    v0.x * v1.y + v0.y * v1.x, //xy
    v0.x * v1.z + v0.z * v1.x, //xz
    T(2.0) * v0.y * v1.y, //yy
    v0.y * v1.z + v0.z * v1.y, //yz
    T(2.0) * v0.z * v1.z  //zz
  );
}


template<typename T>
inline const Tensor3<T> operator +(const Tensor3<T> &t0, const Tensor3<T> &t1)
{
  return Tensor3<T>(
    t0.xx + t1.xx,
    t0.xy + t1.xy,
    t0.xz + t1.xz,
    t0.yy + t1.yy,
    t0.yz + t1.yz,
    t0.zz + t1.zz);
}

template<typename T>
inline const Tensor3<T> operator -(const Tensor3<T> &t0, const Tensor3<T> &t1)
{
  return Tensor3<T>(
    t0.xx - t1.xx,
    t0.xy - t1.xy,
    t0.xz - t1.xz,
    t0.yy - t1.yy,
    t0.yz - t1.yz,
    t0.zz - t1.zz);
}

template<typename T>
inline const Tensor3<T> operator *(const Tensor3<T> &t, const T &s)
{
  return Tensor3<T>(
    t.xx * s,
    t.xy * s,
    t.xz * s,
    t.yy * s,
    t.yz * s,
    t.zz * s);
}

template <class T>
inline T DoubleConvolution(const Tensor3<T> &v0, const Tensor3<T> &v1)
{
  return 
    v0.xx * v1.xx + 
    v0.yy * v1.yy + 
    v0.zz * v1.zz + 
    T(2.0) * v0.xy * v1.xy +
    T(2.0) * v0.xz * v1.xz +
    T(2.0) * v0.yz * v1.yz;
}

typedef Tensor3<float>	Tensor3f;
typedef Tensor3<double> Tensor3d;
