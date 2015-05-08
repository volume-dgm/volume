#pragma once
#include "VectorFunctors.h"
#include "SourceFunctors.h"

#include "../../Maths/Spaces.h"
#include "ElasticSystemCommon.h"

#include "Eigen/Dense"

template <typename Space>
struct ElasticSystem;

template <>
struct ElasticSystem<Space2>: public ElasticSystemCommon<Space2>
{
  SPACE2_TYPEDEFS
  typedef Space2                                        Space;
  typedef ElasticSystemCommon<Space2>::ValueTypeCommon  ValueTypeCommon;
  typedef ElasticSystemCommon<Space2>::MediumParameters MediumParameters;

  struct ValueType: public ValueTypeCommon
  {
    using ValueTypeCommon::SetTension;
    using ValueTypeCommon::values;
    void SetTension(Scalar xx, Scalar yy, Scalar xy);
    Scalar GetPressure() const;
    Scalar GetDeviatorSquare() const;
    Vector GetForce(const Vector& normal) const;
    Tensor GetTension() const;
  };

  struct ElasticSpace: public Space2
  {
    SPACE2_TYPEDEFS
    typedef          Space2               SpaceType;
    typedef          ValueType            Elastic;
    typedef          MediumParameters     MediumParametersType;
  };

  void BuildEdgeTransformMatrix(Vector edgeVertices[2], Eigen::Matrix<Scalar, dimsCount, dimsCount>& transformMatrix);
  void BuildEdgeTransformMatrixInv(Vector edgeVertices[2], Eigen::Matrix<Scalar, dimsCount, dimsCount>& transformMatrixInv);
};

template <>
struct ElasticSystem<Space3>: public ElasticSystemCommon<Space3>
{
  SPACE3_TYPEDEFS
  typedef Space3                                        Space;
  using ElasticSystemCommon<Space3>::ValueTypeCommon;
  using ElasticSystemCommon<Space>::MediumParameters;

  struct ValueType: public ValueTypeCommon
  {
    using ValueTypeCommon::SetTension;
    using ValueTypeCommon::values;
    void SetTension(Scalar xx, Scalar yy, Scalar zz, Scalar xy, Scalar yz, Scalar xz);
    Tensor GetTension() const;

    Scalar GetZZ() const;
    Scalar GetYZ() const;
    Scalar GetXZ() const;
    Scalar GetPressure() const;
    Scalar GetDeviatorSquare() const;
    Vector GetForce(const Vector& normal) const;
  };

  struct ElasticSpace: public Space3
  {
    SPACE3_TYPEDEFS
    typedef          Space3               SpaceType;
    typedef          ValueType            Elastic;
    typedef          MediumParameters     MediumParametersType;
  };

  void BuildFaceTransformMatrix(Vector faceVertices[3], Eigen::Matrix<Scalar, dimsCount, dimsCount>& transformMatrix);
  void BuildFaceTransformMatrixInv(Vector faceVertices[3], Eigen::Matrix<Scalar, dimsCount, dimsCount>& transformMatrixInv);
  void BuildZMatrix(const MediumParameters& mediumParameters, Eigen::Matrix<Scalar, dimsCount, dimsCount>& zMatrix);
};

#include "ElasticSystem.inl"
