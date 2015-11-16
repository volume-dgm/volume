#pragma once

#include "../../Maths/Spaces.h"
#include "BoundaryInfos.h"
#include "ContactInfos.h"
#include "PointSources.h"

#include "../VolumeMethod/FunctionGetters/BoundaryFunctionGetter.h"
#include <vector>

#include "Eigen/Dense"

template <typename Space>
struct ElasticSystemBase;

template <>
struct ElasticSystemBase<Space2>
{
  static const int dimsCount = 5;
};

template <>
struct ElasticSystemBase<Space3>
{
  static const int dimsCount = 9;
};

template <typename Space>
struct ElasticSystemCommon: public ElasticSystemBase<Space>
{
  SPACE_TYPEDEFS
  typedef SourceFunctor<Space> SourceFunctorT;
  using ElasticSystemBase<Space>::dimsCount;

  struct ValueTypeCommon
  {
    Scalar values[ElasticSystemBase<Space>::dimsCount];

    ValueTypeCommon& operator=(const ValueTypeCommon& other)
    {
      std::copy(other.values, other.values + ElasticSystemBase<Space>::dimsCount, values);
      return *this;
    }

    void SetZeroValues();
    Vector GetVelocity() const;
    void SetVelocity(const Vector& velocity);

    Tensor GetTension() const;
    void SetTension(const Tensor& tension);

    void MakeDimension(Scalar tensionDimensionlessMult, Scalar velocityDimensionlessMult);

    virtual Vector GetForce(const Vector& normal) const = 0;

    Scalar GetXX() const;
    Scalar GetYY() const;
    Scalar GetXY() const;
  };

  struct MediumParameters
  {
    SPACE_TYPEDEFS

    MediumParameters();
    MediumParameters(Scalar lambda, Scalar mju, Scalar invRho);

    void SetFromVelocities(Scalar pSpeed, Scalar sSpeed, Scalar rho);
    void SetFromLameParameters(Scalar lambda, Scalar mju, Scalar rho);
    void SetFromYoungModuleAndNju(Scalar youngModule, Scalar nju, Scalar rho);

    Scalar GetPSpeed() const;
    Scalar GetSSpeed() const;

    Scalar GetZp() const;
    Scalar GetZs() const;

    void MakeDimensionless(Scalar tensionDimensionlessMult, Scalar velocityDimensionlessMult);

    static const int ParamsCount = 3;
    union 
    {
      struct 
      {
        Scalar lambda;
        Scalar mju;
        Scalar invRho;
      };
      Scalar params[ParamsCount];
    };

    bool IsZero() const;

    Vector flowVelocity;

    bool destroyed;
    bool dimensionless;

    struct Plasticity
    {
      // plasticity: 
      // s:s < 2 * k^2 
      // k = k0 + a * pressure
      Plasticity(): k0(std::numeric_limits<Scalar>::infinity()), a(0) /* without plasticity */, brittle(false),
        maxPlasticWork(std::numeric_limits<Scalar>::infinity())
      {}
      Scalar k0;
      Scalar a;
      bool brittle;
      Scalar maxPlasticWork;
    };
    Plasticity plasticity;
    bool fixed;
  };
  
  void BuildXMatrix(const MediumParameters& mediumParameters, Eigen::Matrix<Scalar, dimsCount, dimsCount>& xMatrix);
  
  void BuildYMatrix(const MediumParameters& mediumParameters, Eigen::Matrix<Scalar, dimsCount, dimsCount>& yMatrix);
  
  void BuildRMatrix(const MediumParameters& interiorMediumParameters, 
                    const MediumParameters& exteriorMediumParameters, 
                    Eigen::Matrix<Scalar, dimsCount, dimsCount>& rMatrix);

  void BuildXnAuxMatrix(const MediumParameters& interiorMediumParameters, 
                        const MediumParameters& exteriorMediumParameters,
                        Eigen::Matrix<Scalar, dimsCount, dimsCount>& xnAuxMatrix);

  void BuildXnInteriorMatrix(const MediumParameters& interiorMediumParameters, 
                             const MediumParameters& exteriorMediumParameters,
                             const Vector& edgeNormal, Eigen::Matrix<Scalar, dimsCount, dimsCount>& xnInteriorMatrix);

  void BuildXnExteriorMatrix(const MediumParameters& interiorMediumParameters, 
                             const MediumParameters& exteriorMediumParameters,
                             const Vector& edgeNormal, Eigen::Matrix<Scalar, dimsCount, dimsCount>& xnExteriorMatrix);

  void BuildBoundaryMatrix(IndexType interactionType, Eigen::Matrix<Scalar, 1, dimsCount>& boundaryMatrix);

  void BuildContactMatrices(IndexType interactionType,
    Eigen::Matrix<Scalar, 1, dimsCount>& leftContactMatrix,
    Eigen::Matrix<Scalar, 1, dimsCount>& rightContactMatrix);
  
  Scalar GetMaxWaveSpeed(const MediumParameters& mediumParameters) const;

  bool IsProperContact(const ValueTypeCommon& value0, const ValueTypeCommon& value1, 
    const Vector& contactNormal);

  BoundaryInfoFunctor<Space>* GetBoundaryInfoFunctor(IndexType interactionType);

  void SetFreeBoundary(IndexType interactionType, Scalar tensionDimensionlessMult, Scalar reflectionCoeff,
    VectorFunctor<Space>* externalForceFunctor = 0, IndexType dynamicContactInteractionType = IndexType(-1));

  void SetFixedBoundary(IndexType interactionType, Scalar velocityDimensionlessMult, Scalar reflectionCoeff,
    VectorFunctor<Space>* externalVelocityFunctor = 0);

  void SetAbsorbBoundary(IndexType interactiontype);
  void SetSymmetryBoundary(IndexType interactionType);
  void SetAntiSymmetryBoundary(IndexType interactionType);
  BoundaryConditions::Types GetBoundaryType(IndexType interactionType) const;

  void SetGlueContact(IndexType interactionType);
  void SetGlueContact(IndexType interactionType, Scalar maxShearStress, Scalar maxLongitudinalStress,
    IndexType dynamicBoundaryInteractionType);
  void SetGlideContact(IndexType interactionType);
  void SetFrictionContact(IndexType interactionType, Scalar frictionCoeff, IndexType dynamicBoundaryInteractionType);
  void SetFrictionContact(IndexType interactionType);
  void GetFrictionContactInfo(IndexType interactionType, Scalar& frictionCoeff);

  ContactConditions::Types GetContactType(IndexType interactionType) const;

  IndexType GetContactDynamicBoundaryType(IndexType contactInteractionType ) const;
  IndexType GetBoundaryDynamicContactType(IndexType boundaryInteractionType) const;

  void GetContactCriticalInfo(IndexType interactionType, 
    Scalar& maxShearStress, Scalar& maxLongitudinalStress);

  SourceFunctorT* GetSourceFunctor() const;
  void SetSourceFunctor(SourceFunctorT* sourceFunctor);
  void SetPointSources(const std::vector<PointSource<Space>* >& pointSources);
  std::vector< PointSource<Space>* > pointSources;

  class ZeroBoundaryInfoFunctor: public BoundaryInfoFunctor<Space> //default functor
  {
  public:
    SPACE_TYPEDEFS
    void operator()(const Vector&, const Vector&, const Scalar, Scalar* values)
    {
      std::fill_n(values, dimsCount, Scalar(0.0));
    }
  };

  class FreeBoundaryInfoFunctor: public BoundaryInfoFunctor<Space>
  {
    SPACE_TYPEDEFS
  public:
    FreeBoundaryInfoFunctor(VectorFunctor<Space>* forceFunctor, Scalar tensionDimensionlessMult)
    {
      this->forceFunctor             = forceFunctor;
      this->tensionDimensionlessMult = tensionDimensionlessMult;
    }

    void operator()(const Vector& globalPoint, const Vector& externalNormal, const Scalar time, Scalar* values);

    void SetCurrentVelocity(const Vector& velocity)
    {
      forceFunctor->SetCurrentVelocity(velocity);
    }

  private:
    VectorFunctor<Space>* forceFunctor;
    Scalar tensionDimensionlessMult; 
  };

  class FixedBoundaryInfoFunctor: public BoundaryInfoFunctor<Space>
  {
    SPACE_TYPEDEFS
  public:
    FixedBoundaryInfoFunctor(VectorFunctor<Space>* velocityFunctor, Scalar velocityDimensionlessMult)
    {
      this->velocityFunctor           = velocityFunctor;
      this->velocityDimensionlessMult = velocityDimensionlessMult;
    }

    void operator()(const Vector& globalPoint, const Vector& externalNormal, const Scalar time, Scalar* values)
    {
      std::fill(values, values + dimsCount, 0);
      Vector externalVelocity = (*velocityFunctor)(globalPoint, externalNormal, time) / velocityDimensionlessMult;
      SetVelocity(externalVelocity, values);
    }

  private:
    void SetVelocity(const Vector& externalVelocity, Scalar* values);

    VectorFunctor<Space>* velocityFunctor;
    Scalar velocityDimensionlessMult;
  };

protected:
  void SetBoundaryDescription(IndexType interactionType, BoundaryDescription boundaryDescription);
  void SetContactDescription(IndexType interactionType, ContactDescription contactDescription);

  std::vector< FreeBoundaryInfo<Space>  > freeBoundaryInfos;
  std::vector< FixedBoundaryInfo<Space> > fixedBoundaryInfos;
  std::vector< BoundaryDescription >      boundaryDescriptions;

  std::vector< GlueContactInfo<Space> >     glueContactInfos;
  std::vector< FrictionContactInfo<Space> > frictionContactInfos;
  std::vector< ContactDescription >         contactDescriptions;
  SourceFunctorT* sourceFunctor; 

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#include "ElasticSystemCommon.inl"
