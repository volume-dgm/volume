#pragma once

#include "../../Maths/Spaces.h"
#include "BoundaryInfos.h"
#include "ContactInfos.h"
#include "PointSources.h"

#include "../VolumeMethod/FunctionGetters/BoundaryFunctionGetter.h"
#include <vector>

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
      for (IndexType valueIndex = 0; valueIndex < ElasticSystemBase<Space>::dimsCount; ++valueIndex)
      {
        values[valueIndex] = other.values[valueIndex];
      }
      return *this;
    }

    void SetZeroValues();
    Vector GetVelocity() const;
    void SetVelocity(const Vector& velocity);
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

    Scalar lambda;
    Scalar mju;
    Scalar invRho;
    Vector flowVelocity;

    bool destroyed;
    bool dimensionless;

    struct Plasticity
    {
      /* 
          plasticity: 
          s:s < 2 * k^2 
          k = k0 + a * pressure
      */
      Plasticity(): k0(std::numeric_limits<Scalar>::infinity()), a(0) // without plasticity
      {}
      Scalar k0;
      Scalar a;
    } plasticity;
  };

  void BuildXMatrix(const MediumParameters& mediumParameters, Scalar* xMatrix);
  void BuildYMatrix(const MediumParameters& mediumParameters, Scalar* yMatrix);

  void BuildRMatrix(const MediumParameters& interiorMediumParameters, 
                    const MediumParameters& exteriorMediumParameters, 
                    Scalar *rMatrix);

  void BuildXnAuxMatrix(const MediumParameters& interiorMediumParameters, 
                        const MediumParameters& exteriorMediumParameters,
                        Scalar* xnAuxMatrix);

  void BuildXnInteriorMatrix(const MediumParameters& interiorMediumParameters, 
                             const MediumParameters& exteriorMediumParameters,
                             const Vector& edgeNormal, Scalar* xnInteriorMatrix);

  void BuildXnExteriorMatrix(const MediumParameters& interiorMediumParameters, 
                             const MediumParameters& exteriorMediumParameters,
                             const Vector& edgeNormal, Scalar* xnExteriorMatrix);

  void BuildBoundaryMatrix(IndexType interactionType, Scalar *boundaryMatrix);

  void BuildContactMatrices(IndexType interactionType,
                            Scalar *leftContactMatrix,
                            Scalar *rightContactMatrix);

  void CorrectContact(Scalar* values, Vector contactNormal, IndexType dynamicInteractionType);

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

  IndexType GetContactDynamicBoundaryType(IndexType contactInteractionType );
  IndexType GetBoundaryDynamicContactType(IndexType boundaryInteractionType);

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

    void operator()(const Vector& globalPoint, const Vector& externalNormal, const Scalar time, Scalar* values)
    {
      std::fill_n(values, dimsCount, 0);
      Vector externalForce = (*forceFunctor)(globalPoint, externalNormal, time) / tensionDimensionlessMult;

      // TODO
      // 2d case
      //Vector force(externalForce * externalNormal, externalNormal ^ externalForce);

      // 3d case ?
      //Vector force(externalForce * externalNormal, externalNormal ^ externalForce);
      Vector force = Space::Vector::zeroVector();

      /*const Scalar frictionCoeff = Scalar(0.002);
      force.y *= -frictionCoeff * externalForce.GetNorm() * externalNormal;*/
      SetTension(externalNormal, force, values);
    }
  private:
    void SetTension(const Vector& externalNormal, 
      const Vector& force, Scalar* values);

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
      std::fill_n(values, dimsCount, 0);
      Vector externalVelocity = (*velocityFunctor)(globalPoint, externalNormal, time) / velocityDimensionlessMult;
      SetVelocity(externalVelocity, values);
    }

  private:
    void SetVelocity(const Vector& externalVelocity, Scalar* values);

    VectorFunctor<Space>* velocityFunctor;
    Scalar velocityDimensionlessMult;
  };

private:
  void SetBoundaryDescription(IndexType interactionType, BoundaryDescription boundaryDescription);
  void SetContactDescription(IndexType interactionType, ContactDescription contactDescription);

  std::vector< FreeBoundaryInfo<Space>  > freeBoundaryInfos;
  std::vector< FixedBoundaryInfo<Space> > fixedBoundaryInfos;
  std::vector< BoundaryDescription >      boundaryDescriptions;

  std::vector< GlueContactInfo<Space> >     glueContactInfos;
  std::vector< FrictionContactInfo<Space> > frictionContactInfos;
  std::vector< ContactDescription >         contactDescriptions;
  SourceFunctorT* sourceFunctor;
};

#include "ElasticSystemCommon.inl"
