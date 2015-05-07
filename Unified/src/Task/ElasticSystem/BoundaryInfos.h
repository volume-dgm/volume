#pragma once
#include "../../Maths/Spaces.h"

// boundaries
struct BoundaryConditions
{
  enum Types
  {
    Absorb = 0, Free = 1, Fixed, Symmetry, AntiSymmetry
  };
};

struct BoundaryDescription
{
  BoundaryDescription(): infoIndex(-1),
    reflectionCoeff(1.0)
  {}
  BoundaryConditions::Types type;
  int infoIndex;
  double reflectionCoeff;
};

template <typename Space>
struct BoundaryInfoFunctor
{
  SPACE_TYPEDEFS

  virtual void operator()(const Vector& globalPoint,
                          const Vector& externalNormal,
                          const Scalar currTime,
                          Scalar* values) = 0;
};

template <typename Space>
struct FreeBoundaryInfo
{
  SPACE_TYPEDEFS
  FreeBoundaryInfo(): 
    boundaryFunctor(0), 
    dynamicContactInteractionType(IndexType(-1))
  {
  }
  BoundaryInfoFunctor<Space>* boundaryFunctor;
  IndexType dynamicContactInteractionType;
};

template <typename Space>
struct FixedBoundaryInfo
{
  SPACE_TYPEDEFS

  FixedBoundaryInfo(): boundaryFunctor(0)
  {
  }
  BoundaryInfoFunctor<Space>* boundaryFunctor;
};
