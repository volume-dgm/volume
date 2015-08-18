#pragma once

#include "../../../Maths/Spaces.h"
#include "../../../Maths/MatrixMaths.h"
#include "../../../Maths/QuadraturePrecomputer.h"
#include <assert.h>

template<typename Space, typename FunctionSpaceT>
struct QuadratureDecomposer : public FunctionSpaceT
{
  SPACE_TYPEDEFS

  typedef FunctionSpaceT FunctionSpace;

  using FunctionSpace::functionsCount;
  using FunctionSpace::order;

  QuadratureDecomposer()
  {
    //building quadrature points and weights
    QuadraturePrecomputer::BuildQuadrature<Space>(order, weights, points);

    //building function integrals
    std::fill_n(cellVolumeIntegrals, functionsCount * functionsCount, 0);

    for (IndexType functionIndex0 = 0; functionIndex0 < functionsCount; functionIndex0++)
    {

      for (IndexType functionIndex1 = 0; functionIndex1 < functionsCount; functionIndex1++)
      {
        for (IndexType pointIndex = 0; pointIndex < points.size(); pointIndex++)
        {
          Scalar functionValue0 = this->GetBasisFunctionValue(points[pointIndex], functionIndex0);
          Scalar functionValue1 = this->GetBasisFunctionValue(points[pointIndex], functionIndex1);

          cellVolumeIntegrals[functionIndex0 * functionsCount + functionIndex1] += functionValue0 * functionValue1 * weights[pointIndex];
        }
      }
    }

    MatrixInverse<Scalar, IndexType>(cellVolumeIntegrals, cellVolumeIntegralsInv, functionsCount);
  }

  template<typename Function, int ValuesCount>
  void Decompose(Function& func, Scalar* coords)
  {
    Scalar correlations[functionsCount * ValuesCount]; // correlations == inner products
    std::fill_n(correlations, ValuesCount * functionsCount, 0);

    for (IndexType pointIndex = 0; pointIndex < points.size(); pointIndex++)
    {
      Scalar values[ValuesCount];
      func(points[pointIndex], values);

      for (IndexType functionIndex = 0; functionIndex < functionsCount; functionIndex++)
      {
        Scalar basisFunctionValue = this->GetBasisFunctionValue(points[pointIndex], functionIndex);
        for (IndexType valueIndex = 0; valueIndex < ValuesCount; valueIndex++)
        {
          correlations[valueIndex * functionsCount + functionIndex] += basisFunctionValue * values[valueIndex] * weights[pointIndex];
        }
      }
    }

    std::fill_n(coords, ValuesCount * functionsCount, 0);
    for (IndexType valueIndex = 0; valueIndex < ValuesCount; valueIndex++)
      for (IndexType i = 0; i < functionsCount; i++)
        for (IndexType j = 0; j < functionsCount; j++)
          coords[valueIndex * functionsCount + i] +=
          correlations[valueIndex * functionsCount + j] * cellVolumeIntegralsInv[j * functionsCount + i];
  }

  /*Scalar GetBasisFunctionValue(Vector point, IndexType functionIndex) const
  {
    return space.GetBasisFunctionValue(point, functionIndex);
  }*/

  // for quadrature integration
  std::vector<Scalar> weights;
  std::vector<Vector> points;

protected:
  Scalar cellVolumeIntegrals[FunctionSpaceT::functionsCount * FunctionSpaceT::functionsCount];
  Scalar cellVolumeIntegralsInv[FunctionSpaceT::functionsCount * FunctionSpaceT::functionsCount];
};

template<typename Space, typename PolynomialSpaceT>
struct NodalDecomposer : public PolynomialSpaceT
{
  SPACE_TYPEDEFS

  typedef PolynomialSpaceT PolynomialSpace;

  using PolynomialSpace::functionsCount;
  using PolynomialSpace::order;

  NodalDecomposer()
  {
  }

  template<typename Function, int ValuesCount>
  void Decompose(Function& func, Scalar* coords)
  {
    for(IndexType functionIndex = 0; functionIndex < functionsCount; functionIndex++)
    {
      Vector basisPoint = this->GetBasisPoint(functionIndex);
      Scalar values[ValuesCount];
      func(basisPoint, values);

      for(IndexType valueIndex = 0; valueIndex < ValuesCount; valueIndex++)
      {
        coords[valueIndex * functionsCount + functionIndex] = values[valueIndex];
      }
    }
  }
};
