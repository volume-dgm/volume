#pragma once

#include "../../../Maths/Spaces.h"
#include "../../../Maths/MatrixMaths.h"
#include "../../../Maths/QuadraturePrecomputer.h"
#include <assert.h>

template<typename Space, typename PolynomialSpace>
struct Decomposer
{
  SPACE_TYPEDEFS

  const static IndexType functionsCount = PolynomialSpace::functionsCount;
  const static IndexType order = PolynomialSpace::order;
  PolynomialSpace space;

  virtual Scalar ComputeCellVolumeIntegral(IndexType functionIndex0, IndexType functionIndex1) = 0;

  Decomposer()
  {
    QuadraturePrecomputer::BuildQuadrature<Space>(order, weights, points);
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
        Scalar basisFunctionValue = space.GetBasisFunctionValue(points[pointIndex], functionIndex);
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

  Scalar GetBasisFunctionValue(Vector point, IndexType functionIndex) const
  {
    return space.GetBasisFunctionValue(point, functionIndex);
  }


  // for quadrature integration
  std::vector<Scalar> weights;
  std::vector<Vector> points;

protected:
  Scalar cellVolumeIntegrals[functionsCount * functionsCount];
  Scalar cellVolumeIntegralsInv[functionsCount * functionsCount];

  void ComputeVolumeIntegrals()
  {
    for (IndexType functionIndex0 = 0; functionIndex0 < functionsCount; functionIndex0++)
    {
      for (IndexType functionIndex1 = 0; functionIndex1 < functionsCount; functionIndex1++)
      {
        cellVolumeIntegrals[functionIndex0 * functionsCount + functionIndex1] =
          ComputeCellVolumeIntegral(functionIndex1, functionIndex0);
      }
    }
    MatrixInverse<Scalar, IndexType>(cellVolumeIntegrals, cellVolumeIntegralsInv, functionsCount);
  }
};

template<typename Space, typename PolynomialSpace>
struct DummyPrecomputer: public Decomposer<Space, PolynomialSpace>
{
  SPACE_TYPEDEFS

  using Decomposer<Space, PolynomialSpace>::functionsCount;
  using Decomposer<Space, PolynomialSpace>::order;
  using Decomposer<Space, PolynomialSpace>::space;

  DummyPrecomputer()
  {
    Decomposer<Space, PolynomialSpace>::ComputeVolumeIntegrals();
  }

  inline Scalar ComputeCellVolumeIntegral(IndexType functionIndex0, IndexType functionIndex1) // <Ô1, Ô2>
  {
    return space.ComputeCellVolumeIntegral(functionIndex0, functionIndex1);
  }

  inline Scalar GetBasisFunctionValue(const Vector& point, IndexType functionIndex) const
  {
    return space.GetBasisFunctionValue(point, functionIndex);
  }

  inline Scalar ComputeCellVolumeIntegral(IndexType functionIndex)
  {
    return space.ComputeCellVolumeIntegral(functionIndex);
  }

  inline Vector ComputeDerivativeVolumeIntegral(IndexType functionIndex0, IndexType functionIndex1) // {<dÔ1/dx, Ô2>, <dÔ1/dy, Ô2>}
  {
    return space.ComputeDerivativeVolumeIntegral(functionIndex0, functionIndex1);
  }

  inline Scalar ComputeOutgoingFlux(IndexType srcEdgeNumber, IndexType functionIndex0, IndexType functionIndex1)
  {
    return space.ComputeOutgoingFlux(srcEdgeNumber, functionIndex0, functionIndex1);
  }

  inline Scalar ComputeIncomingFlux(IndexType srcEdgeNumber, IndexType dstEdgeNumber, IndexType functionIndex0, IndexType functionIndex1)
  {
    return space.ComputeIncomingFlux(srcEdgeNumber, dstEdgeNumber, functionIndex0, functionIndex1);
  }

  inline Scalar ComputeEdgeFlux(IndexType edgeNumber, IndexType functionIndex)
  {
    return space.ComputeEdgeFlux(edgeNumber, functionIndex);
  }
};
