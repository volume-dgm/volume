#pragma once

#include "../../../Maths/Spaces.h"
#include "../../../Maths/MatrixMaths.h"
#include "../../../Maths/QuadraturePrecomputer.h"
#include <assert.h>

template <typename Space, typename Decomposer, typename Function, typename FuncArgType>
struct GetterCommon
{
  GetterCommon(Decomposer* decomposer, Function& func) : decomposer(decomposer), func(func)
  {}
  Decomposer* decomposer;
  Function& func;
};

template <typename Space, typename Decomposer, typename Function, typename FuncArgType>
struct Getter;

template <typename Space, typename Decomposer, typename Function>
struct Getter<Space, Decomposer, Function, typename Space::IndexType> : public GetterCommon < Space, Decomposer, Function, typename Space::IndexType>
{
  SPACE_TYPEDEFS
  Getter(Decomposer* decomposer, Function& func) : GetterCommon(decomposer, func)
  {}
  void operator()(IndexType basisPointIndex, Scalar* values)
  {
    func(basisPointIndex, values);
  }
};

template <typename Space, typename Decomposer, typename Function>
struct Getter< Space, Decomposer, Function, typename Space::Vector> : public GetterCommon < Space, Decomposer, Function, typename Space::Vector>
{
  SPACE_TYPEDEFS
  Getter(Decomposer* decomposer, Function& func) : GetterCommon(decomposer, func)
  {}
  void operator()(IndexType basisPointIndex, Scalar* values)
  {
    func(decomposer->GetBasisPoints()[basisPointIndex], values);
  }
};

template<typename Space, typename FunctionSpaceT>
struct QuadratureDecomposer : public FunctionSpaceT
{
  SPACE_TYPEDEFS

  typedef FunctionSpaceT FunctionSpace;
  typedef QuadratureDecomposer<Space, FunctionSpace> DecomposerType;

  using FunctionSpace::functionsCount;
  using FunctionSpace::order;

  QuadratureDecomposer()
  {
    //building quadrature points and weights
    QuadraturePrecomputer::BuildQuadrature<Space>(2 * order, weights, points);

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

    basisFunctionValues.resize(functionsCount * points.size());

    for (IndexType pointIndex = 0; pointIndex < points.size(); pointIndex++)
    {
      for (IndexType functionIndex = 0; functionIndex < functionsCount; functionIndex++)
      {
        basisFunctionValues[pointIndex * functionsCount + functionIndex] = this->GetBasisFunctionValue(points[pointIndex], functionIndex);
      }
    }

    MatrixInverse<Scalar, IndexType>(cellVolumeIntegrals, cellVolumeIntegralsInv, functionsCount);
  }

  template<typename Function, int ValuesCount>
  void DecomposePrecomputed(Function& func, Scalar* coords)
  {
    DecomposeBase<Function, ValuesCount, IndexType>(func, coords);
  }

  template<typename Function, int ValuesCount>
  void Decompose(Function& func, Scalar* coords)
  {
    DecomposeBase<Function, ValuesCount, Vector>(func, coords);
  }

  const std::vector<Vector>& GetBasisPoints() const
  {
    return points;
  }

private:
  // for quadrature integration
  std::vector<Scalar> weights;
  std::vector<Vector> points;

  std::vector<Scalar> basisFunctionValues;

  Scalar cellVolumeIntegrals[FunctionSpaceT::functionsCount * FunctionSpaceT::functionsCount];
  Scalar cellVolumeIntegralsInv[FunctionSpaceT::functionsCount * FunctionSpaceT::functionsCount];

  template<typename Function, int ValuesCount, typename FuncArgType>
  void DecomposeBase(Function& func, Scalar* coords)
  {
    Scalar correlations[functionsCount * ValuesCount]; // correlations == inner products
    std::fill_n(correlations, ValuesCount * functionsCount, 0);

    Scalar values[ValuesCount];
    Getter<Space, DecomposerType, Function, FuncArgType> getter(this, func);

    for (IndexType pointIndex = 0; pointIndex < points.size(); pointIndex++)
    {
      getter(pointIndex, values);

      for (IndexType functionIndex = 0; functionIndex < functionsCount; functionIndex++)
      {
        Scalar basisFunctionValue = basisFunctionValues[pointIndex * functionsCount + functionIndex];
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
};

template<typename Space, typename PolynomialSpaceT>
struct NodalDecomposer : public PolynomialSpaceT
{
  SPACE_TYPEDEFS

  typedef PolynomialSpaceT PolynomialSpace;
  typedef NodalDecomposer<Space, PolynomialSpace> DecomposerType;

  using PolynomialSpace::functionsCount;
  using PolynomialSpace::order;

  NodalDecomposer()
  {
    basisPoints.reserve(functionsCount);
    for (IndexType functionIndex = 0; functionIndex < functionsCount; functionIndex++)
    {
      basisPoints.push_back(this->GetBasisPoint(functionIndex));
    }
  }

  template<typename Function, int ValuesCount>
  void Decompose(Function& func, Scalar* coords)
  {
    DecomposeBase<Function, ValuesCount, Vector>(func, coords);
  }

  template<typename Function, int ValuesCount>
  void DecomposePrecomputed(Function& func, Scalar* coords)
  {
    DecomposeBase<Function, ValuesCount, IndexType>(func, coords);
  }

  const std::vector<Vector>& GetBasisPoints() const
  {
    return basisPoints;
  }

private:
  std::vector<Vector> basisPoints;

  template<typename Function, int ValuesCount, typename FuncArgType>
  void DecomposeBase(Function& func, Scalar* coords)
  {
    Scalar values[ValuesCount];
    Getter<Space, DecomposerType, Function, FuncArgType> getter(this, func);

    for (IndexType functionIndex = 0; functionIndex < functionsCount; functionIndex++)
    {
      getter(functionIndex, values);

      for (IndexType valueIndex = 0; valueIndex < ValuesCount; valueIndex++)
      {
        coords[valueIndex * functionsCount + functionIndex] = values[valueIndex];
      }
    }
  }
};
