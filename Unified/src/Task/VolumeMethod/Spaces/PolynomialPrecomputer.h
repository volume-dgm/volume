#pragma once

#include "../Cell.h"
#include "../../../Maths/Spaces.h"

#include "../../../../3rdparty/quadrature_integration/triangle_fekete_rule.hpp"
#include "../../../../3rdparty/quadrature_integration/tetrahedron_arbq_rule.hpp"
#include <assert.h>

template<typename Space, typename PolynomialSpace>
struct PolynomialPrecomputerCommon
{
  SPACE_TYPEDEFS

  const static IndexType functionsCount = PolynomialSpace::functionsCount;
  const static IndexType order = PolynomialSpace::order;
  PolynomialSpace space;

  virtual Scalar ComputeCellVolumeIntegral(IndexType functionIndex0, IndexType functionIndex1) = 0;

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

protected:
  // for quadrature integration
  std::vector<Scalar> weights;
  std::vector<Vector> points;

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
struct PolynomialPrecomputer;

template<typename PolynomialSpace>
struct PolynomialPrecomputer<Space2, PolynomialSpace>: public PolynomialPrecomputerCommon<Space2, PolynomialSpace>
{
  SPACE2_TYPEDEFS
  typedef Space2 Space;

  using PolynomialPrecomputerCommon<Space, PolynomialSpace>::functionsCount;
  using PolynomialPrecomputerCommon<Space, PolynomialSpace>::order;
  using PolynomialPrecomputerCommon<Space, PolynomialSpace>::space;

  PolynomialPrecomputer()
  {
    // this constant sets number of gauss-legendre-lobatto points, i.e. accuracy of integration
    // http://people.sc.fsu.edu/~jburkardt/cpp_src/triangle_fekete_rule/triangle_fekete_rule.html/

    const int Precisions[] = {3, 6, 9, 12, 12, 15, 18};
    int rule = 1;

    while (rule <= 7 && order > Precisions[rule - 1])
      ++rule;

    /* Rule Precision 
       1    3
       2    6
       3    9
       4-5  12
       6    15
       7    18 */

    int rule_num = ::fekete_rule_num();
    assert(rule < rule_num);

    int order_num = ::fekete_order_num(rule);

    Scalar* xy = new Scalar[2 * order_num];
    Scalar* w = new Scalar[order_num];

    ::fekete_rule(rule, order_num, xy, w);

    for (int i = 0; i < order_num; ++i)
    {
      weights.push_back(w[i] * Scalar(0.5) /* area of unit triangle */);
      points.push_back(Vector(xy[2 * i], xy[2 * i + 1]));
    }

    delete[] w;
    delete[] xy;
    ComputeVolumeIntegrals();
  }

  Scalar GetBasisFunctionDerivative(Vector point, IndexVector derivatives, IndexType functionIndex)
  {
    Polynomial<Scalar, IndexType, 2> function = space.GetBasisPolynomial(functionIndex);

    Scalar coords[2];
    coords[0] = point.x;
    coords[1] = point.y;

    IndexType derivativeDegrees[2];
    derivativeDegrees[0] = derivatives.x;
    derivativeDegrees[1] = derivatives.y;
    function.Differentiate(derivativeDegrees);

    return function.GetValue(coords);
  }

  Scalar ComputeCellVolumeIntegral(IndexType functionIndex) 
  {
    Polynomial<Scalar, IndexType, 2> function = space.GetBasisPolynomial(functionIndex);
    return function.ComputeSubspaceIntegral(2);
  }

  Scalar ComputeCellVolumeIntegral(IndexType functionIndex0, IndexType functionIndex1) // <Ф1, Ф2>
  {
    Polynomial<Scalar, IndexType, 2> function0 = space.GetBasisPolynomial(functionIndex0);
    Polynomial<Scalar, IndexType, 2> function1 = space.GetBasisPolynomial(functionIndex1);

    Polynomial<Scalar, IndexType, 2> correlation = function0 * function1;
    return correlation.ComputeSubspaceIntegral(2);
    /*

    Scalar res = 0;

    Scalar step = 5e-3f;
    for(Scalar x = 0; x < Scalar(1.0); x += step)
    {
      for(Scalar y = 0; y < Scalar(1.0); y += step)
      {
        for(Scalar z = 0; z < Scalar(1.0); z += step)
        {
          if(x + y + z < Scalar(1.0))
          {
            res += GetBasisFunctionValue(Vector(x, y, z), functionIndex0) *
                   GetBasisFunctionValue(Vector(x, y, z), functionIndex1) *
                   step * step * step;
          }
        }
      }
    }
    return res;*/
  }

  Vector ComputeDerivativeVolumeIntegral(IndexType functionIndex0, IndexType functionIndex1) // {<Ф2, dФ1/dx>, <Ф2, dФ1/dy>, <Ф2, dФ1/dz>}
  {
    Polynomial<Scalar, IndexType, 2> function0 = space.GetBasisPolynomial(functionIndex0);
    Polynomial<Scalar, IndexType, 2> function1 = space.GetBasisPolynomial(functionIndex1);

    Vector res;

    IndexType derivatives[2];

    derivatives[0] = 1;
    derivatives[1] = 0;
    Polynomial<Scalar, IndexType, 2> dfunction1dx = function1.Differentiate(derivatives);

    derivatives[0] = 0;
    derivatives[1] = 1;
    Polynomial<Scalar, IndexType, 2> dfunction1dy = function1.Differentiate(derivatives);

    return Vector((function0 * dfunction1dx).ComputeSubspaceIntegral(2), (function0 * dfunction1dy).ComputeSubspaceIntegral(2));
    /*
    Vector res(0, 0, 0);

    Scalar step = 1e-2f;
    for(Scalar x = 0; x < Scalar(1.0); x += step)
    {
      for(Scalar y = 0; y < Scalar(1.0); y += step)
      {
        for(Scalar z = 0; z < Scalar(1.0); z += step)
        {
          if(x + y + z < Scalar(1.0))
          {
            res.x +=  GetBasisFunctionValue(Vector(x, y, z), functionIndex0) *
                      GetBasisFunctionDerivative(Vector(x, y, z), IndexVector(1, 0, 0), functionIndex1) *
                      step * step * step;
            res.y +=  GetBasisFunctionValue(Vector(x, y, z), functionIndex0) *
                      GetBasisFunctionDerivative(Vector(x, y, z), IndexVector(0, 1, 0), functionIndex1) *
                      step * step * step;
            res.z +=  GetBasisFunctionValue(Vector(x, y, z), functionIndex0) *
                      GetBasisFunctionDerivative(Vector(x, y, z), IndexVector(0, 0, 1), functionIndex1) *
                      step * step * step;
          }
        }
      }
    }
    return res;*/
  }

  Scalar ComputeOutgoingFlux(IndexType refEdgeNumber, IndexType basisFunction0, IndexType basisFunction1)
  {
    /*Scalar step = Scalar(1e-5) * 10;
    Scalar res = 0;
    for(Scalar x = 0; x < Scalar(1.0); x += step)
    {
      res += 
        GetBasisFunctionValue(TriCell::GetEdgeRefCoords<GeomSpace>(refEdgeNumber, x), basisFunction1) *
        GetBasisFunctionValue(TriCell::GetEdgeRefCoords<GeomSpace>(refEdgeNumber, x), basisFunction0) * step;
    } 
    return res;*/

    Polynomial<Scalar, IndexType, 2> function0 = space.GetBasisPolynomial(basisFunction0);
    Polynomial<Scalar, IndexType, 2> function1 = space.GetBasisPolynomial(basisFunction1);


    Polynomial<Scalar, IndexType, 1> edgeToRef[2];
    Cell<Space2>::GetEdgeRefTransitionPolynomials(refEdgeNumber, edgeToRef);

    typename Polynomial<Scalar, IndexType, 1>::Pows pows;

    Polynomial<Scalar, IndexType, 1> intToEdge[1];
    pows.pows[0] = 1;
    intToEdge[0].AddTerm(pows, 1);

    Polynomial<Scalar, IndexType, 1> intToRef[2];

    for(IndexType coordIndex = 0; coordIndex < 2; coordIndex++)
    {
      intToRef[coordIndex] = edgeToRef[coordIndex].Substitute(intToEdge);
    }

    Polynomial<Scalar, IndexType, 1> intToSrcFunc;
    Polynomial<Scalar, IndexType, 1> intToDstFunc;

    intToSrcFunc = function1.Substitute(intToRef);
    intToDstFunc = function0.Substitute(intToRef);

    Scalar integrationResult = (intToSrcFunc * intToDstFunc).ComputeSubspaceIntegral(1);
    //if(fabs(res - integrationResult) > 1e-3) printf("flux poo");
    return integrationResult;
  }

  Scalar ComputeIncomingFlux(
    IndexType refEdgeNumber, IndexType correspondingEdgeNumber,
    IndexType basisFunction0, IndexType basisFunction1)
  {
    /*if(refEdgeNumber == 0 && correspondingEdgeNumber == 1 && basisFunction0 == 0 && basisFunction1 == 1)
    {
      int pp = 1;
    }
    Scalar step = Scalar(1e-5) * 10;
    Scalar res = 0;
    for(Scalar x = 0; x < Scalar(1.0); x += step)
    {
      Scalar srcEdgeCoords = x;
      Scalar dstEdgeCoords = TriCell::GetCorrespondingEdgeCoord<GeomSpace>(srcEdgeCoords);

      Vector2 srcRefCoords   = TriCell::GetEdgeRefCoords<GeomSpace>(refEdgeNumber,           srcEdgeCoords);
      Vector2 dstRefCoords   = TriCell::GetEdgeRefCoords<GeomSpace>(correspondingEdgeNumber, dstEdgeCoords);

      res += 
        GetBasisFunctionValue(dstRefCoords, basisFunction0) *
        GetBasisFunctionValue(srcRefCoords, basisFunction1) * step;
    } */
    //return res;

    Polynomial<Scalar, IndexType, 2> function0 = space.GetBasisPolynomial(basisFunction0);
    Polynomial<Scalar, IndexType, 2> function1 = space.GetBasisPolynomial(basisFunction1);


    Polynomial<Scalar, IndexType, 1> srcEdgeToRef[2];
    Cell<Space2>::GetEdgeRefTransitionPolynomials(refEdgeNumber          , srcEdgeToRef);

    Polynomial<Scalar, IndexType, 1> dstEdgeToRef[2];
    Cell<Space2>::GetEdgeRefTransitionPolynomials(correspondingEdgeNumber, dstEdgeToRef);


    typename Polynomial<Scalar, IndexType, 1>::Pows pows;

    Polynomial<Scalar, IndexType, 1> intToSrcEdge[1];
    pows.pows[0] = 1;
    intToSrcEdge[0].AddTerm(pows, 1);

    Polynomial<Scalar, IndexType, 1> srcEdgeToDstEdge[1];
    Cell<Space2>::GetCorrespondingEdgeTransitionPolynomials(srcEdgeToDstEdge);

    Polynomial<Scalar, IndexType, 1> intToDstEdge[1];
    intToDstEdge[0] = srcEdgeToDstEdge[0].Substitute(intToSrcEdge);

    Polynomial<Scalar, IndexType, 1> intToSrcRef[2];
    Polynomial<Scalar, IndexType, 1> intToDstRef[2];

    for(IndexType coordIndex = 0; coordIndex < 2; coordIndex++)
    {
      intToSrcRef[coordIndex] = srcEdgeToRef[coordIndex].Substitute(intToSrcEdge);
      intToDstRef[coordIndex] = dstEdgeToRef[coordIndex].Substitute(intToDstEdge);
    }

    Polynomial<Scalar, IndexType, 1> intToSrcFunc;
    Polynomial<Scalar, IndexType, 1> intToDstFunc;

    intToSrcFunc = function1.Substitute(intToSrcRef);
    intToDstFunc = function0.Substitute(intToDstRef);

    Scalar integrationResult = (intToSrcFunc * intToDstFunc).ComputeSubspaceIntegral(1);
    //if(fabs(res - integrationResult) > 1e-3) printf("flux poo");
    return integrationResult;
  }

  Scalar ComputeEdgeFlux(IndexType refEdgeNumber, IndexType basisFunction)
  {
    Polynomial<Scalar, IndexType, 2> function = space.GetBasisPolynomial(basisFunction);

    Polynomial<Scalar, IndexType, 1> edgeToRef[2];
    Cell<Space2>::GetEdgeRefTransitionPolynomials(refEdgeNumber, edgeToRef);

    typename Polynomial<Scalar, IndexType, 1>::Pows pows;

    Polynomial<Scalar, IndexType, 1> intToEdge[1];
    pows.pows[0] = 1;
    intToEdge[0].AddTerm(pows, 1);

    Polynomial<Scalar, IndexType, 1> intToRef[2];

    for(IndexType coordIndex = 0; coordIndex < 2; coordIndex++)
    {
      intToRef[coordIndex] = edgeToRef[coordIndex].Substitute(intToEdge);
    }

    Polynomial<Scalar, IndexType, 1> intToFunc;

    intToFunc = function.Substitute(intToRef);

    Scalar integrationResult = (intToFunc).ComputeSubspaceIntegral(1);
    //if(fabs(res - integrationResult) > 1e-3) printf("flux poo");
    return integrationResult;
  }
};

template<typename PolynomialSpace>
struct PolynomialPrecomputer<Space3, PolynomialSpace>: public PolynomialPrecomputerCommon<Space3, PolynomialSpace>
{
  SPACE3_TYPEDEFS
  typedef Space3 Space;

  using PolynomialPrecomputerCommon<Space, PolynomialSpace>::functionsCount;
  using PolynomialPrecomputerCommon<Space, PolynomialSpace>::order;
  using PolynomialPrecomputerCommon<Space, PolynomialSpace>::space;

  PolynomialPrecomputer(): PolynomialPrecomputerCommon<Space3, PolynomialSpace>()
  {
    // from http://people.sc.fsu.edu/~jburkardt/c_src/tetrahedron_arbq_rule/tetrahedron_arbq_rule.html

    int order_num = ::tetrahedron_arbq_size(order + 1);

    Scalar* xyz = new Scalar[3 * order_num];
    Scalar* w = new Scalar[order_num];

    ::tetrahedron_arbq(order + 1, order_num, xyz, w);

    for (int i = 0; i < order_num; ++i)
    {
      weights.push_back(w[i] / (4 * sqrt(2.0)));
      Scalar* uvw = ::ref_to_koorn(xyz + 3 * i);
      points.push_back(Vector(uvw[0] + 1, uvw[1] + 1, uvw[2] + 1) * Scalar(0.5));

      free(uvw);
    }

    delete[] w;
    delete[] xyz;
    ComputeVolumeIntegrals();
  }

  Scalar GetBasisFunctionDerivative(Vector point, IndexVector derivatives, IndexType functionIndex)
  {
    Polynomial<Scalar, IndexType, 3> function = space.GetBasisPolynomial(functionIndex);

    Scalar coords[3];
    coords[0] = point.x;
    coords[1] = point.y;
    coords[2] = point.z;

    IndexType derivativeDegrees[3];
    derivativeDegrees[0] = derivatives.x;
    derivativeDegrees[1] = derivatives.y;
    derivativeDegrees[2] = derivatives.z;
    function.Differentiate(derivativeDegrees);

    return function.GetValue(coords);
  }

  Scalar ComputeCellVolumeIntegral(IndexType functionIndex0, IndexType functionIndex1) // <Ô1, Ô2>
  {
    Polynomial<Scalar, IndexType, 3> function0 = space.GetBasisPolynomial(functionIndex0);
    Polynomial<Scalar, IndexType, 3> function1 = space.GetBasisPolynomial(functionIndex1);

    Polynomial<Scalar, IndexType, 3> correlation = function0 * function1;
    return correlation.ComputeSubspaceIntegral(3);

    /*

    Scalar res = 0;

    Scalar step = 5e-3f;
    for(Scalar x = 0; x < Scalar(1.0); x += step)
    {
    for(Scalar y = 0; y < Scalar(1.0); y += step)
    {
    for(Scalar z = 0; z < Scalar(1.0); z += step)
    {
    if(x + y + z < Scalar(1.0))
    {
    res += GetBasisFunctionValue(Vector(x, y, z), functionIndex0) *
    GetBasisFunctionValue(Vector(x, y, z), functionIndex1) *
    step * step * step;
    }
    }
    }
    }
    return res;*/
  }

  Vector ComputeDerivativeVolumeIntegral(IndexType functionIndex0, IndexType functionIndex1) // {<Ô2, dÔ1/dx>, <Ô2, dÔ1/dy>, <Ô2, dÔ1/dz>}
  {
    Polynomial<Scalar, IndexType, 3> function0 = space.GetBasisPolynomial(functionIndex0);
    Polynomial<Scalar, IndexType, 3> function1 = space.GetBasisPolynomial(functionIndex1);

    Vector res;

    IndexType derivatives[3];

    derivatives[0] = 1;
    derivatives[1] = 0;
    derivatives[2] = 0;
    Polynomial<Scalar, IndexType, 3> dfunction1dx = function1.Differentiate(derivatives);

    derivatives[0] = 0;
    derivatives[1] = 1;
    derivatives[2] = 0;
    Polynomial<Scalar, IndexType, 3> dfunction1dy = function1.Differentiate(derivatives);

    derivatives[0] = 0;
    derivatives[1] = 0;
    derivatives[2] = 1;
    Polynomial<Scalar, IndexType, 3> dfunction1dz = function1.Differentiate(derivatives);

    return Vector((function0 * dfunction1dx).ComputeSubspaceIntegral(3),
      (function0 * dfunction1dy).ComputeSubspaceIntegral(3),
      (function0 * dfunction1dz).ComputeSubspaceIntegral(3));
    /*
    Vector res(0, 0, 0);

    Scalar step = 1e-2f;
    for(Scalar x = 0; x < Scalar(1.0); x += step)
    {
    for(Scalar y = 0; y < Scalar(1.0); y += step)
    {
    for(Scalar z = 0; z < Scalar(1.0); z += step)
    {
    if(x + y + z < Scalar(1.0))
    {
    res.x +=  GetBasisFunctionValue(Vector(x, y, z), functionIndex0) *
    GetBasisFunctionDerivative(Vector(x, y, z), IndexVector(1, 0, 0), functionIndex1) *
    step * step * step;
    res.y +=  GetBasisFunctionValue(Vector(x, y, z), functionIndex0) *
    GetBasisFunctionDerivative(Vector(x, y, z), IndexVector(0, 1, 0), functionIndex1) *
    step * step * step;
    res.z +=  GetBasisFunctionValue(Vector(x, y, z), functionIndex0) *
    GetBasisFunctionDerivative(Vector(x, y, z), IndexVector(0, 0, 1), functionIndex1) *
    step * step * step;
    }
    }
    }
    }
    return res;*/
  }

  Scalar ComputeOutgoingFlux(IndexType refFaceNumber, IndexType basisFunction0, IndexType basisFunction1)
  {
    /*Scalar step = 1e-2f;
    Scalar res = 0;
    for(Scalar x = 0; x < Scalar(1.0); x += step)
    {
    for(Scalar y = 0; y < Scalar(1.0); y += step)
    {
    if(x + y < Scalar(1.0))
    {
    res +=
    functionSpace.GetBasisFunctionValue(TetraCell::GetSurfaceRefCoords(refFaceNumber, Vector2(x, y)), basisFunction1) *
    functionSpace.GetBasisFunctionValue(TetraCell::GetSurfaceRefCoords(refFaceNumber, Vector2(x, y)), basisFunction0) * step * step;
    }
    }
    }
    return res;*/

    Polynomial<Scalar, IndexType, 3> function0 = space.GetBasisPolynomial(basisFunction0);
    Polynomial<Scalar, IndexType, 3> function1 = space.GetBasisPolynomial(basisFunction1);


    Polynomial<Scalar, IndexType, 2> faceToRef[3];
    Cell<Space3>::GetFaceRefTransitionPolynomials(refFaceNumber, faceToRef);

    Polynomial<Scalar, IndexType, 2>::Pows pows;

    Polynomial<Scalar, IndexType, 2> intToFace[2];

    pows.pows[0] = 1;
    pows.pows[1] = 0;
    intToFace[0].AddTerm(pows, 1);

    pows.pows[0] = 0;
    pows.pows[1] = 1;
    intToFace[1].AddTerm(pows, 1);

    Polynomial<Scalar, IndexType, 2> intToRef[3];

    for (IndexType coordIndex = 0; coordIndex < 3; coordIndex++)
    {
      intToRef[coordIndex] = faceToRef[coordIndex].Substitute(intToFace);
    }

    Polynomial<Scalar, IndexType, 2> intToSrcFunc;
    Polynomial<Scalar, IndexType, 2> intToDstFunc;

    intToSrcFunc = function1.Substitute(intToRef);
    intToDstFunc = function0.Substitute(intToRef);

    Scalar integrationResult = (intToSrcFunc * intToDstFunc).ComputeSubspaceIntegral(2);
    //if(fabs(res - integrationResult) > 1e-3) printf("flux poo");
    return integrationResult;
  }

  Scalar ComputeIncomingFlux(
    IndexType refFaceNumber, IndexType correspondingFaceNumber, IndexType orientation,
    IndexType basisFunction0, IndexType basisFunction1)
  {
    /*Scalar step = 1e-2f;
    Scalar res = 0;
    for(Scalar x = 0; x < Scalar(1.0); x += step)
    {
    for(Scalar y = 0; y < Scalar(1.0); y += step)
    {
    if(x + y < Scalar(1.0))
    {
    Vector2 srcSurfaceCoords = Vector2(x, y);
    Vector2 dstSurfaceCoords = TetraCell::GetCorrespondingSurfaceCoords(srcSurfaceCoords, orientation);

    Vector srcRefCoords   = TetraCell::GetSurfaceRefCoords(refFaceNumber,           srcSurfaceCoords);
    Vector dstRefCoords   = TetraCell::GetSurfaceRefCoords(correspondingFaceNumber, dstSurfaceCoords);

    res +=
    functionSpace.GetBasisFunctionValue(dstRefCoords, basisFunction0) *
    functionSpace.GetBasisFunctionValue(srcRefCoords, basisFunction1) * step * step;
    }
    }
    }
    return res;*/
    /*if(refEdgeNumber == 0 && correspondingEdgeNumber == 1 && basisFunction0 == 0 && basisFunction1 == 1)
    {
    int pp = 1;
    }
    Scalar step = Scalar(1e-5) * 10;
    Scalar res = 0;
    for(Scalar x = 0; x < Scalar(1.0); x += step)
    {
    Scalar srcEdgeCoords = x;
    Scalar dstEdgeCoords = TriCell::GetCorrespondingEdgeCoord<GeomSpace>(srcEdgeCoords);

    Vector2 srcRefCoords   = TriCell::GetEdgeRefCoords<GeomSpace>(refEdgeNumber,           srcEdgeCoords);
    Vector2 dstRefCoords   = TriCell::GetEdgeRefCoords<GeomSpace>(correspondingEdgeNumber, dstEdgeCoords);

    res +=
    GetBasisFunctionValue(dstRefCoords, basisFunction0) *
    GetBasisFunctionValue(srcRefCoords, basisFunction1) * step;
    } */
    //return res;

    Polynomial<Scalar, IndexType, 3> function0 = space.GetBasisPolynomial(basisFunction0);
    Polynomial<Scalar, IndexType, 3> function1 = space.GetBasisPolynomial(basisFunction1);


    Polynomial<Scalar, IndexType, 2> srcFaceToRef[3];
    Cell<Space3>::GetFaceRefTransitionPolynomials(refFaceNumber, srcFaceToRef);

    Polynomial<Scalar, IndexType, 2> dstFaceToRef[3];
    Cell<Space3>::GetFaceRefTransitionPolynomials(correspondingFaceNumber, dstFaceToRef);


    Polynomial<Scalar, IndexType, 2>::Pows pows;

    Polynomial<Scalar, IndexType, 2> intToSrcFace[2];

    pows.pows[0] = 1;
    pows.pows[1] = 0;
    intToSrcFace[0].AddTerm(pows, 1);

    pows.pows[0] = 0;
    pows.pows[1] = 1;
    intToSrcFace[1].AddTerm(pows, 1);

    Polynomial<Scalar, IndexType, 2> srcFaceToDstFace[2];
    Cell<Space3>::GetCorrespondingFaceTransitionPolynomials(orientation, srcFaceToDstFace);

    Polynomial<Scalar, IndexType, 2> intToDstFace[2];
    intToDstFace[0] = srcFaceToDstFace[0].Substitute(intToSrcFace);
    intToDstFace[1] = srcFaceToDstFace[1].Substitute(intToSrcFace);

    Polynomial<Scalar, IndexType, 2> intToSrcRef[3];
    Polynomial<Scalar, IndexType, 2> intToDstRef[3];

    for (IndexType coordIndex = 0; coordIndex < 3; coordIndex++)
    {
      intToSrcRef[coordIndex] = srcFaceToRef[coordIndex].Substitute(intToSrcFace);
      intToDstRef[coordIndex] = dstFaceToRef[coordIndex].Substitute(intToDstFace);
    }

    Polynomial<Scalar, IndexType, 2> intToSrcFunc;
    Polynomial<Scalar, IndexType, 2> intToDstFunc;

    intToSrcFunc = function1.Substitute(intToSrcRef);
    intToDstFunc = function0.Substitute(intToDstRef);

    Scalar integrationResult = (intToSrcFunc * intToDstFunc).ComputeSubspaceIntegral(2);
    //if(fabs(res - integrationResult) > 1e-3) printf("flux poo");
    return integrationResult;
  }
};

