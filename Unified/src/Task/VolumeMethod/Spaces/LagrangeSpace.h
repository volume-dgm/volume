#pragma once
#include "../../../Maths/Spaces.h"
#include "BaseSpace.h"

template<typename Space, int Order>
struct LagrangeSpace;

template<int Order>
struct LagrangeSpace<Space2, Order>: public BaseSpace<Space2, Order>
{
  SPACE2_TYPEDEFS

  using BaseSpace<Space2, Order>::order;
  using BaseSpace<Space2, Order>::functionsCount;

  Scalar functionMult;

  typedef Polynomial<Scalar, IndexType, 2> PolynomialType;

  LagrangeSpace(): BaseSpace<Space2, Order>()
  {
    ComputeCoordPowers();
  }

  Polynomial<Scalar, IndexType, 2> GetBasisPolynomial(IndexType functionIndex)
  {
    return GetBasisFunctionPolynomial(functionIndex);
    /*Vector2i polynomialPows = GetCoordPowers(functionIndex);

    typename Polynomial<Scalar, IndexType, 2>::Pows pows;
    pows.pows[0] = polynomialPows.x;
    pows.pows[1] = polynomialPows.y;

    Polynomial<Scalar, IndexType, 2> res;
    res.AddTerm(pows, Scalar(1.0));
    return res;*/
  }


  PolynomialType GetBasisFunctionPolynomial(IndexType functionIndex)
  {
    Vector localCoordGradients[3];

    Vector vertices[3];
    vertices[0] = Vector(0, 0);
    vertices[1] = Vector(Scalar(1.0), 0);
    vertices[2] = Vector(0, Scalar(1.0));

    Scalar square = Scalar(0.5) * (vertices[1] - vertices[0]) ^ (vertices[2] - vertices[0]);
    Scalar invSquare = Scalar(1.0) / square;

    Vector dirs[3];
    dirs[0] = vertices[2] - vertices[1];
    dirs[1] = vertices[0] - vertices[2];
    dirs[2] = vertices[1] - vertices[0];

    localCoordGradients[0] = Vector(-dirs[0].y, dirs[0].x) * Scalar(0.5) * invSquare;
    localCoordGradients[1] = Vector(-dirs[1].y, dirs[1].x) * Scalar(0.5) * invSquare;
    localCoordGradients[2] = Vector(-dirs[2].y, dirs[2].x) * Scalar(0.5) * invSquare;

    PolynomialType localCoordPolynomials[3];

    typename PolynomialType::Pows defPows;
    defPows.pows[0] = 0;
    defPows.pows[1] = 0;

    localCoordPolynomials[2].AddTerm(defPows, Scalar(1.0));
    for(IndexType coordIndex = 0; coordIndex < 2; coordIndex++)
    {
      typename PolynomialType::Pows termPows;

      termPows.pows[0] = termPows.pows[1] = 0;
      localCoordPolynomials[coordIndex].AddTerm(termPows, -localCoordGradients[coordIndex] * vertices[coordIndex + 1]);

      termPows.pows[0] = 1;
      termPows.pows[1] = 0;
      localCoordPolynomials[coordIndex].AddTerm(termPows, localCoordGradients[coordIndex].x);

      termPows.pows[0] = 0;
      termPows.pows[1] = 1;
      localCoordPolynomials[coordIndex].AddTerm(termPows, localCoordGradients[coordIndex].y);

      localCoordPolynomials[2] -= localCoordPolynomials[coordIndex];
    }

    IndexType coordPowers[3];
    GetCoordPowers(functionIndex, coordPowers);

    Scalar denom(1.0);

    for (IndexType coordIndex = 0; coordIndex < 3; coordIndex++)
    {
      for(IndexType termNumber = 0; termNumber < coordPowers[coordIndex]; termNumber++)
      {
        denom /= Scalar(termNumber + 1) / Scalar(Order);
      }
    }

    PolynomialType resPolynomial;

    resPolynomial.AddTerm(defPows, denom);

    for (IndexType coordIndex = 0; coordIndex < 3; coordIndex++)
    {
      for(IndexType termNumber = 0; termNumber < coordPowers[coordIndex]; termNumber++)
      {
        PolynomialType localTerm;

        typename PolynomialType::Pows localPows;
        localPows.pows[0] = 0;
        localPows.pows[1] = 0;
        localTerm.AddTerm(localPows, -Scalar(termNumber) / Scalar(Order));

        localTerm += localCoordPolynomials[coordIndex];

        resPolynomial = resPolynomial * localTerm;
      }
    }

    return resPolynomial;
  }

  Vector GetBasisPoint(IndexType functionIndex)
  {
    IndexType coordPowers[3];
    GetCoordPowers(functionIndex, coordPowers);

    Vector vertices[3];
    vertices[0] = Vector(0, 0);
    vertices[1] = Vector(Scalar(1.0), 0);
    vertices[2] = Vector(0, Scalar(1.0));

    return vertices[0] 
      + (vertices[1] - vertices[0]) * Scalar(coordPowers[1]) / Scalar(Order)
      + (vertices[2] - vertices[0]) * Scalar(coordPowers[2]) / Scalar(Order);
  }

  Scalar GetBasisFunctionValue(const Vector& point, IndexType functionIndex) const
  {
    IndexType coordPowers[3];
    GetCoordPowers(functionIndex, coordPowers);

    Scalar localCoords[3];
    GetLocalCoords(point, localCoords);

    return GetPolynomeValue(coordPowers, localCoords);

    /*PolynomialType polynomial = GetBasisFunctionPolynomial(functionIndex);
    Scalar polynomialRes = polynomial.GetValue(&(point.x));

    if(fabs(normalRes - polynomialRes) > 1e-5) printf("poo");
    return polynomialRes;*/
  }

  Scalar GetBasisFunctionMaxValue(IndexType functionIndex)
  {
    return Scalar(1.0);
  }

  Vector GetBasisFunctionGradient(const Vector& point, IndexType functionIndex)
  {
    Scalar smallStep(1e-5);
    Vector resGradient(
      GetBasisFunctionValue(point + Vector2( smallStep, 0), functionIndex) - GetBasisFunctionValue(point + Vector2(-smallStep,  0), functionIndex),
      GetBasisFunctionValue(point + Vector2( 0, smallStep), functionIndex) - GetBasisFunctionValue(point + Vector2( 0, -smallStep), functionIndex));
    resGradient *= Scalar(1.0) / (Scalar(2.0) * smallStep);
//    return resGradient;
    Vector localCoordGradients[3];

    Vector vertices[3];
    vertices[0] = Vector(0, 0);
    vertices[1] = Vector(Scalar(1.0), 0);
    vertices[2] = Vector(0, Scalar(1.0));

    Scalar square = Scalar(0.5) * (vertices[1] - vertices[0]) ^ (vertices[2] - vertices[0]);
    Scalar invSquare = Scalar(1.0) / square;

    Vector dirs[3];
    dirs[0] = vertices[2] - vertices[1];
    dirs[1] = vertices[0] - vertices[2];
    dirs[2] = vertices[1] - vertices[0];

    localCoordGradients[0] = Vector(-dirs[0].y, dirs[0].x) * Scalar(0.5) * invSquare;
    localCoordGradients[1] = Vector(-dirs[1].y, dirs[1].x) * Scalar(0.5) * invSquare;
    localCoordGradients[2] = Vector(-dirs[2].y, dirs[2].x) * Scalar(0.5) * invSquare;

    IndexType coordPowers[3];
    GetCoordPowers(functionIndex, coordPowers);

    Scalar localCoords[3];
    GetLocalCoords(point, localCoords);

    Scalar mult(1.0);

    for(IndexType coordIndex = 0; coordIndex < 3; coordIndex++)
    {
      for(IndexType power = 0; power < coordPowers[coordIndex]; power++)
      {
        mult *= Scalar(power + 1) / Scalar(Order);
      }
    }

    Vector gradient(0, 0);
    for(IndexType coordIndex0 = 0; coordIndex0 < 3; coordIndex0++)
    {
      for(IndexType power0 = 0; power0 < coordPowers[coordIndex0]; power0++)
      {
        Scalar term = mult;
        for(IndexType coordIndex1 = 0; coordIndex1 < 3; coordIndex1++)
        {
          for(IndexType power1 = 0; power1 < coordPowers[coordIndex1]; power1++)
          {
            if(coordIndex0 != coordIndex1 || power0 != power1)
            {
              term *= localCoords[coordIndex1] - Scalar(power1) / Scalar(Order);
            }
          }
        }
        gradient += localCoordGradients[coordIndex0] * term;
      }
    }
    return gradient;
  }

  Scalar ComputeCellVolumeIntegral(IndexType functionIndex)
  {
    /*Scalar res = 0;

    Scalar step = Scalar(1e-3);
    for(Scalar x = 0; x < Scalar(1.0); x += step)
    {
      for(Scalar y = 0; y < Scalar(1.0); y += step)
      {
        if(x + y < Scalar(1.0))
        {
          res += GetBasisFunctionValue(Vector2(x, y), functionIndex) *          
                 step * step;
        }
      }
    }*/

    PolynomialType polynomial = GetBasisFunctionPolynomial(functionIndex);   
    Scalar polynomialRes = polynomial.ComputeSubspaceIntegral(2); //?

    //printf("CellVolume1dIntegrals: error - %f, functionIndex - %d", fabs(res - polynomialRes), functionIndex);

    return polynomialRes;
  }

  Scalar ComputeCellVolumeIntegral(IndexType functionIndex0, IndexType functionIndex1) // <Ô1, Ô2>
  {
    /*Scalar res = 0;

    Scalar step = Scalar(1e-3);
    for(Scalar x = 0; x < Scalar(1.0); x += step)
    {
      for(Scalar y = 0; y < Scalar(1.0); y += step)
      {
        if(x + y < Scalar(1.0))
        {
          res += GetBasisFunctionValue(Vector2(x, y), functionIndex0) * 
                 GetBasisFunctionValue(Vector2(x, y), functionIndex1)           
                 step * step;
        }
      }
    } */

    PolynomialType polynomial0 = GetBasisFunctionPolynomial(functionIndex0);
    PolynomialType polynomial1 = GetBasisFunctionPolynomial(functionIndex1);

    PolynomialType resultPolynomial = polynomial0 * polynomial1;

    Scalar polynomialRes = resultPolynomial.ComputeSubspaceIntegral(2);

    return polynomialRes;
  }

  Vector ComputeDerivativeVolumeIntegral(IndexType functionIndex0, IndexType functionIndex1) // {<dÔ1/dx, Ô2>, <dÔ1/dy, Ô2>, <dÔ1/dz, Ô2>}
  {
    /*Vector2 res(0, 0);

    Scalar step = Scalar(1e-3);
    for(Scalar x = 0; x < Scalar(1.0); x += step)
    {
      for(Scalar y = 0; y < Scalar(1.0); y += step)
      {
        if(x + y < Scalar(1.0))
        {
          Vector2 gradient = GetBasisFunctionGradient(Vector2(x, y), functionIndex1);
          res.x +=  GetBasisFunctionValue(Vector2(x, y), functionIndex0) *
                    gradient.x * step * step;
          res.y +=  GetBasisFunctionValue(Vector2(x, y), functionIndex0) *
                    gradient.y * step * step;
        }
      }
    }
    Vector2 normalRes = res;*/

    PolynomialType polynomial0 = GetBasisFunctionPolynomial(functionIndex0);
    PolynomialType polynomial1 = GetBasisFunctionPolynomial(functionIndex1);

    IndexType derivatives[2];

    derivatives[0] = 1;
    derivatives[1] = 0;
    PolynomialType polynomial1xDerivative = polynomial1.Differentiate(derivatives);

    derivatives[0] = 0;
    derivatives[1] = 1;
    PolynomialType polynomial1yDerivative = polynomial1.Differentiate(derivatives);

    PolynomialType resultProductx = polynomial0 * polynomial1xDerivative;
    PolynomialType resultProducty = polynomial0 * polynomial1yDerivative;


    Vector polynomialRes = Vector(resultProductx.ComputeSubspaceIntegral(2), resultProducty.ComputeSubspaceIntegral(2));
    return polynomialRes;
  }

public:
  void GetLocalCoords(const Vector& point, Scalar* localCoords) const
  {
    Vector vertices[3];
    vertices[0] = Vector(0, 0);
    vertices[1] = Vector(Scalar(1.0), 0);
    vertices[2] = Vector(0, Scalar(1.0));
    Scalar square = Scalar(0.5) * (vertices[1] - vertices[0]) ^ (vertices[2] - vertices[0]);
    Scalar invSquare = Scalar(1.0) / square;

    localCoords[0] = ((vertices[2] - vertices[1]) ^ (point - vertices[1])) * invSquare * Scalar(0.5);
    localCoords[1] = ((vertices[0] - vertices[2]) ^ (point - vertices[2])) * invSquare * Scalar(0.5);
    localCoords[2] = ((vertices[1] - vertices[0]) ^ (point - vertices[0])) * invSquare * Scalar(0.5);
  }


  Scalar GetPolynomeValue(IndexType* pows, Scalar* localCoords) const
  {
    Scalar res(1.0);

    for(IndexType coordIndex = 0; coordIndex < 3; coordIndex++)
    {
      for(IndexType power = 0; power < pows[coordIndex]; power++)
      {
        res /= Scalar(power + 1) / Scalar(Order);
      }
    }

    Scalar denom = res;

    for(IndexType coordIndex = 0; coordIndex < 3; coordIndex++)
    {
      for(IndexType power = 0; power < pows[coordIndex]; power++)
      {
        res *= localCoords[coordIndex] - Scalar(power) / Scalar(Order);
      }
    }

    Scalar baseRes = res;
    return res;
  }

  void ComputeCoordPowers()
  {
    IndexType functionIndex = 0;
    for(IndexType i = 0; i <= Order; i++)
    {
      for(IndexType j = 0; j <= Order - i; j++)
      {
        IndexType k = Order - i - j;
        functionPows[functionIndex][0] = i;
        functionPows[functionIndex][1] = j;
        functionPows[functionIndex][2] = k;
        functionIndex++;
      }
    }
  }

  void GetCoordPowers(IndexType functionIndex, IndexType* pows) const
  {
    for(IndexType i = 0; i < 3; i++)
    {
      pows[i] = functionPows[functionIndex][i];
    }
  }
  IndexType functionPows[functionsCount][3];
};

template<int Order>
struct LagrangeSpace<Space3, Order>: public BaseSpace<Space3, Order>
{
  SPACE3_TYPEDEFS

  using BaseSpace<Space3, Order>::order;
  using BaseSpace<Space3, Order>::functionsCount;

  Scalar functionMult;

  typedef Polynomial<Scalar, IndexType, 3> PolynomialType;

  LagrangeSpace()
  {
    ComputeCoordPowers();
  }

  Polynomial<Scalar, IndexType, 3> GetBasisPolynomial(IndexType functionIndex)
  {
    return GetBasisFunctionPolynomial(functionIndex);
  }


  PolynomialType GetBasisFunctionPolynomial(IndexType functionIndex)
  {

    Vector localCoordGradients[4];

    Vector vertices[4];
    vertices[0] = Vector(0, 0, 0);
    vertices[1] = Vector(Scalar(1), 0, 0);
    vertices[2] = Vector(0, Scalar(1), 0);
    vertices[3] = Vector(0, 0, Scalar(1));

    Scalar volume = MixedProduct(vertices[1] - vertices[0], vertices[2] - vertices[0], vertices[3] - vertices[0]) / Scalar(6.0);
    Scalar invVolume = Scalar(1.0) / volume;


    localCoordGradients[0] = (vertices[3] - vertices[1]) ^ (vertices[2] - vertices[1]) / (Scalar(6.0) * volume);
    localCoordGradients[1] = (vertices[2] - vertices[0]) ^ (vertices[3] - vertices[0]) / (Scalar(6.0) * volume);
    localCoordGradients[2] = (vertices[3] - vertices[0]) ^ (vertices[1] - vertices[0]) / (Scalar(6.0) * volume);
    localCoordGradients[3] = (vertices[1] - vertices[0]) ^ (vertices[2] - vertices[0]) / (Scalar(6.0) * volume);

    PolynomialType localCoordPolynomials[4];

    PolynomialType::Pows defPows;
    defPows.pows[0] = 0;
    defPows.pows[1] = 0;
    defPows.pows[2] = 0;

    localCoordPolynomials[3].AddTerm(defPows, Scalar(1.0));
    for (IndexType coordIndex = 0; coordIndex < 3; coordIndex++)
    {
      PolynomialType::Pows termPows;

      termPows.pows[0] = termPows.pows[1] = termPows.pows[2] = 0;
      localCoordPolynomials[coordIndex].AddTerm(termPows, -localCoordGradients[coordIndex] * vertices[coordIndex + 1]);

      termPows.pows[0] = 1;
      termPows.pows[1] = 0;
      termPows.pows[2] = 0;
      localCoordPolynomials[coordIndex].AddTerm(termPows, localCoordGradients[coordIndex].x);

      termPows.pows[0] = 0;
      termPows.pows[1] = 1;
      termPows.pows[2] = 0;
      localCoordPolynomials[coordIndex].AddTerm(termPows, localCoordGradients[coordIndex].y);

      termPows.pows[0] = 0;
      termPows.pows[1] = 0;
      termPows.pows[2] = 1;
      localCoordPolynomials[coordIndex].AddTerm(termPows, localCoordGradients[coordIndex].z);

      localCoordPolynomials[3] -= localCoordPolynomials[coordIndex];
    }

    IndexType coordPowers[4];
    GetCoordPowers(functionIndex, coordPowers);

    Scalar denom(1.0);

    for (IndexType coordIndex = 0; coordIndex < 4; coordIndex++)
    {
      for (IndexType termNumber = 0; termNumber < coordPowers[coordIndex]; termNumber++)
      {
        denom /= Scalar(termNumber + 1) / Scalar(Order);
      }
    }

    PolynomialType resPolynomial;

    resPolynomial.AddTerm(defPows, denom);

    for (IndexType coordIndex = 0; coordIndex < 4; coordIndex++)
    {
      for (IndexType termNumber = 0; termNumber < coordPowers[coordIndex]; termNumber++)
      {
        PolynomialType localTerm;

        PolynomialType::Pows localPows;
        localPows.pows[0] = 0;
        localPows.pows[1] = 0;
        localPows.pows[2] = 0;
        localTerm.AddTerm(localPows, -Scalar(termNumber) / Scalar(Order));

        localTerm += localCoordPolynomials[coordIndex];

        resPolynomial = resPolynomial * localTerm;
      }
    }

    return resPolynomial;
  }

  Vector GetBasisPoint(IndexType functionIndex)
  {
    IndexType coordPowers[4];
    GetCoordPowers(functionIndex, coordPowers);

    Vector vertices[4];
    vertices[0] = Vector(0, 0, 0);
    vertices[1] = Vector(Scalar(1), 0, 0);
    vertices[2] = Vector(0, Scalar(1), 0);
    vertices[3] = Vector(0, 0, Scalar(1));

    return vertices[0]
      + (vertices[1] - vertices[0]) * Scalar(coordPowers[1]) / Scalar(Order)
      + (vertices[2] - vertices[0]) * Scalar(coordPowers[2]) / Scalar(Order)
      + (vertices[3] - vertices[0]) * Scalar(coordPowers[3]) / Scalar(Order);

  }

  Scalar GetBasisFunctionValue(const Vector& point, IndexType functionIndex) const
  {
    IndexType coordPowers[4];
    GetCoordPowers(functionIndex, coordPowers);

    Scalar localCoords[4];
    GetLocalCoords(point, localCoords);

    return GetPolynomeValue(coordPowers, localCoords);

    /*PolynomialType polynomial = GetBasisFunctionPolynomial(functionIndex);
    Scalar polynomialRes = polynomial.GetValue(&(point.x));

    if(fabs(normalRes - polynomialRes) > 1e-5) printf("poo");
    return polynomialRes;*/
  }

  Scalar GetBasisFunctionMaxValue(IndexType functionIndex)
  {
    return Scalar(1.0);
  }

  Vector GetBasisFunctionGradient(const Vector& point, IndexType functionIndex)
  {
    /*Scalar smallStep(1e-5);
    Vector2 resGradient(
    GetBasisFunctionValue(point + Vector2( smallStep, 0), functionIndex) - GetBasisFunctionValue(point + Vector2(-smallStep,  0), functionIndex),
    GetBasisFunctionValue(point + Vector2( 0, smallStep), functionIndex) - GetBasisFunctionValue(point + Vector2( 0, -smallStep), functionIndex));
    resGradient *= Scalar(1.0) / (Scalar(2.0) * smallStep);*/
    //    return resGradient;
    Vector2 localCoordGradients[4];

    Vector vertices[4];
    vertices[0] = Vector(0, 0, 0);
    vertices[1] = Vector(Scalar(1), 0, 0);
    vertices[2] = Vector(0, Scalar(1), 0);
    vertices[3] = Vector(0, 0, Scalar(1));

    Scalar volume = MixedProduct(vertices[1] - vertices[0], vertices[2] - vertices[0], vertices[3] - vertices[0]) / Scalar(6.0);
    Scalar invVolume = Scalar(1.0) / volume;

    localCoordGradients[0] = (vertices[3] - vertices[1]) ^ (vertices[2] - vertices[1]) / (Scalar(6.0) * volume);
    localCoordGradients[1] = (vertices[2] - vertices[0]) ^ (vertices[3] - vertices[0]) / (Scalar(6.0) * volume);
    localCoordGradients[2] = (vertices[3] - vertices[0]) ^ (vertices[1] - vertices[0]) / (Scalar(6.0) * volume);
    localCoordGradients[3] = (vertices[1] - vertices[0]) ^ (vertices[2] - vertices[0]) / (Scalar(6.0) * volume);

    IndexType coordPowers[4];
    GetCoordPowers(functionIndex, coordPowers);

    Scalar localCoords[4];
    GetLocalCoords(point, localCoords);

    Scalar mult(1.0);

    for (IndexType coordIndex = 0; coordIndex < 4; coordIndex++)
    {
      for (IndexType power = 0; power < coordPowers[coordIndex]; power++)
      {
        mult *= Scalar(power + 1) / Scalar(Order);
      }
    }

    Vector2 gradient(0, 0);
    for (IndexType coordIndex0 = 0; coordIndex0 < 4; coordIndex0++)
    {
      for (IndexType power0 = 0; power0 < coordPowers[coordIndex0]; power0++)
      {
        Scalar term = mult;
        for (IndexType coordIndex1 = 0; coordIndex1 < 4; coordIndex1++)
        {
          for (IndexType power1 = 0; power1 < coordPowers[coordIndex1]; power1++)
          {
            if (coordIndex0 != coordIndex1 || power0 != power1)
            {
              term *= localCoords[coordIndex1] - Scalar(power1) / Scalar(Order);
            }
          }
        }
        gradient += localCoordGradients[coordIndex0] * term;
      }
    }
    return gradient;
  }

public:

  void GetLocalCoords(const Vector& point, Scalar *localCoords) const
  {
    Vector vertices[4];
    vertices[0] = Vector(0, 0, 0);
    vertices[1] = Vector(Scalar(1), 0, 0);
    vertices[2] = Vector(0, Scalar(1), 0);
    vertices[3] = Vector(0, 0, Scalar(1));

    Scalar volume = MixedProduct(vertices[1] - vertices[0], vertices[2] - vertices[0], vertices[3] - vertices[0]) / Scalar(6.0);
    Scalar invVolume = Scalar(1.0) / volume;

    Vector localCoordGradients[4];

    localCoordGradients[0] = (vertices[3] - vertices[1]) ^ (vertices[2] - vertices[1]) / (Scalar(6.0) * volume);
    localCoordGradients[1] = (vertices[2] - vertices[0]) ^ (vertices[3] - vertices[0]) / (Scalar(6.0) * volume);
    localCoordGradients[2] = (vertices[3] - vertices[0]) ^ (vertices[1] - vertices[0]) / (Scalar(6.0) * volume);
    localCoordGradients[3] = (vertices[1] - vertices[0]) ^ (vertices[2] - vertices[0]) / (Scalar(6.0) * volume);

    localCoords[0] = localCoordGradients[0] * (point - vertices[1]);
    localCoords[1] = localCoordGradients[1] * (point - vertices[2]);
    localCoords[2] = localCoordGradients[2] * (point - vertices[3]);
    localCoords[3] = localCoordGradients[3] * (point - vertices[0]);
  }

  Scalar GetPolynomeValue(IndexType *pows, Scalar *localCoords) const
  {
    Scalar res(1.0);

    for (IndexType coordIndex = 0; coordIndex < 4; coordIndex++)
    {
      for (IndexType power = 0; power < pows[coordIndex]; power++)
      {
        res /= Scalar(power + 1) / Scalar(Order);
      }
    }

    Scalar denom = res;

    for (IndexType coordIndex = 0; coordIndex < 4; coordIndex++)
    {
      for (IndexType power = 0; power < pows[coordIndex]; power++)
      {
        res *= localCoords[coordIndex] - Scalar(power) / Scalar(Order);
      }
    }

    Scalar baseRes = res;

    return res;
  }

  void ComputeCoordPowers()
  {
    IndexType functionIndex = 0;
    for (IndexType i = 0; i <= Order; i++)
    {
      for (IndexType j = 0; j <= Order - i; j++)
      {
        for (IndexType k = 0; k <= Order - i - j; k++)
        {
          IndexType l = Order - i - j - k;
          functionPows[functionIndex][0] = i;
          functionPows[functionIndex][1] = j;
          functionPows[functionIndex][2] = k;
          functionPows[functionIndex][3] = l;
          functionIndex++;
        }
      }
    }
  }

  void GetCoordPowers(IndexType functionIndex, IndexType *pows) const
  {
    for (IndexType i = 0; i < 4; i++)
    {
      pows[i] = functionPows[functionIndex][i];
    }
  }
  IndexType functionPows[functionsCount][4];
};
