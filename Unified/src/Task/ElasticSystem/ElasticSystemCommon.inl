template <typename Space>
void ElasticSystemCommon<Space>::ValueTypeCommon::SetZeroValues()
{
  std::fill_n(values, dimsCount, Scalar(0.0));
}

template <>
Space2::Vector ElasticSystemCommon<Space2>::ValueTypeCommon::GetVelocity() const
{
  return Space2::Vector(values[3], values[4]);
}

template <>
Space3::Vector ElasticSystemCommon<Space3>::ValueTypeCommon::GetVelocity() const
{
  return Space3::Vector(values[6], values[7], values[8]);
}

template <>
void ElasticSystemCommon<Space2>::ValueTypeCommon::SetVelocity(const Vector& velocity)
{
  values[3] = velocity.x;
  values[4] = velocity.y;
}

template <>
void ElasticSystemCommon<Space3>::ValueTypeCommon::SetVelocity(const Vector& velocity)
{
  values[6] = velocity.x;
  values[7] = velocity.y;
  values[8] = velocity.z;
}

template <>
void ElasticSystemCommon<Space2>::ValueTypeCommon::SetTension(const Tensor& tension)
{
  values[0] = tension.xx;
  values[1] = tension.yy;
  values[2] = tension.xy;
}

template <>
void ElasticSystemCommon<Space3>::ValueTypeCommon::SetTension(const Tensor& tension)
{
  values[0] = tension.xx;
  values[1] = tension.yy;
  values[2] = tension.zz;
  values[3] = tension.xy;
  values[4] = tension.yz;
  values[5] = tension.xz;
}

template <typename Space>
typename Space::Scalar ElasticSystemCommon<Space>::ValueTypeCommon::GetXX() const
{
  return values[0];
}

template <typename Space>
typename Space::Scalar ElasticSystemCommon<Space>::ValueTypeCommon::GetYY() const
{
  return values[1];
}

template <>
Space2::Scalar ElasticSystemCommon<Space2>::ValueTypeCommon::GetXY() const
{
  return values[2];
}

template <>
Space3::Scalar ElasticSystemCommon<Space3>::ValueTypeCommon::GetXY() const
{
  return values[3];
}

template <>
void ElasticSystemCommon<Space2>::ValueTypeCommon::MakeDimension(
  Scalar tensionDimensionlessMult, Scalar velocityDimensionlessMult)
{
  for (IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
  {
    switch (valueIndex)
    {
      case 0: case 1: case 2:
        values[valueIndex] *= tensionDimensionlessMult;
      break;
      case 3: case 4:
        values[valueIndex] *= velocityDimensionlessMult;
      break;
      default:
        assert(0);
      break;
    }
  }
}

template <>
void ElasticSystemCommon<Space3>::ValueTypeCommon::MakeDimension(
  Scalar tensionDimensionlessMult, Scalar velocityDimensionlessMult)
{
  for (IndexType valueIndex = 0; valueIndex < dimsCount; ++valueIndex)
  {
    switch (valueIndex)
    {
      case 0: case 1: case 2: case 3: case 4: case 5:
        values[valueIndex] *= tensionDimensionlessMult;
      break;
      case 6: case 7: case 8:
        values[valueIndex] *= velocityDimensionlessMult;
      break;
      default:
        assert(0);
      break;
    }
  }
}

template <>
void ElasticSystemCommon<Space2>::BuildXMatrix(const MediumParameters& mediumParameters, Scalar* xMatrix)
{
  std::fill_n(xMatrix, dimsCount * dimsCount, Scalar(0.0));

  xMatrix[3 ] = -(mediumParameters.lambda + Scalar(2.0) * mediumParameters.mju);
  xMatrix[8 ] = -mediumParameters.lambda;
  xMatrix[14] = -mediumParameters.mju;
  xMatrix[15] = -mediumParameters.invRho;
  xMatrix[22] = -mediumParameters.invRho;

  for (int row = 0; row < dimsCount; ++row)
  {
    // diagonal elements
    xMatrix[row + row * dimsCount] += mediumParameters.flowVelocity.x;
  }
}

template <>
void ElasticSystemCommon<Space3>::BuildXMatrix(const MediumParameters& mediumParameters, Scalar* xMatrix)
{
  std::fill_n(xMatrix, dimsCount * dimsCount, Scalar(0.0));

  xMatrix[ 6] = -(mediumParameters.lambda + Scalar(2.0) * mediumParameters.mju);
  xMatrix[15] = -mediumParameters.lambda;
  xMatrix[24] = -mediumParameters.lambda;
  xMatrix[34] = -mediumParameters.mju;
  xMatrix[53] = -mediumParameters.mju;
  xMatrix[54] = -mediumParameters.invRho;
  xMatrix[66] = -mediumParameters.invRho;
  xMatrix[77] = -mediumParameters.invRho;

  for (int row = 0; row < dimsCount; ++row)
  {
    // diagonal elements
    xMatrix[row + row * dimsCount] += mediumParameters.flowVelocity.x;
  }
}

template <>
void ElasticSystemCommon<Space2>::BuildYMatrix(const MediumParameters& mediumParameters, Scalar* yMatrix)
{
  std::fill_n(yMatrix, dimsCount * dimsCount, Scalar(0.0));

  yMatrix[4 ] = -mediumParameters.lambda;
  yMatrix[9 ] = -(mediumParameters.lambda + Scalar(2.0) * mediumParameters.mju);
  yMatrix[13] = -mediumParameters.mju;
  yMatrix[17] = -mediumParameters.invRho;
  yMatrix[21] = -mediumParameters.invRho;

  for (int row = 0; row < dimsCount; ++row)
  {
    // diagonal elements
    yMatrix[row + row * dimsCount] += mediumParameters.flowVelocity.y;
  }
}

template <>
void ElasticSystemCommon<Space3>::BuildYMatrix(const MediumParameters& mediumParameters, Scalar* yMatrix)
{
  std::fill_n(yMatrix, dimsCount * dimsCount, Scalar(0.0));

  yMatrix[ 7] = -mediumParameters.lambda;
  yMatrix[16] = -(mediumParameters.lambda + Scalar(2.0) * mediumParameters.mju);
  yMatrix[25] = -mediumParameters.lambda;
  yMatrix[33] = -mediumParameters.mju;
  yMatrix[44] = -mediumParameters.mju;
  yMatrix[57] = -mediumParameters.invRho;
  yMatrix[64] = -mediumParameters.invRho;
  yMatrix[76] = -mediumParameters.invRho;

  for (int row = 0; row < dimsCount; ++row)
  {
    // diagonal elements
    yMatrix[row + row * dimsCount] += mediumParameters.flowVelocity.y;
  }
}

template <>
void ElasticSystemCommon<Space2>::BuildRMatrix(
  const MediumParameters& interiorMediumParameters, 
  const MediumParameters& exteriorMediumParameters, 
  Scalar *rMatrix)
{
  std::fill_n(rMatrix, dimsCount * dimsCount, Scalar(0.0));
  rMatrix[0]  = interiorMediumParameters.lambda + 2 * interiorMediumParameters.mju;
  rMatrix[5]  = interiorMediumParameters.lambda;
  rMatrix[15] = interiorMediumParameters.GetPSpeed();

  rMatrix[11] = interiorMediumParameters.mju;
  rMatrix[21] = interiorMediumParameters.GetSSpeed();

  rMatrix[7]  = 1;

  rMatrix[13] = exteriorMediumParameters.mju;
  rMatrix[23] = -exteriorMediumParameters.GetSSpeed();

  rMatrix[4]  = exteriorMediumParameters.lambda + 2 * exteriorMediumParameters.mju;
  rMatrix[9]  = exteriorMediumParameters.lambda;
  rMatrix[19] = -exteriorMediumParameters.GetPSpeed();
}

template <>
void ElasticSystemCommon<Space3>::BuildRMatrix(
  const MediumParameters& interiorMediumParameters, 
  const MediumParameters& exteriorMediumParameters, 
  Scalar *rMatrix)
{
  std::fill(rMatrix, rMatrix + dimsCount * dimsCount, Scalar(0.0));

  rMatrix[0]  = interiorMediumParameters.lambda + 2 * interiorMediumParameters.mju;
  rMatrix[9]  = interiorMediumParameters.lambda;
  rMatrix[18] = interiorMediumParameters.lambda;
  rMatrix[54] = interiorMediumParameters.GetPSpeed();

  rMatrix[28] = interiorMediumParameters.mju;
  rMatrix[64] = interiorMediumParameters.GetSSpeed();

  rMatrix[47] = interiorMediumParameters.mju;
  rMatrix[74] = interiorMediumParameters.GetSSpeed();

  rMatrix[39] = Scalar(1.0);
  rMatrix[13] = Scalar(1.0);
  rMatrix[23] = Scalar(1.0);

  rMatrix[51] =  exteriorMediumParameters.mju;
  rMatrix[78] = -exteriorMediumParameters.GetSSpeed();

  rMatrix[34] = exteriorMediumParameters.mju;
  rMatrix[70] = -exteriorMediumParameters.GetSSpeed();

  rMatrix[8] = exteriorMediumParameters.lambda + 2 * exteriorMediumParameters.mju;
  rMatrix[17] = exteriorMediumParameters.lambda;
  rMatrix[26] = exteriorMediumParameters.lambda;
  rMatrix[62] = -exteriorMediumParameters.GetPSpeed();
}

template <>
void ElasticSystemCommon<Space2>::BuildXnAuxMatrix(
  const MediumParameters& interiorMediumParameters, 
  const MediumParameters& exteriorMediumParameters,
  Scalar* xnAuxMatrix)
{
/*  
  Scalar rMatrix[dimsCount * dimsCount];
  BuildRMatrix(interiorMediumParameters, exteriorMediumParameters, rMatrix);

  Scalar invRMatrix[dimsCount * dimsCount];
  MatrixInverse(rMatrix, invRMatrix, dimsCount);

  Scalar lambdaMatrix[dimsCount * dimsCount];
  std::fill(lambdaMatrix, lambdaMatrix + dimsCount * dimsCount, Scalar(0.0));
  lambdaMatrix[0]  = interiorMediumParameters.GetPSpeed();
  lambdaMatrix[6]  = interiorMediumParameters.GetSSpeed();

  Scalar tmpMatrix[dimsCount * dimsCount];
  MatrixMulMatrix(rMatrix, lambdaMatrix, tmpMatrix, dimsCount, dimsCount, dimsCount);
  MatrixMulMatrix(tmpMatrix, invRMatrix, xnAuxMatrix, dimsCount, dimsCount, dimsCount);
*/

  std::fill_n(xnAuxMatrix, dimsCount * dimsCount, Scalar(0.0));
  Scalar cpInt = interiorMediumParameters.GetPSpeed();
  Scalar cpExt = exteriorMediumParameters.GetPSpeed();
  Scalar csInt = interiorMediumParameters.GetSSpeed();
  Scalar csExt = exteriorMediumParameters.GetSSpeed();

  Scalar rhoInt = 1 / interiorMediumParameters.invRho;
  Scalar rhoExt = 1 / exteriorMediumParameters.invRho;

  Scalar zpInt = cpInt * rhoInt;
  Scalar zsInt = csInt * rhoInt;
  Scalar zpExt = cpExt * rhoExt;
  Scalar zsExt = csExt * rhoExt;

  Scalar coeffP = 1 / (zpInt + zpExt);
  Scalar coeffS = Scalar(1.0) / (rhoInt + rhoExt);
  
  if (zsInt + zsExt > std::numeric_limits<Scalar>::epsilon())
  {
    coeffS = csInt / (zsInt + zsExt);
  }

  xnAuxMatrix[0]  = coeffP * zpInt * cpInt;
  xnAuxMatrix[3]  = coeffP * zpInt * zpExt * cpInt;
  xnAuxMatrix[5]  = coeffP * interiorMediumParameters.lambda;
  xnAuxMatrix[8]  = coeffP * interiorMediumParameters.lambda * zpExt;
  xnAuxMatrix[12] = coeffS * zsInt;
  xnAuxMatrix[14] = coeffS * zsInt * zsExt;
  xnAuxMatrix[15] = coeffP * cpInt;
  xnAuxMatrix[18] = coeffP * cpInt * zpExt;
  xnAuxMatrix[22] = coeffS;
  xnAuxMatrix[24] = coeffS * zsExt;
}

template <>
void ElasticSystemCommon<Space3>::BuildXnAuxMatrix(
  const MediumParameters& interiorMediumParameters, 
  const MediumParameters& exteriorMediumParameters,
  Scalar* xnAuxMatrix)
{
/*
  Scalar rMatrix[dimsCount * dimsCount];
  BuildRMatrix(interiorMediumParameters, exteriorMediumParameters, rMatrix);

  Scalar invRMatrix[dimsCount * dimsCount];
  MatrixInverse(rMatrix, invRMatrix, dimsCount);

  Scalar lambdaMatrix[dimsCount * dimsCount];
  std::fill_n(lambdaMatrix, dimsCount * dimsCount, Scalar(0.0));
  lambdaMatrix[0]  = interiorMediumParameters.GetPSpeed();
  lambdaMatrix[10] = interiorMediumParameters.GetSSpeed();
  lambdaMatrix[20] = interiorMediumParameters.GetSSpeed();

  Scalar tmpMatrix[dimsCount * dimsCount];
  MatrixMulMatrix(rMatrix, lambdaMatrix, tmpMatrix, dimsCount, dimsCount, dimsCount);
  MatrixMulMatrix(tmpMatrix, invRMatrix, xnAuxMatrix, dimsCount, dimsCount, dimsCount);
*/

  std::fill_n(xnAuxMatrix, dimsCount * dimsCount, Scalar(0.0));
  Scalar cpInt = interiorMediumParameters.GetPSpeed();
  Scalar cpExt = exteriorMediumParameters.GetPSpeed();
  Scalar csInt = interiorMediumParameters.GetSSpeed();
  Scalar csExt = exteriorMediumParameters.GetSSpeed();

  Scalar rhoInt = 1 / interiorMediumParameters.invRho;
  Scalar rhoExt = 1 / exteriorMediumParameters.invRho;

  Scalar zpInt = cpInt * rhoInt;
  Scalar zsInt = csInt * rhoInt;
  Scalar zpExt = cpExt * rhoExt;
  Scalar zsExt = csExt * rhoExt;

  Scalar coeffP = 1 / (zpInt + zpExt);
  Scalar coeffS = Scalar(1.0) / (rhoInt + rhoExt);

  if (zsInt + zsExt > std::numeric_limits<Scalar>::epsilon())
  {
    coeffS = csInt / (zsInt + zsExt);
  }

  xnAuxMatrix[ 0] = coeffP * zpInt * cpInt;
  xnAuxMatrix[ 9] = 
  xnAuxMatrix[18] = coeffP * interiorMediumParameters.lambda;
  xnAuxMatrix[54] = coeffP * cpInt;

  xnAuxMatrix[30] =
  xnAuxMatrix[50] = coeffS * zsInt;

  xnAuxMatrix[66] = 
  xnAuxMatrix[77] = coeffS;

  xnAuxMatrix[ 6] = coeffP * zpInt * zpExt * cpInt;
  xnAuxMatrix[15] = 
  xnAuxMatrix[24] = coeffP * zpExt * interiorMediumParameters.lambda;

  xnAuxMatrix[60] = coeffP * zpExt * cpInt;
  xnAuxMatrix[34] = 
  xnAuxMatrix[53] = coeffS * zsInt * zsExt;
  xnAuxMatrix[70] = 
  xnAuxMatrix[80] = coeffS * zsExt;
}

template <>
void ElasticSystemCommon<Space2>::BuildXnInteriorMatrix(
  const MediumParameters& interiorMediumParameters, 
  const MediumParameters& exteriorMediumParameters,
  const Vector& edgeNormal, Scalar* xnInteriorMatrix)
{
  BuildXnAuxMatrix(interiorMediumParameters, exteriorMediumParameters, xnInteriorMatrix);
  
  xnInteriorMatrix[3 ] -= interiorMediumParameters.lambda + Scalar(2.0) * interiorMediumParameters.mju;
  xnInteriorMatrix[8 ] -= interiorMediumParameters.lambda;
  xnInteriorMatrix[14] -= interiorMediumParameters.mju;
  xnInteriorMatrix[15] -= interiorMediumParameters.invRho;
  xnInteriorMatrix[22] -= interiorMediumParameters.invRho;

  for (int row = 0; row < dimsCount; ++row)
  {
    // diagonal elements
    Scalar un = interiorMediumParameters.flowVelocity * edgeNormal;
    xnInteriorMatrix[row + row * dimsCount] += (un + fabs(un)) * Scalar(0.5);
  }
}

template <>
void ElasticSystemCommon<Space3>::BuildXnInteriorMatrix(
  const MediumParameters& interiorMediumParameters, 
  const MediumParameters& exteriorMediumParameters,
  const Vector& faceNormal, Scalar* xnInteriorMatrix)
{
  BuildXnAuxMatrix(interiorMediumParameters, exteriorMediumParameters, xnInteriorMatrix);

  xnInteriorMatrix[ 6] -= interiorMediumParameters.lambda + Scalar(2.0) * interiorMediumParameters.mju;
  xnInteriorMatrix[15] -= interiorMediumParameters.lambda;
  xnInteriorMatrix[24] -= interiorMediumParameters.lambda;
  xnInteriorMatrix[34] -= interiorMediumParameters.mju;
  xnInteriorMatrix[53] -= interiorMediumParameters.mju;
  xnInteriorMatrix[54] -= interiorMediumParameters.invRho;
  xnInteriorMatrix[66] -= interiorMediumParameters.invRho;
  xnInteriorMatrix[77] -= interiorMediumParameters.invRho;

  for (int row = 0; row < dimsCount; ++row)
  {
    // diagonal elements
    Scalar un = interiorMediumParameters.flowVelocity * faceNormal;
    xnInteriorMatrix[row + row * dimsCount] += (un + fabs(un)) * Scalar(0.5);
  }
}

template <typename Space>
void ElasticSystemCommon<Space>::BuildXnExteriorMatrix(
  const MediumParameters& interiorMediumParameters, 
  const MediumParameters& exteriorMediumParameters,
  const Vector& normal, Scalar* xnExteriorMatrix)
{
  BuildXnAuxMatrix(interiorMediumParameters, exteriorMediumParameters, xnExteriorMatrix);
  MatrixMulScalar(xnExteriorMatrix, Scalar(-1.0), dimsCount, dimsCount);

  for (int row = 0; row < dimsCount; ++row)
  {
    // diagonal elements
    Scalar un = interiorMediumParameters.flowVelocity * normal;
    xnExteriorMatrix[row + row * dimsCount] += (un - fabs(un)) * Scalar(0.5);
  }
}

template <>
void ElasticSystemCommon<Space2>::BuildBoundaryMatrix(IndexType interactionType, Scalar* boundaryMatrix)
{
  std::fill_n(boundaryMatrix, dimsCount * dimsCount, Scalar(0.0));

  IndexType boundaryType = boundaryDescriptions[interactionType].type;
  Scalar reflectionCoeff = boundaryDescriptions[interactionType].reflectionCoeff;

  switch(boundaryType)
  {
    case BoundaryConditions::Absorb:
    {
      boundaryMatrix[0 ] = 0;
      boundaryMatrix[6 ] = 0;
      boundaryMatrix[12] = 0;
      boundaryMatrix[18] = 0;
      boundaryMatrix[24] = 0;
    } break;
    case BoundaryConditions::Free:
    {
      boundaryMatrix[ 0] = -Scalar(reflectionCoeff);
      boundaryMatrix[ 6] =  Scalar(reflectionCoeff);
      boundaryMatrix[12] = -Scalar(reflectionCoeff);
      boundaryMatrix[18] =  Scalar(reflectionCoeff);
      boundaryMatrix[24] =  Scalar(reflectionCoeff);
    } break;
    case BoundaryConditions::Fixed:
    {
      boundaryMatrix[ 0] =  Scalar(reflectionCoeff);
      boundaryMatrix[ 6] =  Scalar(reflectionCoeff);
      boundaryMatrix[12] =  Scalar(reflectionCoeff);
      boundaryMatrix[18] = -Scalar(reflectionCoeff);
      boundaryMatrix[24] = -Scalar(reflectionCoeff);
    } break;
    case BoundaryConditions::Symmetry:
    {
      boundaryMatrix[0 ] =  Scalar(1.0); 
      boundaryMatrix[6 ] =  Scalar(1.0);
      boundaryMatrix[12] = -Scalar(1.0); // sigma.xy
      boundaryMatrix[18] = -Scalar(1.0); // v.x
      boundaryMatrix[24] =  Scalar(1.0); 
    } break;
    case BoundaryConditions::AntiSymmetry:
    {
      boundaryMatrix[ 0] = -Scalar(1.0); // sigma.xx
      boundaryMatrix[ 6] = Scalar(1.0);
      boundaryMatrix[12] = Scalar(1.0); 
      boundaryMatrix[18] = Scalar(1.0); 
      boundaryMatrix[24] = -Scalar(1.0); // v.y
    } break;
    default:
      assert(0);
    break;
  }
}

template <>
void ElasticSystemCommon<Space3>::BuildBoundaryMatrix(IndexType interactionType, Scalar *boundaryMatrix)
{
  std::fill_n(boundaryMatrix, dimsCount * dimsCount, Scalar(0.0));

  IndexType boundaryType = boundaryDescriptions[interactionType].type;
  Scalar reflectionCoeff = boundaryDescriptions[interactionType].reflectionCoeff;

  switch (boundaryType)
  {
    case BoundaryConditions::Absorb:
    {
      boundaryMatrix[ 0] = 0;
      boundaryMatrix[10] = 0;
      boundaryMatrix[20] = 0;
      boundaryMatrix[30] = 0;
      boundaryMatrix[40] = 0;
      boundaryMatrix[50] = 0;
      boundaryMatrix[60] = 0;
      boundaryMatrix[70] = 0;
      boundaryMatrix[80] = 0;
    } break;
    case BoundaryConditions::Free:
    {
      boundaryMatrix[ 0] = -Scalar(reflectionCoeff);
      boundaryMatrix[10] =  Scalar(reflectionCoeff);
      boundaryMatrix[20] =  Scalar(reflectionCoeff);
      boundaryMatrix[30] = -Scalar(reflectionCoeff);
      boundaryMatrix[40] =  Scalar(reflectionCoeff);
      boundaryMatrix[50] = -Scalar(reflectionCoeff);
      boundaryMatrix[60] =  Scalar(reflectionCoeff);
      boundaryMatrix[70] =  Scalar(reflectionCoeff);
      boundaryMatrix[80] =  Scalar(reflectionCoeff);
    } break;
    case BoundaryConditions::Fixed:
    {
      boundaryMatrix[ 0] =  Scalar(reflectionCoeff);
      boundaryMatrix[10] =  Scalar(reflectionCoeff);
      boundaryMatrix[20] =  Scalar(reflectionCoeff);
      boundaryMatrix[30] =  Scalar(reflectionCoeff);
      boundaryMatrix[40] =  Scalar(reflectionCoeff);
      boundaryMatrix[50] =  Scalar(reflectionCoeff);
      boundaryMatrix[60] = -Scalar(reflectionCoeff);
      boundaryMatrix[70] = -Scalar(reflectionCoeff);
      boundaryMatrix[80] = -Scalar(reflectionCoeff);
    } break;
    case BoundaryConditions::Symmetry:
    {
      boundaryMatrix[ 0] =  Scalar(reflectionCoeff);
      boundaryMatrix[10] =  Scalar(reflectionCoeff);
      boundaryMatrix[20] =  Scalar(reflectionCoeff);
      boundaryMatrix[30] = -Scalar(reflectionCoeff); // sigma xy
      boundaryMatrix[40] = -Scalar(reflectionCoeff); // sigma yz
      boundaryMatrix[50] = -Scalar(reflectionCoeff); // sigma xz
      boundaryMatrix[60] = -Scalar(reflectionCoeff); // v.x
      boundaryMatrix[70] =  Scalar(reflectionCoeff);
      boundaryMatrix[80] =  Scalar(reflectionCoeff);
    } break;
    case BoundaryConditions::AntiSymmetry:
    {
      boundaryMatrix[ 0] = -Scalar(reflectionCoeff); // sigma.xx
      boundaryMatrix[10] = Scalar(reflectionCoeff);
      boundaryMatrix[20] = Scalar(reflectionCoeff);
      boundaryMatrix[30] = Scalar(reflectionCoeff); 
      boundaryMatrix[40] = Scalar(reflectionCoeff); 
      boundaryMatrix[50] = Scalar(reflectionCoeff); 
      boundaryMatrix[60] = Scalar(reflectionCoeff); 
      boundaryMatrix[70] = -Scalar(reflectionCoeff); // v.y
      boundaryMatrix[80] = -Scalar(reflectionCoeff); // v.z
    } break;
    default:
      assert(0);
    break;
  }
}

template <>
void ElasticSystemCommon<Space2>::BuildContactMatrices(IndexType interactionType, 
  Scalar* leftContactMatrix, Scalar* rightContactMatrix)
{
  std::fill_n(leftContactMatrix,  IndexType(dimsCount * dimsCount), Scalar(0.0));
  std::fill_n(rightContactMatrix, IndexType(dimsCount * dimsCount), Scalar(0.0));

  IndexType contactType = contactDescriptions[interactionType].type;

  switch(contactType)
  {
    case ContactConditions::Glue:
    {
      leftContactMatrix[0 ] = 0;
      leftContactMatrix[6 ] = 0;
      leftContactMatrix[12] = 0;
      leftContactMatrix[18] = 0;
      leftContactMatrix[24] = 0;

      rightContactMatrix[0 ] = 1;
      rightContactMatrix[6 ] = 1;
      rightContactMatrix[12] = 1;
      rightContactMatrix[18] = 1;
      rightContactMatrix[24] = 1;
    } break;
    case ContactConditions::Glide:
    {
      // it`s correct when materials on the both sides of contact are the same 
      leftContactMatrix[0 ] = 0;
      leftContactMatrix[6 ] = 1;
      leftContactMatrix[12] = -1;
      leftContactMatrix[18] = 0;
      leftContactMatrix[24] = 1;

      rightContactMatrix[0 ] = 1;
      rightContactMatrix[6 ] = 0;
      rightContactMatrix[12] = 0;
      rightContactMatrix[18] = 1;
      rightContactMatrix[24] = 0;
    } break;
    case ContactConditions::Friction:
    {
      leftContactMatrix[0 ] = 0;
      leftContactMatrix[6 ] = 0;
      leftContactMatrix[12] = 0;
      leftContactMatrix[18] = 0;
      leftContactMatrix[24] = 0;

      rightContactMatrix[0 ] = 1;
      rightContactMatrix[6 ] = 1;
      rightContactMatrix[12] = 1;
      rightContactMatrix[18] = 1;
      rightContactMatrix[24] = 1;
    } break;
  }
}

template <>
void ElasticSystemCommon<Space3>::BuildContactMatrices(IndexType interactionType,
  Scalar* leftContactMatrix, Scalar* rightContactMatrix)
{
  std::fill_n(leftContactMatrix, IndexType(dimsCount * dimsCount), Scalar(0.0));
  std::fill_n(rightContactMatrix, IndexType(dimsCount * dimsCount), Scalar(0.0));

  IndexType contactType = contactDescriptions[interactionType].type;

  switch (contactType)
  {
    case ContactConditions::Glue:
    {
      leftContactMatrix[0 ] = 0;
      leftContactMatrix[10] = 0;
      leftContactMatrix[20] = 0;
      leftContactMatrix[30] = 0;
      leftContactMatrix[40] = 0;
      leftContactMatrix[50] = 0;
      leftContactMatrix[60] = 0;
      leftContactMatrix[70] = 0;
      leftContactMatrix[80] = 0;

      rightContactMatrix[0 ] = 1;
      rightContactMatrix[10] = 1;
      rightContactMatrix[20] = 1;
      rightContactMatrix[30] = 1;
      rightContactMatrix[40] = 1;
      rightContactMatrix[50] = 1;
      rightContactMatrix[60] = 1;
      rightContactMatrix[70] = 1;
      rightContactMatrix[80] = 1;
    } break;
    case ContactConditions::Glide:
    {
      leftContactMatrix[ 0] = 0;
      leftContactMatrix[10] = 1;
      leftContactMatrix[20] = 1;
      leftContactMatrix[30] = -1;
      leftContactMatrix[40] = -1;
      leftContactMatrix[50] = -1;
      leftContactMatrix[60] = 0;
      leftContactMatrix[70] = 1;
      leftContactMatrix[80] = 1;

      rightContactMatrix[ 0] = 1;
      rightContactMatrix[10] = 0;
      rightContactMatrix[20] = 0;
      rightContactMatrix[30] = 0;
      rightContactMatrix[40] = 0;
      rightContactMatrix[50] = 0;
      rightContactMatrix[60] = 1;
      rightContactMatrix[70] = 0;
      rightContactMatrix[80] = 0;
    } break;
    case ContactConditions::Friction:
    {
      leftContactMatrix[ 0] = 0;
      leftContactMatrix[10] = 0;
      leftContactMatrix[20] = 0;
      leftContactMatrix[30] = 0;
      leftContactMatrix[40] = 0;
      leftContactMatrix[50] = 0;
      leftContactMatrix[60] = 0;
      leftContactMatrix[70] = 0;
      leftContactMatrix[80] = 0;

      rightContactMatrix[ 0] = 1;
      rightContactMatrix[10] = 1;
      rightContactMatrix[20] = 1;
      rightContactMatrix[30] = 1;
      rightContactMatrix[40] = 1;
      rightContactMatrix[50] = 1;
      rightContactMatrix[60] = 1;
      rightContactMatrix[70] = 1;
      rightContactMatrix[80] = 1;
    } break;
    default:
      assert(0);
    break;
  }
}

template <>
void ElasticSystemCommon<Space2>::CorrectContact(Scalar* values, 
  Vector contactNormal, IndexType dynamicInteractionType)
{
  IndexType contactType = contactDescriptions[dynamicInteractionType].type;
  contactNormal.Normalize();
  switch (contactType)
  {
    case ContactConditions::Glue:
    {
    } break;
    case ContactConditions::Glide:
    {
    } break;
    case ContactConditions::Friction:
    {
      Scalar frictionCoeff;
      GetFrictionContactInfo(dynamicInteractionType, frictionCoeff);
      // transform into local coordinate system
      Scalar sigmaNormal = values[0] * Sqr(contactNormal.x) + 
                            values[1] * Sqr(contactNormal.y) + 
                            values[2] * 2 * contactNormal.x * contactNormal.y;

      if (sigmaNormal > 0) // sigma_xx > 0
      {
        Scalar sigmaTangent = values[0] * Sqr(contactNormal.y) + 
                              values[1] * Sqr(contactNormal.x) - 
                              values[2] * 2 * contactNormal.x * contactNormal.y;

        Scalar tau          = -(values[0] - values[1]) * contactNormal.x * contactNormal.y + 
                                            values[2]  * (Sqr(contactNormal.x) - Sqr(contactNormal.y));

        tau = Sgn(tau) * std::min(fabs(tau), frictionCoeff * sigmaNormal); // sigma_xy correction

        values[0] = sigmaNormal * Sqr(contactNormal.x) + sigmaTangent * Sqr(contactNormal.y) - tau * 2 * contactNormal.x * contactNormal.y;
        values[1] = sigmaNormal * Sqr(contactNormal.y) + sigmaTangent * Sqr(contactNormal.x) + tau * 2 * contactNormal.x * contactNormal.y;
        values[2] = (sigmaNormal - sigmaTangent) * contactNormal.x * contactNormal.y + tau * (Sqr(contactNormal.x) - Sqr(contactNormal.y));
      }
    } break;
    default:
      assert(0);
    break;
  }
}

template <>
void ElasticSystemCommon<Space3>::CorrectContact(Scalar* value, 
  Vector contactNormal, IndexType dynamicInteractionType)
{
  // TODO
}

template<typename Space>
void ElasticSystemCommon<Space>::SetAbsorbBoundary(IndexType interactionType)
{
  BoundaryDescription newBoundary;
  newBoundary.type = BoundaryConditions::Absorb;
  newBoundary.infoIndex = -1;
  SetBoundaryDescription(interactionType, newBoundary);
}

template<typename Space>
void ElasticSystemCommon<Space>::SetSymmetryBoundary(IndexType interactionType)
{
  BoundaryDescription newBoundary;
  newBoundary.type = BoundaryConditions::Symmetry;
  newBoundary.infoIndex = -1;
  SetBoundaryDescription(interactionType, newBoundary);
}

template<typename Space>
void ElasticSystemCommon<Space>::SetAntiSymmetryBoundary(IndexType interactionType)
{
  BoundaryDescription newBoundary;
  newBoundary.type = BoundaryConditions::AntiSymmetry;
  newBoundary.infoIndex = -1;
  SetBoundaryDescription(interactionType, newBoundary);
}

template<typename Space>
void ElasticSystemCommon<Space>::SetGlueContact(IndexType interactionType)
{
  ContactDescription newContact;
  newContact.type = ContactConditions::Glue;
  newContact.infoIndex = -1;
  SetContactDescription(interactionType, newContact);
}

template<typename Space>
void ElasticSystemCommon<Space>::SetGlueContact(IndexType interactionType, 
  Scalar maxShearStress, Scalar maxLongitudinalStress, IndexType dynamicBoundaryInteractionType)
{
  ContactDescription newContact;
  newContact.type = ContactConditions::Glue;
  newContact.infoIndex = glueContactInfos.size();
  SetContactDescription(interactionType, newContact);

  glueContactInfos.push_back(
    GlueContactInfo<Space>(maxShearStress, maxLongitudinalStress, dynamicBoundaryInteractionType));
}

template<typename Space>
void ElasticSystemCommon<Space>::SetFrictionContact(IndexType interactionType, Scalar frictionCoeff, IndexType dynamicBoundaryInteractionType)
{
  ContactDescription newContact;
  newContact.type = ContactConditions::Friction;
  newContact.infoIndex = frictionContactInfos.size();
  SetContactDescription(interactionType, newContact);

  frictionContactInfos.push_back(FrictionContactInfo<Space>(frictionCoeff, dynamicBoundaryInteractionType));
}

template<typename Space>
void ElasticSystemCommon<Space>::SetFrictionContact(IndexType interactionType)
{
  ContactDescription newContact;
  newContact.type = ContactConditions::Friction;
  newContact.infoIndex = frictionContactInfos.size();
  SetContactDescription(interactionType, newContact);
}

template<typename Space>
void ElasticSystemCommon<Space>::SetGlideContact(IndexType interactionType)
{
  ContactDescription newContact;
  newContact.type = ContactConditions::Glide;
  newContact.infoIndex = -1;
  SetContactDescription(interactionType, newContact);
}

template<typename Space>
typename Space::IndexType ElasticSystemCommon<Space>::GetContactDynamicBoundaryType(IndexType contactInteractionType )
{
  if(contactDescriptions[contactInteractionType].type == ContactConditions::Glue &&
     contactDescriptions[contactInteractionType].infoIndex != IndexType(-1))
  {
    IndexType glueInfoIndex = contactDescriptions[contactInteractionType].infoIndex;
    return glueContactInfos[glueInfoIndex].dynamicBoundaryInteractionType;
  }
  return IndexType(-1);
}

template<typename Space>
typename Space::IndexType ElasticSystemCommon<Space>::GetBoundaryDynamicContactType(IndexType boundaryInteractionType)
{
  if(boundaryDescriptions[boundaryInteractionType].type == BoundaryConditions::Free &&
     boundaryDescriptions[boundaryInteractionType].infoIndex != IndexType(-1))
  {
    IndexType freeInfoIndex = boundaryDescriptions[boundaryInteractionType].infoIndex;
    return freeBoundaryInfos[freeInfoIndex].dynamicContactInteractionType;
  }
  return IndexType(-1);
}

template<typename Space>
void ElasticSystemCommon<Space>::GetFrictionContactInfo(IndexType interactionType, Scalar &frictionCoeff)
{
  if(contactDescriptions[interactionType].type == ContactConditions::Friction &&
     contactDescriptions[interactionType].infoIndex != IndexType(-1))
  {
    frictionCoeff = frictionContactInfos[contactDescriptions[interactionType].infoIndex].frictionCoeff;
    return;
  }
  assert(0);
  frictionCoeff = 0; //should never happen. hopefully.
}

template<typename Space>
void ElasticSystemCommon<Space>::GetContactCriticalInfo(IndexType interactionType, 
  Scalar &maxShearStress, Scalar &maxLongitudinalStress)
{
  assert(GetContactDynamicBoundaryType(interactionType) != IndexType(-1));
  if(contactDescriptions[interactionType].type == ContactConditions::Glue &&
     contactDescriptions[interactionType].infoIndex != IndexType(-1))
  {
    maxShearStress        = glueContactInfos[contactDescriptions[interactionType].infoIndex].maxShearStress;
    maxLongitudinalStress = glueContactInfos[contactDescriptions[interactionType].infoIndex].maxLongitudinalStress;
    return;
  }
  assert(0);
  maxShearStress = 0; //should never happen. hopefully.
  maxLongitudinalStress = 0;
}

template<typename Space>
BoundaryConditions::Types ElasticSystemCommon<Space>::GetBoundaryType(IndexType interactionType) const
{
  return boundaryDescriptions[interactionType].type;
}

template<typename Space>
ContactConditions::Types ElasticSystemCommon<Space>::GetContactType(IndexType interactionType) const
{
  return contactDescriptions[interactionType].type;
}

template<typename Space>
typename ElasticSystemCommon<Space>::SourceFunctorT* ElasticSystemCommon<Space>::GetSourceFunctor() const
{
  return sourceFunctor;
}

template<typename Space>
void ElasticSystemCommon<Space>::SetSourceFunctor(SourceFunctorT* sourceFunctor)
{
  this->sourceFunctor = sourceFunctor;
}

template<typename Space>
void ElasticSystemCommon<Space>::SetPointSources(const std::vector<PointSource<Space>* >& pointSources)
{
  this->pointSources = pointSources;
}

template<typename Space>
typename Space::Scalar ElasticSystemCommon<Space>::GetMaxWaveSpeed(const MediumParameters& mediumParameters) const
{
  return mediumParameters.GetPSpeed();
}

template<typename Space>
void ElasticSystemCommon<Space>::SetBoundaryDescription(IndexType interactionType, BoundaryDescription boundaryDescription)
{
  if(boundaryDescriptions.size() < interactionType + 1)
  {
    boundaryDescriptions.resize(interactionType + 1);
  }
  boundaryDescriptions[interactionType] = boundaryDescription;
}

template<typename Space>
void ElasticSystemCommon<Space>::SetContactDescription(IndexType interactionType, ContactDescription contactDescription)
{
  if(contactDescriptions.size() < interactionType + 1)
  {
    contactDescriptions.resize(interactionType + 1);
  }
  contactDescriptions[interactionType] = contactDescription;
}

template <typename Space>
ElasticSystemCommon<Space>::MediumParameters::MediumParameters(): 
  lambda(Scalar(2.0)), 
  mju(Scalar(1.0)), 
  invRho(Scalar(1.0)),
  destroyed(false),
  dimensionless(false)
{
}

template <typename Space>
ElasticSystemCommon<Space>::MediumParameters::MediumParameters(
  Scalar lambda, Scalar mju, Scalar invRho): 
  lambda(lambda), 
  mju(mju), 
  invRho(invRho),
  destroyed(false),
  dimensionless(false)
{
}

template <typename Space>
void ElasticSystemCommon<Space>::MediumParameters::SetFromVelocities(Scalar pSpeed, Scalar sSpeed, Scalar rho)
{
  mju = Sqr(sSpeed) * rho;
  lambda = (Sqr(pSpeed) - 2.0 * Sqr(sSpeed)) * rho;
  invRho = 1 / rho;
}

template <typename Space>
void ElasticSystemCommon<Space>::MediumParameters::SetFromLameParameters(Scalar lambda, Scalar mju, Scalar rho)
{
  this->lambda = lambda;
  this->mju    = mju;
  invRho = 1 / rho;
}

template <typename Space>
void ElasticSystemCommon<Space>::MediumParameters::SetFromYoungModuleAndNju(Scalar youngModule, Scalar nju, Scalar rho)
{
  this->lambda = nju * youngModule / (1 + nju) / (1 - 2 * nju);
  this->mju    = youngModule / Scalar(2.0) / (1 + nju);
  invRho = 1 / rho;
}

template <typename Space>
typename Space::Scalar ElasticSystemCommon<Space>::MediumParameters::GetPSpeed() const
{
  return sqrt((lambda + Scalar(2.0) * mju) * invRho);
}

template <typename Space>
typename Space::Scalar ElasticSystemCommon<Space>::MediumParameters::GetSSpeed() const
{
  return sqrt(mju * invRho);
}

template <typename Space>
typename Space::Scalar ElasticSystemCommon<Space>::MediumParameters::GetZp() const
{
  return GetPSpeed() / invRho;
}

template <typename Space>
typename Space::Scalar ElasticSystemCommon<Space>::MediumParameters::GetZs() const
{
  return GetSSpeed() / invRho;
}

template <typename Space>
void ElasticSystemCommon<Space>::MediumParameters::MakeDimensionless(
  Scalar tensionDimensionlessMult, Scalar velocityDimensionlessMult)
{
  if (!dimensionless)
  {
    dimensionless = true;
    lambda        /= tensionDimensionlessMult;
    mju           /= tensionDimensionlessMult;
    invRho        *= tensionDimensionlessMult / Sqr(velocityDimensionlessMult);
    plasticity.k0 /= tensionDimensionlessMult;
    flowVelocity  /= velocityDimensionlessMult;
  }
}

template<typename Space>
bool ElasticSystemCommon<Space>::IsProperContact(const ValueTypeCommon& value0, 
                                                 const ValueTypeCommon& value1, const Vector& contactNormal)
{
  // return true;
  Scalar eps = std::numeric_limits<Scalar>::epsilon();
  if (value1.GetForce(-contactNormal) * contactNormal - value0.GetForce(contactNormal) * contactNormal > -eps &&
    (value0.GetVelocity() - value1.GetVelocity()) * contactNormal > 0) return true;
  return false;
}

template<typename Space>
BoundaryInfoFunctor<Space>*
  ElasticSystemCommon<Space>::GetBoundaryInfoFunctor(IndexType interactionType)
{
  int infoIndex = boundaryDescriptions[interactionType].infoIndex;

  BoundaryFunctionGetter<Space>* getter = 0;

  switch(boundaryDescriptions[interactionType].type)
  {
    case BoundaryConditions::Absorb:
    {
      return 0;
    }break;
    case BoundaryConditions::Free:
    {
      assert(infoIndex != -1);
      return freeBoundaryInfos[infoIndex].boundaryFunctor;
    }break;
    case BoundaryConditions::Fixed:
    {
      assert(infoIndex != -1);
      return fixedBoundaryInfos[infoIndex].boundaryFunctor;
    }break;
    case BoundaryConditions::Symmetry:
    {
      return 0;
    }break;
    case BoundaryConditions::AntiSymmetry:
    {
      return 0;
    }break;
    default:
      assert(0);
    break;
  }
  assert(0); //unknown interation type
  return 0;
}

template<typename Space>
void ElasticSystemCommon<Space>::SetFreeBoundary(IndexType interactionType, Scalar tensionDimensionlessMult, Scalar reflectionCoeff,
  VectorFunctor<Space>* externalForceFunctor, 
  IndexType dynamicContactInteractionType)
{
  BoundaryDescription newBoundary;
  newBoundary.type = BoundaryConditions::Free;
  newBoundary.infoIndex = freeBoundaryInfos.size();
  newBoundary.reflectionCoeff = reflectionCoeff;
  SetBoundaryDescription(interactionType, newBoundary);

  FreeBoundaryInfo<Space> freeInfo;

  if (externalForceFunctor)
    freeInfo.boundaryFunctor = externalForceFunctor ? new FreeBoundaryInfoFunctor(externalForceFunctor, tensionDimensionlessMult) : 0;
  freeInfo.dynamicContactInteractionType = dynamicContactInteractionType;
  freeBoundaryInfos.push_back(freeInfo);
}

template<typename Space>
void ElasticSystemCommon<Space>::SetFixedBoundary(IndexType interactionType, Scalar velocityDimensionlessMult, Scalar reflectionCoeff,
  VectorFunctor<Space>* externalVelocityFunctor)
{
  BoundaryDescription newBoundary;
  newBoundary.type = BoundaryConditions::Fixed;
  newBoundary.infoIndex = fixedBoundaryInfos.size();
  newBoundary.reflectionCoeff = reflectionCoeff;
  SetBoundaryDescription(interactionType, newBoundary);
  
  FixedBoundaryInfo<Space> fixedInfo;
  fixedInfo.boundaryFunctor = externalVelocityFunctor ? new FixedBoundaryInfoFunctor(externalVelocityFunctor, velocityDimensionlessMult) : 0;

  fixedBoundaryInfos.push_back(fixedInfo);
}

template <>
void ElasticSystemCommon<Space2>::FreeBoundaryInfoFunctor::SetTension(
  const Vector& externalNormal, 
  const Vector& force, Scalar* values)
{
  values[0] = force.x * Sqr(externalNormal.x) - 2 * force.y * externalNormal.x * externalNormal.y;
  values[1] = force.x * Sqr(externalNormal.y) + 2 * force.y * externalNormal.x * externalNormal.y;
  values[2] = force.x * externalNormal.x * externalNormal.y + force.y * (Sqr(externalNormal.x) - Sqr(externalNormal.y));
}

template <>
void ElasticSystemCommon<Space3>::FreeBoundaryInfoFunctor::SetTension(
  const Vector& externalNormal, 
  const Vector& force, Scalar* values)
{
  // TODO
}

template <>
void ElasticSystemCommon<Space2>::FixedBoundaryInfoFunctor::SetVelocity(
  const Vector& externalVelocity, Scalar* values)
{
  values[3] = externalVelocity.x;
  values[4] = externalVelocity.y;
}

template <>
void ElasticSystemCommon<Space3>::FixedBoundaryInfoFunctor::SetVelocity(
  const Vector& externalVelocity, Scalar* values)
{
  values[6] = externalVelocity.x;
  values[7] = externalVelocity.y;
  values[8] = externalVelocity.z;
}

