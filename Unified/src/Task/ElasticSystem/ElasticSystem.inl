void ElasticSystem<Space2>::ValueType::SetTension(Scalar xx, Scalar yy, Scalar xy)
{
  values[0] = xx;
  values[1] = yy;
  values[2] = xy;
}

void ElasticSystem<Space3>::ValueType::SetTension(Scalar xx, Scalar yy, Scalar zz, Scalar xy, Scalar yz, Scalar xz)
{
  values[0] = xx;
  values[1] = yy;
  values[2] = zz;
  values[3] = xy;
  values[4] = yz;
  values[5] = xz;
}

Space3::Scalar ElasticSystem<Space3>::ValueType::GetPressure() const
{
  return -(GetXX() + GetYY() + GetZZ()) / Scalar(3.0);
}

Space3::Scalar ElasticSystem<Space3>::ValueType::GetDeviatorSquare() const
{
  // wiki http://en.wikipedia.org/wiki/Cauchy_stress_tensor
  return Scalar(1.0 / 3) * (Sqr(GetXX() - GetYY()) + Sqr(GetXX() - GetZZ()) + Sqr(GetYY() - GetZZ())) + 
          Scalar(2.0) * (Sqr(GetXY()) + Sqr(GetXZ()) + Sqr(GetYZ()));
}

Space2::Scalar ElasticSystem<Space2>::ValueType::GetPressure() const
{
  return -(GetXX() + GetYY()) / Scalar(2.0);
}

Space2::Scalar ElasticSystem<Space2>::ValueType::GetDeviatorSquare() const
{
  return Scalar(0.5) * (Sqr(GetXX()) + Sqr(GetYY())) - 
          GetXX() * GetYY() + 
          2.0 * Sqr(GetXY());
}

Space2::Vector ElasticSystem<Space2>::ValueType::GetForce(const Vector& normal) const
{
  return Space2::Vector(
    GetXX() * normal.x + GetXY() * normal.y, 
    GetXY() * normal.x + GetYY() * normal.y);
}

Space3::Vector ElasticSystem<Space3>::ValueType::GetForce(const Vector& normal) const
{
  return Space3::Vector(
    GetXX() * normal.x + GetXY() * normal.y + GetXZ() * normal.z, 
    GetXY() * normal.x + GetYY() * normal.y + GetYZ() * normal.z,
    GetXZ() * normal.x + GetYZ() * normal.y + GetZZ() * normal.z);
}

Space2::Tensor ElasticSystem<Space2>::ValueType::GetTension() const
{
  return Space2::Tensor(GetXX(), GetXY(), GetYY());
}

Space3::Scalar ElasticSystem<Space3>::ValueType::GetZZ() const
{
  return values[2];
}

Space3::Scalar ElasticSystem<Space3>::ValueType::GetYZ() const
{
  return values[4];
}

Space3::Scalar ElasticSystem<Space3>::ValueType::GetXZ() const
{
  return values[5];
}

Space3::Tensor ElasticSystem<Space3>::ValueType::GetTension() const
{
  return Space3::Tensor(values[0], values[3], values[5], values[1], values[4], values[2]);
}

void ElasticSystem<Space2>::BuildEdgeTransformMatrix(Vector edgeVertices[2], Scalar *transformMatrix)
{
  Vector normal  = Vector((edgeVertices[1] - edgeVertices[0]).y, -(edgeVertices[1] - edgeVertices[0]).x).GetNorm();
  Vector tangent = Vector(-normal.y, normal.x);

  std::fill_n(transformMatrix, dimsCount * dimsCount, Scalar(0.0));

  transformMatrix[0 ] = Sqr(normal.x);
  transformMatrix[1 ] = Sqr(normal.y);
  transformMatrix[5 ] = Sqr(normal.y);
  transformMatrix[6 ] = Sqr(normal.x);

  transformMatrix[2 ] = -Scalar(2.0) * normal.x * normal.y;
  transformMatrix[7 ] =  Scalar(2.0) * normal.x * normal.y;

  transformMatrix[10] =  normal  .x  * normal  .y;
  transformMatrix[11] = -normal  .x  * normal  .y;

  transformMatrix[12] = Sqr(normal.x) - Sqr(normal.y);

  transformMatrix[18] =  normal.x;
  transformMatrix[19] = -normal.y;
  transformMatrix[23] =  normal.y;
  transformMatrix[24] =  normal.x;
}

void ElasticSystem<Space2>::BuildEdgeTransformMatrixInv(Vector edgeVertices[2], Scalar *transformMatrixInv)
{
  Vector normal  = Vector((edgeVertices[1] - edgeVertices[0]).y, -(edgeVertices[1] - edgeVertices[0]).x).GetNorm();
  Vector tangent = Vector(-normal.y, normal.x);

  std::fill_n(transformMatrixInv, dimsCount * dimsCount, Scalar(0.0));

  transformMatrixInv[0 ] = Sqr(normal.x);
  transformMatrixInv[1 ] = Sqr(normal.y);
  transformMatrixInv[5 ] = Sqr(normal.y);
  transformMatrixInv[6 ] = Sqr(normal.x);

  transformMatrixInv[2 ] =  Scalar(2.0) * normal.x * normal.y;
  transformMatrixInv[7 ] = -Scalar(2.0) * normal.x * normal.y;

  transformMatrixInv[10] = -normal  .x  * normal  .y;
  transformMatrixInv[11] =  normal  .x  * normal  .y;

  transformMatrixInv[12] = Sqr(normal.x) - Sqr(normal.y);

  transformMatrixInv[18] =  normal.x;
  transformMatrixInv[19] =  normal.y;
  transformMatrixInv[23] = -normal.y;
  transformMatrixInv[24] =  normal.x;
}

void ElasticSystem<Space3>::BuildFaceTransformMatrix(Vector faceVertices[3], Scalar* transformMatrix)
{
  Vector normal   = ((faceVertices[1] - faceVertices[0]) ^ (faceVertices[2] - faceVertices[0])).GetNorm();
  Vector tangent0 = (faceVertices[1] - faceVertices[0]).GetNorm();
  Vector tangent1 = (normal ^ tangent0).GetNorm();

  std::fill_n(transformMatrix, dimsCount * dimsCount, 0);

  transformMatrix[0 ] = Sqr(normal  .x);
  transformMatrix[1 ] = Sqr(tangent0.x);
  transformMatrix[2 ] = Sqr(tangent1.x);
  transformMatrix[9 ] = Sqr(normal  .y);
  transformMatrix[10] = Sqr(tangent0.y);
  transformMatrix[11] = Sqr(tangent1.y);
  transformMatrix[18] = Sqr(normal  .z);
  transformMatrix[19] = Sqr(tangent0.z);
  transformMatrix[20] = Sqr(tangent1.z);

  transformMatrix[3 ] = Scalar(2.0) * normal  .x  * tangent0.x;
  transformMatrix[4 ] = Scalar(2.0) * tangent0.x  * tangent1.x;
  transformMatrix[5 ] = Scalar(2.0) * normal  .x  * tangent1.x;
  transformMatrix[12] = Scalar(2.0) * normal  .y  * tangent0.y;
  transformMatrix[13] = Scalar(2.0) * tangent0.y  * tangent1.y;
  transformMatrix[14] = Scalar(2.0) * normal  .y  * tangent1.y;
  transformMatrix[21] = Scalar(2.0) * normal  .z  * tangent0.z;
  transformMatrix[22] = Scalar(2.0) * tangent0.z  * tangent1.z;
  transformMatrix[23] = Scalar(2.0) * normal  .z  * tangent1.z;

  transformMatrix[27] = normal  .y  * normal  .x;
  transformMatrix[28] = tangent0.y  * tangent0.x;
  transformMatrix[29] = tangent1.y  * tangent1.x;
  transformMatrix[36] = normal  .z  * normal  .y;
  transformMatrix[37] = tangent0.z  * tangent0.y;
  transformMatrix[38] = tangent1.z  * tangent1.y;
  transformMatrix[45] = normal  .z  * normal  .x;
  transformMatrix[46] = tangent0.z  * tangent0.x;
  transformMatrix[47] = tangent1.z  * tangent1.x;

  transformMatrix[30] = normal  .y  * tangent0.x + normal  .x  * tangent0.y;
  transformMatrix[31] = tangent0.y  * tangent1.x + tangent0.x  * tangent1.y;
  transformMatrix[32] = normal  .y  * tangent1.x + normal  .x  * tangent1.y;
  transformMatrix[39] = normal  .z  * tangent0.y + normal  .y  * tangent0.z;
  transformMatrix[40] = tangent0.z  * tangent1.y + tangent0.y  * tangent1.z;
  transformMatrix[41] = normal  .z  * tangent1.y + normal  .y  * tangent1.z;
  transformMatrix[48] = normal  .z  * tangent0.x + normal  .x  * tangent0.z;
  transformMatrix[49] = tangent0.z  * tangent1.x + tangent0.x  * tangent1.z;
  transformMatrix[50] = normal  .z  * tangent1.x + normal  .x  * tangent1.z;

  transformMatrix[60] = normal  .x;
  transformMatrix[61] = tangent0.x;
  transformMatrix[62] = tangent1.x;
  transformMatrix[69] = normal  .y;
  transformMatrix[70] = tangent0.y;
  transformMatrix[71] = tangent1.y;
  transformMatrix[78] = normal  .z;
  transformMatrix[79] = tangent0.z;
  transformMatrix[80] = tangent1.z;
}

void ElasticSystem<Space3>::BuildFaceTransformMatrixInv(Vector faceVertices[3], Scalar *transformMatrixInv)
{
  Vector normal = ((faceVertices[1] - faceVertices[0]) ^ (faceVertices[2] - faceVertices[0])).GetNorm();
  Vector tangent0 = (faceVertices[1] - faceVertices[0]).GetNorm();
  Vector tangent1 = (normal ^ tangent0).GetNorm();
  std::swap(tangent0.x, normal.y);
  std::swap(tangent1.x, normal.z);
  std::swap(tangent1.y, tangent0.z);

  std::fill_n(transformMatrixInv, dimsCount * dimsCount, 0);

  transformMatrixInv[0 ] = Sqr(normal  .x);
  transformMatrixInv[1 ] = Sqr(tangent0.x);
  transformMatrixInv[2 ] = Sqr(tangent1.x);
  transformMatrixInv[9 ] = Sqr(normal  .y);
  transformMatrixInv[10] = Sqr(tangent0.y);
  transformMatrixInv[11] = Sqr(tangent1.y);
  transformMatrixInv[18] = Sqr(normal  .z);
  transformMatrixInv[19] = Sqr(tangent0.z);
  transformMatrixInv[20] = Sqr(tangent1.z);

  transformMatrixInv[3 ] = Scalar(2.0) * normal  .x  * tangent0.x;
  transformMatrixInv[4 ] = Scalar(2.0) * tangent0.x  * tangent1.x;
  transformMatrixInv[5 ] = Scalar(2.0) * normal  .x  * tangent1.x;
  transformMatrixInv[12] = Scalar(2.0) * normal  .y  * tangent0.y;
  transformMatrixInv[13] = Scalar(2.0) * tangent0.y  * tangent1.y;
  transformMatrixInv[14] = Scalar(2.0) * normal  .y  * tangent1.y;
  transformMatrixInv[21] = Scalar(2.0) * normal  .z  * tangent0.z;
  transformMatrixInv[22] = Scalar(2.0) * tangent0.z  * tangent1.z;
  transformMatrixInv[23] = Scalar(2.0) * normal  .z  * tangent1.z;

  transformMatrixInv[27] = normal  .y  * normal  .x;
  transformMatrixInv[28] = tangent0.y  * tangent0.x;
  transformMatrixInv[29] = tangent1.y  * tangent1.x;
  transformMatrixInv[36] = normal  .z  * normal  .y;
  transformMatrixInv[37] = tangent0.z  * tangent0.y;
  transformMatrixInv[38] = tangent1.z  * tangent1.y;
  transformMatrixInv[45] = normal  .z  * normal  .x;
  transformMatrixInv[46] = tangent0.z  * tangent0.x;
  transformMatrixInv[47] = tangent1.z  * tangent1.x;

  transformMatrixInv[30] = normal  .y  * tangent0.x + normal  .x  * tangent0.y;
  transformMatrixInv[31] = tangent0.y  * tangent1.x + tangent0.x  * tangent1.y;
  transformMatrixInv[32] = normal  .y  * tangent1.x + normal  .x  * tangent1.y;
  transformMatrixInv[39] = normal  .z  * tangent0.y + normal  .y  * tangent0.z;
  transformMatrixInv[40] = tangent0.z  * tangent1.y + tangent0.y  * tangent1.z;
  transformMatrixInv[41] = normal  .z  * tangent1.y + normal  .y  * tangent1.z;
  transformMatrixInv[48] = normal  .z  * tangent0.x + normal  .x  * tangent0.z;
  transformMatrixInv[49] = tangent0.z  * tangent1.x + tangent0.x  * tangent1.z;
  transformMatrixInv[50] = normal  .z  * tangent1.x + normal  .x  * tangent1.z;

  transformMatrixInv[60] = normal  .x;
  transformMatrixInv[61] = tangent0.x;
  transformMatrixInv[62] = tangent1.x;
  transformMatrixInv[69] = normal  .y;
  transformMatrixInv[70] = tangent0.y;
  transformMatrixInv[71] = tangent1.y;
  transformMatrixInv[78] = normal  .z;
  transformMatrixInv[79] = tangent0.z;
  transformMatrixInv[80] = tangent1.z;
}

void ElasticSystem<Space3>::BuildZMatrix(const MediumParameters& mediumParameters, Scalar* zMatrix)
{
  std::fill_n(zMatrix, dimsCount * dimsCount, Scalar(0.0));

  zMatrix[ 8] = -mediumParameters.lambda;
  zMatrix[17] = -mediumParameters.lambda;
  zMatrix[26] = -(mediumParameters.lambda + Scalar(2.0) * mediumParameters.mju);
  zMatrix[43] = -mediumParameters.mju;
  zMatrix[51] = -mediumParameters.mju;
  zMatrix[59] = -mediumParameters.invRho;
  zMatrix[67] = -mediumParameters.invRho;
  zMatrix[74] = -mediumParameters.invRho;

  for (int row = 0; row < dimsCount; ++row)
  {
    // diagonal elements
    zMatrix[row + row * dimsCount] += mediumParameters.flowVelocity.z;
  }
}
