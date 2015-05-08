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

void ElasticSystem<Space2>::BuildEdgeTransformMatrix(Vector edgeVertices[Space::NodesPerEdge], 
  Eigen::Matrix<Scalar, dimsCount, dimsCount>& transformMatrix)
{
  Vector normal  = Vector((edgeVertices[1] - edgeVertices[0]).y, -(edgeVertices[1] - edgeVertices[0]).x).GetNorm();

  transformMatrix <<
          Sqr(normal.x),         Sqr(normal.y), -Scalar(2.0) * normal.x * normal.y,         0,         0,
          Sqr(normal.y),         Sqr(normal.x),  Scalar(2.0) * normal.x * normal.y,         0,         0,
    normal.x * normal.y, -normal.x  * normal.y,      Sqr(normal.x) - Sqr(normal.y),         0,         0,
                      0,                     0,                                   0, normal.x, -normal.y,
                      0,                     0,                                   0, normal.y,  normal.x;
}

void ElasticSystem<Space2>::BuildEdgeTransformMatrixInv(Vector edgeVertices[Space::NodesPerEdge], 
  Eigen::Matrix<Scalar, dimsCount, dimsCount>& transformMatrixInv)
{
  Vector normal  = Vector((edgeVertices[1] - edgeVertices[0]).y, -(edgeVertices[1] - edgeVertices[0]).x).GetNorm();

  transformMatrixInv <<
            Sqr(normal.x),        Sqr(normal.y),  Scalar(2.0) * normal.x * normal.y,         0,        0,
            Sqr(normal.y),        Sqr(normal.x), -Scalar(2.0) * normal.x * normal.y,         0,        0,
    -normal.x  * normal.y, normal.x  * normal.y,      Sqr(normal.x) - Sqr(normal.y),         0,        0,
                        0,                    0,                                  0,  normal.x, normal.y,
                        0,                    0,                                  0, -normal.y, normal.x;
}

void ElasticSystem<Space3>::BuildFaceTransformMatrix(Vector faceVertices[Space::NodesPerFace], 
  Eigen::Matrix<Scalar, dimsCount, dimsCount>& transformMatrix)
{
  Vector normal   = ((faceVertices[1] - faceVertices[0]) ^ (faceVertices[2] - faceVertices[0])).GetNorm();
  Vector tangent0 = (faceVertices[1] - faceVertices[0]).GetNorm();
  Vector tangent1 = (normal ^ tangent0).GetNorm();

  transformMatrix <<
          Sqr(normal.x),         Sqr(tangent0.x),         Sqr(tangent1.x),           Scalar(2.0) * normal.x * tangent0.x,             Scalar(2.0) * tangent0.x * tangent1.x,           Scalar(2.0) * normal.x * tangent1.x,        0,          0,          0,
          Sqr(normal.y),         Sqr(tangent0.y),         Sqr(tangent1.y),           Scalar(2.0) * normal.y * tangent0.y,             Scalar(2.0) * tangent0.y * tangent1.y,           Scalar(2.0) * normal.y * tangent1.y,        0,          0,          0,
          Sqr(normal.z),         Sqr(tangent0.z),         Sqr(tangent1.z),           Scalar(2.0) * normal.z * tangent0.z,             Scalar(2.0) * tangent0.z * tangent1.z,           Scalar(2.0) * normal.z * tangent1.z,        0,          0,          0,
    normal.y * normal.x, tangent0.y * tangent0.x, tangent1.y * tangent1.x, normal.y * tangent0.x + normal.x * tangent0.y, tangent0.y * tangent1.x + tangent0.x * tangent1.y, normal.y * tangent1.x + normal.x * tangent1.y,        0,          0,          0,
    normal.z * normal.y, tangent0.z * tangent0.y, tangent1.z * tangent1.y, normal.z * tangent0.y + normal.y * tangent0.z, tangent0.z * tangent1.y + tangent0.y * tangent1.z, normal.z * tangent1.y + normal.y * tangent1.z,        0,          0,          0,
    normal.z * normal.x, tangent0.z * tangent0.x, tangent1.z * tangent1.x, normal.z * tangent0.x + normal.x * tangent0.z, tangent0.z * tangent1.x + tangent0.x * tangent1.z, normal.z * tangent1.x + normal.x * tangent1.z,        0,          0,          0,
                      0,                       0,                       0,                                             0,                                                 0,                                             0, normal.x, tangent0.x, tangent1.x,
                      0,                       0,                       0,                                             0,                                                 0,                                             0, normal.y, tangent0.y, tangent1.y,
                      0,                       0,                       0,                                             0,                                                 0,                                             0, normal.z, tangent0.z, tangent1.z;
}

void ElasticSystem<Space3>::BuildFaceTransformMatrixInv(Vector faceVertices[Space::NodesPerFace], 
  Eigen::Matrix<Scalar, dimsCount, dimsCount>& transformMatrixInv)
{
  Vector normal = ((faceVertices[1] - faceVertices[0]) ^ (faceVertices[2] - faceVertices[0])).GetNorm();
  Vector tangent0 = (faceVertices[1] - faceVertices[0]).GetNorm();
  Vector tangent1 = (normal ^ tangent0).GetNorm();
  std::swap(tangent0.x, normal.y);
  std::swap(tangent1.x, normal.z);
  std::swap(tangent1.y, tangent0.z);

  transformMatrixInv << 
          Sqr(normal.x),         Sqr(tangent0.x),         Sqr(tangent1.x),          Scalar(2.0) * normal.x  * tangent0.x,            Scalar(2.0) * tangent0.x  * tangent1.x,          Scalar(2.0) * normal.x  * tangent1.x,        0,          0,          0,
          Sqr(normal.y),         Sqr(tangent0.y),         Sqr(tangent1.y),          Scalar(2.0) * normal.y  * tangent0.y,            Scalar(2.0) * tangent0.y  * tangent1.y,          Scalar(2.0) * normal.y  * tangent1.y,        0,          0,          0,
          Sqr(normal.z),         Sqr(tangent0.z),         Sqr(tangent1.z),          Scalar(2.0) * normal.z  * tangent0.z,            Scalar(2.0) * tangent0.z  * tangent1.z,          Scalar(2.0) * normal.z  * tangent1.z,        0,          0,          0,
    normal.y * normal.x, tangent0.y * tangent0.x, tangent1.y * tangent1.x, normal.y * tangent0.x + normal.x * tangent0.y, tangent0.y * tangent1.x + tangent0.x * tangent1.y, normal.y * tangent1.x + normal.x * tangent1.y,        0,          0,          0,
    normal.z * normal.y, tangent0.z * tangent0.y, tangent1.z * tangent1.y, normal.z * tangent0.y + normal.y * tangent0.z, tangent0.z * tangent1.y + tangent0.y * tangent1.z, normal.z * tangent1.y + normal.y * tangent1.z,        0,          0,          0,
    normal.z * normal.x, tangent0.z * tangent0.x, tangent1.z * tangent1.x, normal.z * tangent0.x + normal.x * tangent0.z, tangent0.z * tangent1.x + tangent0.x * tangent1.z, normal.z * tangent1.x + normal.x * tangent1.z,        0,          0,          0,
                      0,                       0,                       0,                                             0,                                                 0,                                             0, normal.x, tangent0.x, tangent1.x,
                      0,                       0,                       0,                                             0,                                                 0,                                             0, normal.y, tangent0.y, tangent1.y, 
                      0,                       0,                       0,                                             0,                                                 0,                                             0, normal.z, tangent0.z, tangent1.z;
}

void ElasticSystem<Space3>::BuildZMatrix(const MediumParameters& mediumParameters, 
  Eigen::Matrix<Scalar, dimsCount, dimsCount>& zMatrix)
{
  zMatrix <<
    mediumParameters.flowVelocity.z,                               0,                               0,                               0,                               0,                               0,                               0,                               0,  -mediumParameters.lambda                                      ,
                                  0, mediumParameters.flowVelocity.z,                               0,                               0,                               0,                               0,                               0,                               0,  -mediumParameters.lambda                                      ,
                                  0,                               0, mediumParameters.flowVelocity.z,                               0,                               0,                               0,                               0,                               0, -(mediumParameters.lambda + Scalar(2.0) * mediumParameters.mju),
                                  0,                               0,                               0, mediumParameters.flowVelocity.z,                               0,                               0,                               0,                               0,                                                               0,
                                  0,                               0,                               0,                               0, mediumParameters.flowVelocity.z,                               0,                               0,           -mediumParameters.mju,                                                               0,
                                  0,                               0,                               0,                               0,                               0, mediumParameters.flowVelocity.z,           -mediumParameters.mju,                               0,                                                               0,
                                  0,                               0,                               0,                               0,                               0,        -mediumParameters.invRho, mediumParameters.flowVelocity.z,                               0,                                                               0,
                                  0,                               0,                               0,                               0,        -mediumParameters.invRho,                               0,                               0, mediumParameters.flowVelocity.z,                                                               0,
                                  0,                               0,        -mediumParameters.invRho,                               0,                               0,                               0,                               0,                               0,                                 mediumParameters.flowVelocity.z;
}
