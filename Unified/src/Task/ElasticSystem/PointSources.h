#pragma once
#include <limits>

template <typename Space>
class PointSource
{
public:
  SPACE_TYPEDEFS

  PointSource(Scalar tensionDimensionlessMult, Scalar velocityDimensionlessMult):
    tensionDimensionlessMult(tensionDimensionlessMult),
    velocityDimensionlessMult(velocityDimensionlessMult)
  {}

  virtual ~PointSource() {}
  virtual Vector GetPoint() const = 0;
  virtual void operator()(Scalar time, Scalar* values) const = 0;

protected:
  Scalar tensionDimensionlessMult;
  Scalar velocityDimensionlessMult;
};

template <typename ElasticSystem>
class RiekerPointSource: public PointSource<typename ElasticSystem::Space>
{
public:
  typedef typename ElasticSystem::Space Space;
  SPACE_TYPEDEFS
  const static IndexType dimsCount = ElasticSystem::dimsCount;

  RiekerPointSource(const Vector& point,
    Scalar peakFrequency, Vector acceleration, Scalar latency,
    Scalar tensionDimensionlessMult, Scalar velocityDimensionlessMult): 
    PointSource<Space>(tensionDimensionlessMult, velocityDimensionlessMult), 
    point(point), peakFrequency(peakFrequency), acceleration(acceleration), latency(latency)
  {
    latency = std::max(latency, sqrt(1.5) / (pi * peakFrequency) * 3);
  }

  Vector GetPoint() const
  {
    return point;
  }

  void operator()(Scalar time, Scalar* values) const
  {
    Scalar arg = Sqr(pi * peakFrequency * (time - latency));
    Scalar mult = (1 - 2 * arg) * exp(-arg);
    std::fill_n(values, dimsCount, Scalar(0.0));
    
    SetVelocity(Overload<Space>(), mult, values);
  }

private:
  void SetVelocity(Overload<Space2>, typename Space2::Scalar mult, typename Space2::Scalar* values) const
  {
    values[3] = mult * acceleration.x;
    values[4] = mult * acceleration.y;
  }

  void SetVelocity(Overload<Space3>, typename Space3::Scalar mult, typename Space3::Scalar* values) const
  {
    values[6] = mult * acceleration.x;
    values[7] = mult * acceleration.y;
    values[8] = mult * acceleration.z;
  }

  Vector point; 
  Scalar peakFrequency;
  Vector acceleration;
  Scalar latency; 
};
