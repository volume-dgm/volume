#pragma once

template <typename VolumeMesh>
class SourceFunctionGetter
{
public:
  typedef typename VolumeMesh::Space Space;
  typedef typename VolumeMesh::SystemT System;
  SPACE_TYPEDEFS

    SourceFunctionGetter(typename System::SourceFunctorT* functor,
    const Scalar time,
    const VolumeMesh* const volumeMesh,
    IndexType cellIndex):
    functor(functor),
    time(time),
    volumeMesh(volumeMesh)
  {
    volumeMesh->GetCellVertices(cellIndex, cellVertices);
  }

  void operator()(const Vector& refPoint, Scalar* values)
  {
    Vector globalPoint = volumeMesh->RefToGlobalVolumeCoords(refPoint, cellVertices);
    if (functor)
    {
      functor->operator()(globalPoint, time, values);
    } else
    {
      for (IndexType valueIndex = 0; valueIndex < VolumeMesh::dimsCount; ++valueIndex)
      {
        values[valueIndex] = Scalar(0.0);
      }
    }
  }

private:
  typename System::SourceFunctorT* functor;
  Scalar time;
  const VolumeMesh* const volumeMesh;
  Vector cellVertices[Space::NodesPerCell];
};
