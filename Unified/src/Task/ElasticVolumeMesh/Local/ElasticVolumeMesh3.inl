template<typename FunctionSpace>
ElasticVolumeMesh<Space3, FunctionSpace>::ElasticVolumeMesh(DifferentialSolver<Scalar>* solver, Scalar tolerance, int hierarchyLevelsCount):
  ElasticVolumeMeshCommon<Space3, FunctionSpace>(solver, tolerance, hierarchyLevelsCount)
{
}

template<typename FunctionSpace>
void ElasticVolumeMesh<Space3, FunctionSpace>::UnfoldMesh(Scalar minHeight, IndexType iterationsCount)
{
}

template<typename FunctionSpace>
void ElasticVolumeMesh<Space3, FunctionSpace>::
  FindDestructions(FacePairIndices* contactFaces, IndexType* contactFacesCount, IndexType contactTypesCount,
  std::vector<bool>* isContactBroken, std::vector<bool>* isCellBroken)
{
  if (allowContinuousDestruction)
  {
    HandleContinuousDestruction();
    for (IndexType cellIndex = 0; cellIndex < volumeMesh.cellMediumParameters.size(); ++cellIndex)
    {
      (*isCellBroken)[cellIndex] = volumeMesh.cellMediumParameters[cellIndex].destroyed;
    }
  }

  if (allowDiscreteDestruction)
  {
    IndexType offset = 0;
    for (IndexType contactTypeIndex = 0; contactTypeIndex < contactTypesCount; ++contactTypeIndex)
    {
      Scalar maxLongitudinalStress;
      Scalar maxShearStress;

      IndexType dynamicBoundaryType = volumeMesh.system.GetContactDynamicBoundaryType(contactTypeIndex);
      if (dynamicBoundaryType != IndexType(-1))
      {
        volumeMesh.system.GetContactCriticalInfo(contactTypeIndex, maxShearStress, maxLongitudinalStress);
      }
      else
      {
        continue;
      }

      struct DetectionApproaches
      {
        enum Types
        {
          Volumetric = 0x01, Surface = 0x02, Combined = 0x01 | 0x02
        };
      };

      const IndexType detectionApproach =
        DetectionApproaches::Volumetric;
      // DetectionApproaches::Surface;
      // DetectionApproaches::Combined;

      for (IndexType contactFaceIndex = 0; contactFaceIndex < contactFacesCount[contactTypeIndex]; ++contactFaceIndex)
      {
        bool isCurrentContactBroken = false;
        FaceLocationPair contactFacesLocation = volumeMesh.GetFaceLocation(contactFaces[offset + contactFaceIndex]);

        for (IndexType faceIndex = 0; faceIndex < 2; ++faceIndex)
        {
          IndexType cellIndex = contactFacesLocation.faces[faceIndex].cellIndex;
          IndexType faceNumber = contactFacesLocation.faces[faceIndex].faceNumber;

          // average value in the cell was used, because of presence of spurious oscillations
          Elastic elastic = GetAverageCellElastic(cellIndex);

          Scalar principalStresses[Space3::Dimension];
          Vector principalNormals[Space3::Dimension];
          elastic.GetTension().GetEigenValues(principalStresses, principalNormals);

          IndexType maxStressIndex = IndexType(-1);

          for (IndexType stressIndex = 0; stressIndex < Space3::Dimension; ++stressIndex)
          {
            /* only extension */
            if (principalStresses[stressIndex] > fabs(maxLongitudinalStress))
            {
              if (maxStressIndex == IndexType(-1) || principalStresses[stressIndex] > principalStresses[maxStressIndex])
              {
                maxStressIndex = stressIndex;
              }
            }
          }

          if (maxStressIndex == IndexType(-1)) continue;
          const Vector& mainNormal = principalNormals[maxStressIndex];

          // minimun cos of angle between normal to crack and edge tangent
          Scalar maxSinOfAngle = 0;
          IndexType bestFaceNumber = IndexType(-1);

          for (IndexType face = 0; face < Space3::FacesPerCell; face++)
          {
            Vector faceNormal = volumeMesh.GetFaceExternalNormal(cellIndex, face).GetNorm();

            if (maxSinOfAngle < fabs(mainNormal * faceNormal))
            {
              maxSinOfAngle = fabs(mainNormal * faceNormal);
              bestFaceNumber = face;
            }
          }

          if (bestFaceNumber == faceNumber)
          {
            isCurrentContactBroken = true;
            break;
          }
        }

        if (isCurrentContactBroken && !(*isContactBroken)[offset + contactFaceIndex])
        {
          (*isContactBroken)[offset + contactFaceIndex] = true;
          for (IndexType faceIndex = 0; faceIndex < 2; ++faceIndex)
          {
            IndexType cellIndex = contactFacesLocation.faces[faceIndex].cellIndex;
            IndexType faceNumber = contactFacesLocation.faces[faceIndex].faceNumber;
            if (volumeMesh.cellMediumParameters[cellIndex].destroyed) continue;

            // if this cell is deleted due to erosion
            if (!volumeMesh.isCellAvailable[cellIndex]) continue;

            IndexType correspondingCellIndex = volumeMesh.GetCorrespondingCellIndex(cellIndex, faceNumber);
            IndexType correspondingFaceNumber = volumeMesh.GetCorrespondingFaceNumber(cellIndex, faceNumber);
            if (correspondingCellIndex == IndexType(-1) ||
              correspondingFaceNumber == IndexType(-1) ||
              volumeMesh.cellMediumParameters[correspondingCellIndex].destroyed) continue;

            // check if it is last non-destructed edge for this cell
            IndexType incidentCellsCount = 0;
            for (IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; ++faceNumber)
            {
              if (volumeMesh.GetCorrespondingCellIndex(cellIndex, faceNumber) != IndexType(-1) &&
                volumeMesh.GetCorrespondingFaceNumber(cellIndex, faceNumber) != IndexType(-1))
              {
                incidentCellsCount++;
              }
            }

            if (incidentCellsCount == 1)
            {
              DestroyCellMaterial(cellIndex, volumeMesh.cellMediumParameters[cellIndex].plasticity.powderShearMult);
              (*isCellBroken)[cellIndex] = true;
              continue;
            }

            if (volumeMesh.isCellAvailable[cellIndex])
            {
              //  volumeMesh.AddToAABBTree(cellIndex);
            }

            if (volumeMesh.isCellAvailable[correspondingCellIndex])
            {
              //  volumeMesh.AddToAABBTree(correspondingCellIndex);
            }

            DestroyFace(cellIndex, faceNumber, dynamicBoundaryType);
          }
        }
      }
      offset += contactFacesCount[contactTypeIndex];
    }
  }

  HandleMaterialErosion();
}

template<typename MeshType>
struct ContinuousDestructionCorrector3
{
  typedef typename MeshType::Scalar    Scalar;
  typedef typename MeshType::Vector    Vector;
  typedef typename MeshType::Tensor    Tensor;
  typedef typename MeshType::IndexType IndexType;
  typedef typename MeshType::Elastic   Elastic;

  ContinuousDestructionCorrector3(MeshType* mesh, IndexType cellIndex)
  {
    this->mesh = mesh;
    this->cellIndex = cellIndex;
  }

  void operator()(const Vector& refPoint, Scalar* values) const
  {
    bool brittle = mesh->volumeMesh.cellMediumParameters[cellIndex].plasticity.brittle;
    Elastic elastic = mesh->InterpolateElasticRef(cellIndex, refPoint);

    if (brittle)
    {
      const Scalar k0 = mesh->volumeMesh.cellMediumParameters[cellIndex].plasticity.k0;
      const Scalar  a = mesh->volumeMesh.cellMediumParameters[cellIndex].plasticity.a;
      const Scalar pressure = elastic.GetPressure();
      const Scalar k = k0 + a * pressure;

      if (mesh->ProcessPlasticity(k, elastic, false))
      {
        mesh->DestroyCellMaterial(cellIndex, mesh->volumeMesh.cellMediumParameters[cellIndex].plasticity.powderShearMult);
      }
    }

    if (mesh->volumeMesh.cellMediumParameters[cellIndex].destroyed)
    {
      Scalar principalStresses[Space3::Dimension];
      Vector principalNormals[Space3::Dimension];
      elastic.GetTension().GetEigenValues(principalStresses, principalNormals);

      for (IndexType stressIndex = 0; stressIndex < Space3::Dimension; ++stressIndex)
      {
        principalStresses[stressIndex] = std::min<Scalar>(0, principalStresses[stressIndex]);
      }

      Tensor t;
      t.xx = principalStresses[0];
      t.yy = principalStresses[1];
      t.zz = principalStresses[2];

      std::swap(principalNormals[0].y, principalNormals[1].x);
      std::swap(principalNormals[0].z, principalNormals[2].x);
      std::swap(principalNormals[1].z, principalNormals[2].y);
      
      t = t.RotateAxes(principalNormals[0], principalNormals[1], principalNormals[2]);
      // rotate from local principal tensions axes to global coordinates
      elastic.SetTension(t);
    }

    std::copy(elastic.values, elastic.values + mesh->volumeMesh.dimsCount, values);
  }
  MeshType* mesh;
  IndexType cellIndex;
};

template<typename FunctionSpace>
void ElasticVolumeMesh<Space3, FunctionSpace>::HandleContinuousDestruction()
{
  typedef ContinuousDestructionCorrector3< ElasticVolumeMesh<Space3, FunctionSpace> > DestructionCorrectorType;
  ApplyCorrector< ElasticVolumeMesh<Space3, FunctionSpace>, DestructionCorrectorType>(this);
}

template<typename FunctionSpace>
void ElasticVolumeMesh<Space3, FunctionSpace>::MakeSnapshot(Elastic* destData,
  const Vector& origin, const Vector& spacing,
  const Vector& boxPoint1, const Vector& boxPoint2,
  bool halfStepSolution)
{
  IntAABB globalSnapshotArea;
  AABB snapshotArea;
  snapshotArea.Set(boxPoint1, boxPoint2);

  IndexVector snapshotSize = ::ComputeSnapshotAABBSize<Space3>(origin, spacing, snapshotArea,
    &globalSnapshotArea);

  int globalIxmin = globalSnapshotArea.boxPoint1.x;
  int globalIxmax = globalSnapshotArea.boxPoint2.x;
  int globalIymin = globalSnapshotArea.boxPoint1.y;
  int globalIymax = globalSnapshotArea.boxPoint2.y;
  int globalIzmin = globalSnapshotArea.boxPoint1.z;
  int globalIzmax = globalSnapshotArea.boxPoint2.z;

  for (IndexType elasticIndex = 0; elasticIndex < snapshotSize.GetVolume(); ++elasticIndex)
  {
    std::fill_n(destData[elasticIndex].values, volumeMesh.dimsCount, Scalar(0.0));
  }

  for (IndexType cellIndex = 0; cellIndex < volumeMesh.cells.size(); ++cellIndex)
  {
    Vector points[Space::NodesPerCell];
    volumeMesh.GetCellVertices(cellIndex, points);

    Scalar
      xmin = points[0].x, xmax = points[0].x,
      ymin = points[0].y, ymax = points[0].y,
      zmin = points[0].z, zmax = points[0].z;

    for (IndexType nodeNumber = 0; nodeNumber < Space::NodesPerCell; ++nodeNumber)
    {
      if (points[nodeNumber].x < xmin) xmin = points[nodeNumber].x;
      if (points[nodeNumber].x > xmax) xmax = points[nodeNumber].x;
      if (points[nodeNumber].y < ymin) ymin = points[nodeNumber].y;
      if (points[nodeNumber].y > ymax) ymax = points[nodeNumber].y;
      if (points[nodeNumber].z < zmin) zmin = points[nodeNumber].z;
      if (points[nodeNumber].z > zmax) zmax = points[nodeNumber].z;
    }

    int ixmin = int((xmin - origin.x) / spacing.x) - 1;
    int ixmax = int((xmax - origin.x) / spacing.x) + 1;

    int iymin = int((ymin - origin.y) / spacing.y) - 1;
    int iymax = int((ymax - origin.y) / spacing.y) + 1;

    int izmin = int((zmin - origin.z) / spacing.z) - 1;
    int izmax = int((zmax - origin.z) / spacing.z) + 1;

    for (int x = ixmin; x <= ixmax; x++)
    {
      for (int y = iymin; y <= iymax; y++)
      {
        for (int z = izmin; z <= izmax; z++)
        {
          if (x - globalIxmin >= 0 && x - globalIxmin < int(snapshotSize.x) &&
              y - globalIymin >= 0 && y - globalIymin < int(snapshotSize.y) &&              
              z - globalIzmin >= 0 && z - globalIzmin < int(snapshotSize.z))
          {
            Vector point(origin.x + spacing.x * Scalar(x),
                         origin.y + spacing.y * Scalar(y),
                         origin.z + spacing.z * Scalar(z));

            if (PointInCell<Scalar>(points, point) && snapshotArea.Includes(point))
            {
              Elastic e = InterpolateElastic(cellIndex, point, halfStepSolution);
              e.MakeDimension(tensionDimensionlessMult, velocityDimensionlessMult);
              destData[(z - globalIzmin) * snapshotSize.x * snapshotSize.y +
                       (y - globalIymin) * snapshotSize.x + 
                       (x - globalIxmin)] = e;
            }
          }
        }
      }
    }
  }
}

template<typename FunctionSpace>
void ElasticVolumeMesh<Space3, FunctionSpace>::FindDestructions(std::vector<bool>* isCellBroken)
{
  if (allowContinuousDestruction)
  {
    HandleContinuousDestruction();
    for (IndexType cellIndex = 0; cellIndex < volumeMesh.cellMediumParameters.size(); ++cellIndex)
    {
      (*isCellBroken)[cellIndex] = volumeMesh.cellMediumParameters[cellIndex].destroyed;
    }
  }

  if (allowDiscreteDestruction)
  {
    for (IndexType cellIndex = 0; cellIndex < volumeMesh.cells.size(); ++cellIndex)
    {
      Scalar maxSinOfAngle = 0;
      IndexType bestFaceNumber = IndexType(-1);
      IndexType dynamicBoundaryTypes[Space::FacesPerCell];

      for (IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; ++faceNumber)
      {
        Scalar maxLongitudinalStress;
        Scalar maxShearStress;

        IndexType correspondingCellIndex  = volumeMesh.GetCorrespondingCellIndex(cellIndex, faceNumber);
        IndexType correspondingFaceNumber = volumeMesh.GetCorrespondingFaceNumber(cellIndex, faceNumber);
        IndexType interactionType = volumeMesh.GetInteractionType(cellIndex, faceNumber);

        if (interactionType == IndexType(-1)) continue;

        dynamicBoundaryTypes[faceNumber] = volumeMesh.system.GetContactDynamicBoundaryType(interactionType);

        if (dynamicBoundaryTypes[faceNumber] != IndexType(-1))
        {
          volumeMesh.system.GetContactCriticalInfo(interactionType, maxShearStress, maxLongitudinalStress);
        } else
        {
          continue;
        }

        Elastic elastic = volumeMesh.GetFaceAverageSolution(cellIndex, faceNumber);
        //Elastic elastic = volumeMesh.GetCellAverageSolution(cellIndex);

        Scalar principalStresses[Space3::Dimension];
        Vector principalNormals[Space3::Dimension];
        elastic.GetTension().GetEigenValues(principalStresses, principalNormals);

        IndexType maxStressIndex = IndexType(-1);

        for (IndexType stressIndex = 0; stressIndex < Space3::Dimension; ++stressIndex)
        {
          if (maxStressIndex == IndexType(-1) || principalStresses[stressIndex] > principalStresses[maxStressIndex])
          {
            maxStressIndex = stressIndex;
          }
        }

        /* only extension */
        if (principalStresses[maxStressIndex] > fabs(maxLongitudinalStress) || correspondingCellIndex == IndexType(-1))
        {
          const Vector& mainNormal = principalNormals[maxStressIndex];
          Vector faceNormal = volumeMesh.GetFaceExternalNormal(cellIndex, faceNumber).GetNorm();

          if (maxSinOfAngle < fabs(mainNormal * faceNormal))
          {
            maxSinOfAngle = fabs(mainNormal * faceNormal);
            bestFaceNumber = faceNumber;
          }
        }
      }

      if (bestFaceNumber != IndexType(-1))
      {
        if (volumeMesh.cellMediumParameters[cellIndex].destroyed) continue;

        // if this cell is deleted due to erosion
        if (!volumeMesh.isCellAvailable[cellIndex]) continue;

        IndexType correspondingCellIndex = volumeMesh.GetCorrespondingCellIndex(cellIndex, bestFaceNumber);
        IndexType correspondingFaceNumber = volumeMesh.GetCorrespondingFaceNumber(cellIndex, bestFaceNumber);
        if (correspondingCellIndex == IndexType(-1) ||
          correspondingFaceNumber == IndexType(-1) ||
          volumeMesh.cellMediumParameters[correspondingCellIndex].destroyed) continue;

        // check if it is last non-destructed edge for this cell
        IndexType incidentCellsCount = 0;
        for (IndexType faceNumber = 0; faceNumber < Space::FacesPerCell; ++faceNumber)
        {
          if (volumeMesh.GetCorrespondingCellIndex(cellIndex, faceNumber) != IndexType(-1) &&
            volumeMesh.GetCorrespondingFaceNumber(cellIndex, faceNumber) != IndexType(-1))
          {
            incidentCellsCount++;
          }
        }

        if (incidentCellsCount == 1)
        {
          DestroyCellMaterial(cellIndex, volumeMesh.cellMediumParameters[cellIndex].plasticity.powderShearMult);
          (*isCellBroken)[cellIndex] = true;
          continue;
        }

        if (volumeMesh.isCellAvailable[cellIndex])
        {
          //  volumeMesh.AddToAABBTree(cellIndex);
        }

        if (volumeMesh.isCellAvailable[correspondingCellIndex])
        {
          //  volumeMesh.AddToAABBTree(correspondingCellIndex);
        }

        DestroyFace(cellIndex, bestFaceNumber, dynamicBoundaryTypes[bestFaceNumber]);
      }
    }
  }

  HandleMaterialErosion();
}
