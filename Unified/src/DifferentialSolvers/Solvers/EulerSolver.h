#include "../DifferentialSolver.h"

template<typename Scalar>
class EulerSolver : public DifferentialSolver<Scalar>
{
public:
  using DifferentialSolver<Scalar>::timeStep;
  using DifferentialSolver<Scalar>::currTime;

  EulerSolver(): DifferentialSolver<Scalar>() {}

  void SetSystem(DifferentialSystem<Scalar>* system)
  {
    this->system = system;
    currCoords  = new Scalar[system->GetMaxDimentionsCount()];
    oldCoords   = (system->GetHierarchyLevelsCount() > 1) ? new Scalar[system->GetMaxDimentionsCount()] : 0;
    nextCoords  = new Scalar[system->GetMaxDimentionsCount()];
    derivatives = new Scalar[system->GetMaxDimentionsCount()];
    currTime = 0;
  }

  ~EulerSolver()
  {
    if (system->GetHierarchyLevelsCount() > 1)
    {
      delete [] oldCoords;
    }
    delete [] currCoords;
    delete [] nextCoords;
    delete [] derivatives;
  }

  int GetPhasesCount()
  {
    return 1;
  }

  void InitStep(Scalar timeStep, Scalar tolerance, bool updateInitialCoords)
  {
    DifferentialSolver<Scalar>::InitStep(timeStep, tolerance, updateInitialCoords);
    if (system->GetHierarchyLevelsCount() == 1)
    {
      system->GetCurrCoords(currTime, currCoords);
    }
  }

  void InitStep(const SolverState& solverState)
  {
    if (system->GetHierarchyLevelsCount() > 1)
    {
      system->GetCurrCoords(currTime, currCoords, oldCoords, solverState);
    }
  }

  bool AdvancePhase(const SolverState& solverState)
  {
    Scalar currStep = timeStep * (1 << solverState.hierarchyLevel);
    system->GetCurrDerivatives(derivatives, solverState);

    #pragma omp parallel for
    for(int coordIndex = 0; coordIndex < system->GetDimentionsCount(solverState); coordIndex++)
    {
      nextCoords[coordIndex] = currCoords[coordIndex] + derivatives[coordIndex] * currStep;
    }
    return true;
  }

  void AdvanceStep(const SolverState& solverState)
  {
    if (solverState.IsPreInitial())
    {
      currTime += timeStep;
    }

    if (system->GetHierarchyLevelsCount() > 1)
    {
      system->SetCurrCoords(currTime, nextCoords, oldCoords, solverState);
    } else
    {
      system->SetCurrCoords(currTime, nextCoords);
    }
  }

  void RevertStep(Scalar)
  {
    // never will be called
  }

  Scalar  GetLastStepError()
  {
    return Scalar(1.0);
  }

  Scalar GetTimeStepPrediction()
  {
    return timeStep;
  }
 
private:
  Scalar* currCoords;
  Scalar* oldCoords;
  Scalar* nextCoords;
  Scalar* derivatives;

  DifferentialSystem<Scalar> *system;
};
