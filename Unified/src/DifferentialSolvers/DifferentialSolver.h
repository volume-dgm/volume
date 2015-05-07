#pragma once
#include "SolverState.h"
#include "DifferentialSystem.h"
#include <limits>

template<typename Scalar>
class DifferentialSolver
{
public:
  DifferentialSolver() {}
  virtual ~DifferentialSolver() {}

  virtual void    SetSystem(DifferentialSystem<Scalar>* system) = 0;
  virtual int     GetPhasesCount() = 0;

  virtual void    InitStep(Scalar timeStep, Scalar tolerance, bool)
  {
    this->timeStep    = timeStep;
    this->tolerance   = tolerance;
    predictedStep = std::numeric_limits<Scalar>::max();
  }

  virtual void    InitStep(const SolverState&) = 0;
  virtual bool    AdvancePhase(const SolverState&) = 0;
  virtual void    AdvanceStep(const SolverState&) = 0;
  virtual void    RevertStep(Scalar currTime) = 0;

  virtual Scalar  GetLastStepError() = 0;
  virtual Scalar  GetTimeStepPrediction() = 0;

  Scalar GetCurrTime() const
  {
    return currTime;
  }

protected:
  Scalar currTime;
  Scalar timeStep;
  Scalar tolerance;
  Scalar stepError;
  Scalar predictedStep;
};