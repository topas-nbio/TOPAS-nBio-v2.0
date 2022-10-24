// Extra Class for use by DaMaRiS
//
// ********************************************************************
// *                                                                  *
// * This file is part of the TOPAS-nBio extensions to the            *
// *   TOPAS Simulation Toolkit.                                      *
// * The TOPAS-nBio extensions are freely available under the license *
// *   agreement set forth at: https://topas-nbio.readthedocs.io/     *
// *                                                                  *
// ********************************************************************
//
// Created by John-William Warmenhoven.
// DaMaRiS is developed at the University of Manchester.
// See README for references.
//
#include "DrWaitingTime.hh"
#include "DrDefinitions.hh"
#include "TsParameterManager.hh"
#include "DrPrecompiler.hh"
#include "DrBreakMolecule.hh"
#include <G4SystemOfUnits.hh>

#ifndef State
#define State(theXInfo) (GetState<WaitingTimeState>()->theXInfo)
#endif

DrWaitingTime::DrWaitingTime() : G4VITDiscreteProcess("DrWaitingTime", fUserDefined) {
  // meaning this class contains a class inheriting from G4ProcessState
  G4VITProcess::SetInstantiateProcessState(false);
  SetProcessSubType(904);
  enableAlongStepDoIt = false;
  enablePostStepDoIt = true;
  enableAtRestDoIt = false;
  fProposesTimeStep = true;
  pParticleChange = &aParticleChange;
}

DrWaitingTime::~DrWaitingTime() {}

DrWaitingTime::WaitingTimeState::WaitingTimeState() : G4ProcessState() {
  fPreviousTimeAtPreStepPoint = -1;
}

void DrWaitingTime::StartTracking(G4Track *track) {
  G4VProcess::StartTracking(track);
  G4VITProcess::fpState.reset(new WaitingTimeState());
  G4VITProcess::StartTracking(track);
}

G4bool DrWaitingTime::IsApplicable(const G4ParticleDefinition &aParticleType) {
  if (aParticleType.GetParticleName().substr(0, 3) != "DSB") return false;
  else return true;
}

G4double DrWaitingTime::GetMeanFreePath(const G4Track &, G4double, G4ForceCondition *) {
  G4double value = DBL_MAX * (-1.);
  return value;
}

G4double DrWaitingTime::PostStepGetPhysicalInteractionLength(const G4Track &track, G4double, G4ForceCondition *condition) {
  *condition = NotForced;

  //@@@@ This is to get the length of the previous time step
  G4double previousTimeStep(-1.);
  if (State(fPreviousTimeAtPreStepPoint) != -1.) {
    previousTimeStep = track.GetGlobalTime() - State(fPreviousTimeAtPreStepPoint);
  }

  State(fPreviousTimeAtPreStepPoint) = track.GetGlobalTime();

  DrBreakMolecule* breakMolecule = (DrBreakMolecule*)(track.GetAuxiliaryTrackInformation(G4PhysicsModelCatalog::GetIndex("DrBreakMolecule")));

  G4double moleculeWaitingTime = breakMolecule->fWaitingTime;
  G4double moleculeDiffusionTime = breakMolecule->fDiffusionTime;
  G4bool isWaiting = breakMolecule->fIsWaiting;
  G4double returnValue{-1};

  //@@@@ This triggers if this is the first step in the simulation (prevTimeStep
  //@@@@ < 0) or if this is the first time the molecule passes through this
  //@@@@ class (moleculeWaitingTime and moleculeDiffusionTime are both
  //@@@@ initialised as -1).
  if (previousTimeStep < 0 || (moleculeWaitingTime < 0 && moleculeDiffusionTime < 0)) {
    GenerateWaitingTime(track);
    breakMolecule->fIsWaiting = true;
    returnValue = breakMolecule->fWaitingTime;
  }
  //@@@@ This triggers if the molecule has just switched from a diffusive state
  //@@@@ to a trapped state in the previous step
  else if (isWaiting && State(fHasJustProc) == true) {
    State(fHasJustProc) = false;
    returnValue = moleculeWaitingTime;
  }
  //@@@@ This triggers if the molecule has just switched from a trapped state to
  //@@@@ a diffusive state in the previous step
  else if (!isWaiting && State(fHasJustProc) == true) {
    State(fHasJustProc) = false;
    returnValue = moleculeDiffusionTime;
  }
  //@@@@ This triggers if the molecule was in a trapped state last step as well
  //@@@@ and needs to update it's trapped time by subtracting the current time
  //@@@@ step.
  else if (isWaiting && State(fHasJustProc) == false) {
    G4double tempReturnValue = moleculeWaitingTime - previousTimeStep;
    if (tempReturnValue < perMillion) {
      tempReturnValue = 0;
    }
    breakMolecule->fWaitingTime = tempReturnValue;
    returnValue = tempReturnValue;
  }
  //@@@@ This triggers if the molecule was in a diffusive state last step as
  //@@@@ well and needs to update it's diffusion time by subtracting the current
  //@@@@ time step. Because of the small diffusion time step this may cause
  //@@@@ returning of negative time steps therefore a safety has been built in
  //@@@@ here and above to ensure that negative times get returned as 0 instead
  //@@@@ (meaning the process is late but will now immediately execute).
  else if (!isWaiting && State(fHasJustProc) == false) {
    G4double tempReturnValue = moleculeDiffusionTime - previousTimeStep;
    if (tempReturnValue <= perMillion) {
      tempReturnValue = 0;
    }
    breakMolecule->fDiffusionTime = tempReturnValue;
    returnValue = tempReturnValue;
  }
  else {
    G4cerr << "ERROR: DrWaitingTime parameters not as expected PSGPIL" << G4endl;
      DrDefinitions::Instance()->GetParameterManager()->AbortSession(1);
  }

  fpState->currentInteractionLength = returnValue;

  G4double value = fpState->currentInteractionLength;

  if (value < 0) {
    G4cerr << "ERROR: Post step interaction length returned by"
           << "DrWaitingTime cannot be negative!" << G4endl;
      DrDefinitions::Instance()->GetParameterManager()->AbortSession(1);
  }

#ifdef DEBUG_DAMARIS
DrBreakTable::Instance()->diffTimeStore.push_back(value);
#endif

  //@@@@ negative lets the stepper know we are returning a time
  return value * -1.0;
}

G4VParticleChange *DrWaitingTime::PostStepDoIt(const G4Track &track, const G4Step &) {
  if ((track.GetTrackStatus() == fStopButAlive) || (track.GetTrackStatus() == fStopAndKill)) {
    aParticleChange.Initialize(track);
    return &aParticleChange;
  }
  else {
    //@@@@ Clean the number of interaction length
    //@@@@ left when the process triggers
    ClearNumberOfInteractionLengthLeft();
    aParticleChange.Initialize(track);
    Process(track);
    return &aParticleChange;
  }
}

void DrWaitingTime::Process(const G4Track &track) {
    DrBreakMolecule* breakMolecule = (DrBreakMolecule*)(track.GetAuxiliaryTrackInformation(G4PhysicsModelCatalog::GetIndex("DrBreakMolecule")));

    G4bool isWaiting = breakMolecule->fIsWaiting;
    if (!isWaiting) {
        breakMolecule->fDiffusionTime = -1.0;
        GenerateWaitingTime(track);
        breakMolecule->fIsWaiting = true;
        State(fHasJustProc) = true;
    }
    else {
        breakMolecule->fDiffusionTime = 1.0 * picosecond;
        breakMolecule->fIsWaiting = false;
        breakMolecule->fWaitingTime = -1.0;
        State(fHasJustProc) = true;
    }
}

void DrWaitingTime::GenerateWaitingTime(const G4Track &track) {
  G4double A = DrDefinitions::Instance()->GetMinumumWaitingTime();
  G4double rand = G4UniformRand();
  G4double alpha = 0.5;
  G4double waitingTime = A / pow((1 - rand), (1 / alpha));
  DrBreakMolecule* breakMolecule = (DrBreakMolecule*)(track.GetAuxiliaryTrackInformation(G4PhysicsModelCatalog::GetIndex("DrBreakMolecule")));
  breakMolecule->fWaitingTime = waitingTime;
}
