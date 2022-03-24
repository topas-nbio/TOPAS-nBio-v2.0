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
#include "DrFRAP.hh"
#include "DrClock.hh"
#include "DrDefinitions.hh"
#include "TsParameterManager.hh"
#include "DrBreakMolecule.hh"
#include <G4MoleculeFinder.hh>
#include <G4Scheduler.hh>
#include <G4SystemOfUnits.hh>


#ifndef State
#define State(theXInfo) (GetState<FRAPState>()->theXInfo)
#endif

DrFRAP::DrFRAP() : G4VITDiscreteProcess("DrFRAP", fUserDefined), fVerbose(0) {
  // meaning this class contains a class inheriting from G4ProcessState
  G4VITProcess::SetInstantiateProcessState(false);
  SetProcessSubType(902);
  enableAlongStepDoIt = false;
  enablePostStepDoIt = true;
  enableAtRestDoIt = false;
  fProposesTimeStep = true;
  pParticleChange = &aParticleChange;
}

DrFRAP::~DrFRAP() {}

DrFRAP::DrFRAP(const DrFRAP &right)
    : G4VITDiscreteProcess("FRAP", fUserDefined) {
  fVerbose = right.fVerbose;
}

DrFRAP::FRAPState::FRAPState() : G4ProcessState() {
  fPreviousTimeAtPreStepPoint = -1;
}

void DrFRAP::StartTracking(G4Track *track) {
  G4VProcess::StartTracking(track);
  G4VITProcess::fpState.reset(new FRAPState());
  G4VITProcess::StartTracking(track);
}

G4bool DrFRAP::IsApplicable(const G4ParticleDefinition &aParticleType) {
  if(aParticleType == *DrClock::Definition()) return true;
  else return false;
}

G4double DrFRAP::GetMeanFreePath(const G4Track &/*track*/, G4double, G4ForceCondition *) {
  G4double bleachAtTime = 30 * second;
  G4double globalTime = G4Scheduler::Instance()->GetGlobalTime();
  if (globalTime > bleachAtTime) return DBL_MAX;
  else return bleachAtTime;
}

G4double DrFRAP::PostStepGetPhysicalInteractionLength(const G4Track &track, G4double, G4ForceCondition *condition) {

  *condition = NotForced;

  //@@@@ This is to get the length of the previous time step
  G4double previousTimeStep(-1.);
  //@@@@ -----------------------------------------------------------
  //@@@@ setting time with G4Scheduler GlobalTime rather than with
  //@@@@ track.GetGlobalTime() as is the case in similar code by M.K.. This is
  //@@@@ done because at larger time scales the track global time diverges
  //@@@@ from the scheduler time. Have been unable to determine why. Interest-
  //@@@@ ingly it seems to incrimentally grow in deviation size by the same
  //@@@@ .something ps amount across all simulations.
  G4double time = G4Scheduler::Instance()->GetGlobalTime();
  //@@@@ Should be:
  //G4double time = track.GetGlobalTime();
  //@@@@ But this throws exception ITStepManager015
  //@@@@ -----------------------------------------------------------

  if (State(fPreviousTimeAtPreStepPoint) != -1) {
    previousTimeStep = time - State(fPreviousTimeAtPreStepPoint);
  }

  State(fPreviousTimeAtPreStepPoint) = time;

  if ((previousTimeStep < 0.0) || (fpState->theNumberOfInteractionLengthLeft < 0.0)) {
    fpState->theNumberOfInteractionLengthLeft = 1.0;
  }

  else if (previousTimeStep > 0.0) {
    SubtractNumberOfInteractionLengthLeft(previousTimeStep);
    if (fpState->theNumberOfInteractionLengthLeft <= perMillion) {
      fpState->theNumberOfInteractionLengthLeft = 0.0;
    }
  }

  fpState->currentInteractionLength = GetMeanFreePath(track, previousTimeStep, condition);

  G4double value = fpState->theNumberOfInteractionLengthLeft * fpState->currentInteractionLength;

  if (value < 0) {
    G4cerr << "ERROR: Post step interaction length returned by"
           << "DrFRAP cannot be negative!" << G4endl;
      DrDefinitions::Instance()->GetParameterManager()->AbortSession(1);
  }

  //@@@@ negative lets the stepper know we are returning a time
  return value * -1;
}

G4VParticleChange *DrFRAP::FRAP(const G4Track &track) {
  aParticleChange.Clear();
  aParticleChange.Initialize(track);

  auto listOfMoleculesWithPKcs = DrDefinitions::Instance()->GetHasPKcs();

  G4TrackManyList *mainList = G4ITTrackHolder::Instance()->GetMainList();
  G4TrackManyList::iterator it = mainList->begin();
  G4TrackManyList::iterator end = mainList->end();

  for (; it != end; ++it) {
    G4Track *molTrack = *it;
    G4Molecule *molecule = GetMolecule(molTrack);
    const G4MoleculeDefinition *moleculeDefinition = molecule->GetDefinition();
    const G4String& moleculeName = moleculeDefinition->GetName();

    for(const auto& NameOfMoleculesWithPKc: listOfMoleculesWithPKcs){
      auto listOfAllMolecules = DrDefinitions::Instance()->GetNameMap();
      if(listOfAllMolecules.find(NameOfMoleculesWithPKc) != listOfAllMolecules.end()){
        if(moleculeName == NameOfMoleculesWithPKc){
            DrBreakMolecule* breakMolecule = (DrBreakMolecule*)(molTrack->GetAuxiliaryTrackInformation(G4PhysicsModelCatalog::GetIndex("DrBreakMolecule")));
          breakMolecule->sBreakEndA->fPKcsIsBleached = true;
          if(breakMolecule->sBreakEndB->fOriginalBreakMoleculeID != -1) breakMolecule->sBreakEndB->fPKcsIsBleached = true;
        }
      }
    }
  }
  //--------------------------------------------------------------------------
  return &aParticleChange;
}

G4VParticleChange *DrFRAP::PostStepDoIt(const G4Track &track, const G4Step &) {
    //@@@@ Clean the number of interaction length
    //@@@@ left when the process triggers
    ClearNumberOfInteractionLengthLeft();
    return FRAP(track);
}
