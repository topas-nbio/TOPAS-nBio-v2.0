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
#include "DrReportMSD.hh"
#include "DrDefinitions.hh"
#include "DrClock.hh"
#include "DrBreakTable.hh"
#include <G4Scheduler.hh>
#include <G4UnitsTable.hh>
#include <G4SystemOfUnits.hh>
#include <G4MoleculeFinder.hh>

#ifndef State
#define State(theXInfo) (GetState<ReportMSDState>()->theXInfo)
#endif

DrReportMSD::DrReportMSD() : G4VITDiscreteProcess("ReportMSD", fUserDefined), fVerbose(0) {
  // meaning this class contains a class inheriting from G4ProcessState
  G4VITProcess::SetInstantiateProcessState(false);
  pParticleChange = &aParticleChange;
  enableAlongStepDoIt = false;
  enableAtRestDoIt = false;
  enablePostStepDoIt = true;
  SetProcessSubType(905);
  fProposesTimeStep = true;
}

DrReportMSD::~DrReportMSD() {}

DrReportMSD::DrReportMSD(const DrReportMSD &right) : G4VITDiscreteProcess("ReportMSD", fUserDefined) {
  fVerbose = right.fVerbose;
}

DrReportMSD::ReportMSDState::ReportMSDState() : G4ProcessState() {
  fPreviousTimeAtPreStepPoint = -1;
}

void DrReportMSD::StartTracking(G4Track *track) {
  G4VProcess::StartTracking(track);
  G4VITProcess::fpState.reset(new ReportMSDState());
  G4VITProcess::StartTracking(track);
}

G4bool DrReportMSD::IsApplicable(const G4ParticleDefinition &aParticleType) {
  if(aParticleType == *DrClock::Definition()) return true;
  if(aParticleType.GetParticleName() == "Clock") return true;
  else return false;
}

G4double DrReportMSD::GetMeanFreePath(const G4Track &, G4double, G4ForceCondition *) {
  G4double MSDRecordingDuration = DrDefinitions::Instance()->GetObserveDurationForMSD();
  G4double MSDRecordingStepSize = DrDefinitions::Instance()->GetObserveStepSizeForMSD();
  G4double GlobalTime = G4Scheduler::Instance()->GetGlobalTime();
  if (GlobalTime <= MSDRecordingDuration) return MSDRecordingStepSize;
  else return DBL_MAX;
}

G4double DrReportMSD::PostStepGetPhysicalInteractionLength(
        const G4Track &track, G4double, G4ForceCondition *condition) {
  // form:
  // http://geant4.slac.stanford.edu/UsersWorkshop/PDF/Marc/AddingNewProcess.pdf
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

  if (State(fPreviousTimeAtPreStepPoint) != -1.) {

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
           << "DrReportMSD cannot be negative!" << G4endl;
    exit(EXIT_FAILURE);
  }

  //@@@@ negative lets the stepper know we are returning a time
  return value * -1;
}

G4VParticleChange *DrReportMSD::RecordMSD(const G4Track &thisTrack) {
  aParticleChange.Clear();


  //@@@@ -----------------------------------------------------------
  //@@@@ setting time with G4Scheduler GlobalTime rather than with
  //@@@@ track.GetGlobalTime() as is the case in similar code by M.K.. This is
  //@@@@ done because at larger time scales the track global time diverges
  //@@@@ from the scheduler time. Have been unable to determine why. Interest-
  //@@@@ ingly it seems to incrimentally grow in deviation size by the same
  //@@@@ .something ps amount across all simulations.
  G4double time = G4Scheduler::Instance()->GetGlobalTime();
  //@@@@ Should be:
  // G4double time = thisTrack.GetGlobalTime();
  //@@@@ But this throws exception ITStepManager015
  //@@@@ -----------------------------------------------------------

  //-----------------------------verbose--------------------------------------
  if (fVerbose > 1) {
    G4cout<< "tick @ "<< std::setprecision(4)
          << G4Scheduler::Instance()->GetGlobalTime()/s
          << " s from molecule ID = "
          << thisTrack.GetTrackID() << G4endl;
  }
  //--------------------------------------------------------------------------

  //@@@@ Get track list for checks below
  DrBreakTable *bTable = DrBreakTable::Instance();

  if(fVerbose > 1){
    G4cout<<"Storing displacement for "<<G4BestUnit(time,"Time") <<G4endl;
  }
  G4TrackManyList *mainList = G4ITTrackHolder::Instance()->GetMainList();
  G4TrackManyList::iterator it = mainList->begin();
  G4TrackManyList::iterator end = mainList->end();

  for (; it != end; ++it) {
    G4Track *track = *it;
    G4Molecule *molecule = GetMolecule(track);
    const G4MoleculeDefinition *moleculeDefinition = molecule->GetDefinition();
    G4String name = moleculeDefinition->GetName();
    if (name.substr(0, 3) == "DSB") {

      DrBreakMolecule *breakMolecule = bTable->GetBreakMolecule(*track, "UserPTSA");

      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      //@@@@ This calculates the displacement of every particle.
      //@@@@ Displacements are pushed onto a static container,
      //@@@@ used in calculating the MSD at end of simulations.
      G4ThreeVector startPosition1;
      G4ThreeVector startPosition2;
      G4double displacement1{0.};
      G4double displacement2{0.};
      //----------------------------------------------------------------------
      //@@@@ MSD Calculation
      if(G4Scheduler::Instance()->GetGlobalTime() < DrDefinitions::Instance()->GetObserveDurationForMSD()){
        if (moleculeDefinition->GetName().substr(0, 6) == "DSBEnd" || moleculeDefinition->GetName().substr(0, 12) == "DSB_Fixed_HR") {

          startPosition1 = breakMolecule->sBreakEndA->fOriginalPosition;
          displacement1 = (startPosition1 - track->GetPosition()).mag() / nm;
          if(fVerbose > 3){
            G4cout
                    <<"Start("<<startPosition1/nm<<") | "
                    <<"Position("<<track->GetPosition()/nm<<") |"
                    <<"Displacement("<<displacement1<<")"
                    <<G4endl;
          }
          bTable->fMSDTrackingStore[round(time/s)].push_back(displacement1);
        }
        else if (moleculeDefinition->GetName().substr(0, 6) == "DSBSyn" || moleculeDefinition->GetName().substr(0, 6) == "DSB_Fi") {
          startPosition1 = breakMolecule->sBreakEndA->fOriginalPosition;
          startPosition2 = breakMolecule->sBreakEndB->fOriginalPosition;
          displacement1 = (startPosition1 - track->GetPosition()).mag() / nm;
          displacement2 = (startPosition2 - track->GetPosition()).mag() / nm;

          if(fVerbose > 3){
            G4cout
                    <<"Start1("<<startPosition1/nm<<") | "
                    <<"Position1("<<track->GetPosition()/nm<<") |"
                    <<"Displacement1("<<displacement1<<")"
                    <<G4endl
                    <<"Start2("<<startPosition2/nm<<") | "
                    <<"Position2("<<track->GetPosition()/nm<<") |"
                    <<"Displacement2("<<displacement2<<")"
                    <<G4endl;
          }
          bTable->fMSDTrackingStore[round(time/s)].push_back(displacement1);
          bTable->fMSDTrackingStore[round(time/s)].push_back(displacement2);
        }
      }
    }
  }
  return &aParticleChange;
}

G4VParticleChange *DrReportMSD::PostStepDoIt(const G4Track &track, const G4Step &) {
  //@@@@ Clean the number of interaction length
  //@@@@ left when the process triggers
  ClearNumberOfInteractionLengthLeft();
  return RecordMSD(track);
}
