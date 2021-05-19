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
#include "DrReportSystem.hh"
#include "DrDefinitions.hh"
#include "DrBreakTable.hh"
#include <G4Scheduler.hh>
#include <G4SystemOfUnits.hh>
#include <G4MoleculeFinder.hh>

#ifndef State
#define State(theXInfo) (GetState<ReportSystemState>()->theXInfo)
#endif

DrReportSystem::DrReportSystem() : G4VITDiscreteProcess("DrReportSystem", fUserDefined), fVerbose(0) {
    // meaning this class contains a class inheriting from G4ProcessState
    G4VITProcess::SetInstantiateProcessState(false);
    pParticleChange = &aParticleChange;
    enableAlongStepDoIt = false;
    enableAtRestDoIt = false;
    enablePostStepDoIt = true;
    SetProcessSubType(903);
    fProposesTimeStep = true;
}

DrReportSystem::~DrReportSystem() {}

DrReportSystem::DrReportSystem(const DrReportSystem &right) : G4VITDiscreteProcess("ReportSystem", fUserDefined) {
    fVerbose = right.fVerbose;
}

DrReportSystem::ReportSystemState::ReportSystemState() : G4ProcessState() {
    fPreviousTimeAtPreStepPoint = -1;
}

void DrReportSystem::StartTracking(G4Track *track) {
    G4VProcess::StartTracking(track);
    G4VITProcess::fpState.reset(new ReportSystemState());
    G4VITProcess::StartTracking(track);
}

G4bool DrReportSystem::IsApplicable(const G4ParticleDefinition &aParticleType) {
    if (aParticleType.GetParticleName() == "Clock") return true;
    else return false;
}

G4double DrReportSystem::GetMeanFreePath(const G4Track &, G4double, G4ForceCondition *) {
    G4int currentBin = DrDefinitions::Instance()->GetCurrentExplicitBinNuber();
    G4double value = DrDefinitions::Instance()->GetExplicitBinDifferences()[currentBin];
    return value;
}

G4double DrReportSystem::PostStepGetPhysicalInteractionLength(
        const G4Track &track, G4double, G4ForceCondition *condition) {
    // form:
    // http://geant4.slac.stanford.edu/UsersWorkshop/PDF/Marc/AddingNewProcess.pdf
    *condition = NotForced;
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

    //@@@@ This is to get the length of the previous time step
    G4double previousTimeStep(-1.);
    if (State(fPreviousTimeAtPreStepPoint) != -1) {
        previousTimeStep = time - State(fPreviousTimeAtPreStepPoint);
    }

    State(fPreviousTimeAtPreStepPoint) = time;

    if ((previousTimeStep < 0.0) || (fpState->theNumberOfInteractionLengthLeft < 0.0)) {
        //@@@@ First step or just after proc
        fpState->theNumberOfInteractionLengthLeft = 1.0;
        DrDefinitions::Instance()->IncrimentCurrentExplicitBinNuber();
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
               << "DrReportSystem cannot be negative!" << G4endl;
        exit(EXIT_FAILURE);
    }

    //@@@@ negative lets the stepper know we are returning a time
    return value * -1;
}

G4VParticleChange *DrReportSystem::RecordSystem(const G4Track &) {

    std::map<const G4MoleculeDefinition*, G4int > moleculeCount;
    auto bTable = DrBreakTable::Instance();
    G4int bleachedCount{0};

    for(auto read: DrDefinitions::Instance()->GetNameMap()){
        moleculeCount[read.second] = 0;
    }

    G4TrackManyList *mainList = G4ITTrackHolder::Instance()->GetMainList();
    G4TrackManyList::iterator it = mainList->begin();
    G4TrackManyList::iterator end = mainList->end();

    for (; it != end; ++it) {
        G4Track *track = *it;
        const G4MoleculeDefinition *molDef = GetMolecule(track)->GetDefinition();

        if (moleculeCount.find(molDef) == moleculeCount.end()) {} // not found molecule, it is not one defined in pathway.in
        else moleculeCount[molDef]++; // found molecule, add one to count

        if(molDef->GetName().substr(0,3) == "DSB"){
            if(bTable->GetBreakMolecule(*track, "DrReportSystem")->sBreakEndA->fPKcsIsBleached) bleachedCount++;
            if(bTable->GetBreakMolecule(*track, "DrReportSystem")->sBreakEndB->fPKcsIsBleached) bleachedCount++;
        }
    }
    //----------------------------------------------------------------------
    // Add the recorded moleculeCount above into TypeTable
    //----------------------------------------------------------------------
    G4double time = round(G4Scheduler::Instance()->GetGlobalTime()/s);
    G4int initialDSBs = DrBreakTable::Instance()->fInitialBreakNumber;
    bTable->UpdateTypeTracking(moleculeCount,time);
    bTable->fBleachedStoreRun[time].push_back(G4double(bleachedCount)/G4double(initialDSBs*2));

    return &aParticleChange;
}

G4VParticleChange *DrReportSystem::PostStepDoIt(const G4Track &track, const G4Step &) {
    //@@@@ Clean the number of interaction length
    //@@@@ left when the process triggers
    ClearNumberOfInteractionLengthLeft();

    return RecordSystem(track);
}
