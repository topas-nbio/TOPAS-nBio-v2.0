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
#include <CLHEP/Random/Stat.h>

#include "DrTransportation_Subdiffusion_CTRW.hh"

#include <G4DNAMolecularMaterial.hh>
#include <G4ITSafetyHelper.hh>
#include <G4NistManager.hh>
#include <G4RandomDirection.hh>
#include <G4SystemOfUnits.hh>
#include <G4TrackingInformation.hh>
#include <G4UnitsTable.hh>
#include <G4Scheduler.hh>
#include <G4MoleculeFinder.hh>
#include "DrBreakTable.hh"
#include "DrDefinitions.hh"
#include "DrPrecompiler.hh"

using namespace std;

#ifndef State
#define State(theXInfo) (GetState<DrITSubdiffusionState>()->theXInfo)
#endif

#ifdef USE_COLOR
#define RED "\033[0;31m"
#define LIGHT_RED "\33[1;31m"
#define GREEN "\033[32;40m"
#define GREEN_ON_BLUE "\033[1;32;44m"
#define RESET_COLOR "\033[0m"
#else
#define RED ""
#define LIGHT_RED ""
#define GREEN ""
#define GREEN_ON_BLUE ""
#define RESET_COLOR ""
#endif

#ifdef DEBUG_MEM
#include "G4MemStat.hh"
using namespace G4MemStat;
using G4MemStat::MemStat;
#endif

static double InvErf(double x) { return CLHEP::HepStat::inverseErf(x); }

static double InvErfc(double x) { return CLHEP::HepStat::inverseErf(1. - x); }

static double Erfc(double x) { return 1 - CLHEP::HepStat::erf(1. - x); }

#ifndef State
#define State(theXInfo) (GetState<G4ITTransportationState>()->theXInfo)
#endif

DrTransportation_Subdiffusion_CTRW::DrTransportation_Subdiffusion_CTRW(const G4String &aName, G4int verbosity)
        : G4ITTransportation(aName, verbosity) {

    fVerboseLevel = 0;

    fpState.reset(new DrITSubdiffusionState());

    // ctor
    SetProcessSubType(92);

    fNistWater = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
    fpWaterDensity = nullptr;

    fUseMaximumTimeBeforeReachingBoundary = true;
    fUseSchedulerMinTimeSteps = false;
    fSpeedMeUp = true;

    fInternalMinTimeStep = 1 * picosecond;
    fpSubDiffusiveAction = nullptr;

    //@@@@
    fTrappedDiffCoef = DrDefinitions::Instance()->GetTrapDiff();
    fJumpDiffCoef = DrDefinitions::Instance()->GetJumpDiff();
    //@@@@
}

DrTransportation_Subdiffusion_CTRW::~DrTransportation_Subdiffusion_CTRW() {
    if (fpSubDiffusiveAction) delete fpSubDiffusiveAction;
}

DrTransportation_Subdiffusion_CTRW::DrTransportation_Subdiffusion_CTRW(const DrTransportation_Subdiffusion_CTRW &right)
        : G4ITTransportation(right) {
    // copy ctor
    SetProcessSubType(61);
    fUseMaximumTimeBeforeReachingBoundary = right.fUseMaximumTimeBeforeReachingBoundary;
    fUseSchedulerMinTimeSteps = right.fUseSchedulerMinTimeSteps;
    fNistWater = right.fNistWater;
    fpWaterDensity = right.fpWaterDensity;
    fInternalMinTimeStep = right.fInternalMinTimeStep;
    fSpeedMeUp = right.fSpeedMeUp;
    fpSubDiffusiveAction = right.fpSubDiffusiveAction;
}

DrTransportation_Subdiffusion_CTRW &DrTransportation_Subdiffusion_CTRW:: operator=(const DrTransportation_Subdiffusion_CTRW &rhs) {
    if (this == &rhs) return *this; // handle self assignment
    // assignment operator
    return *this;
}

DrTransportation_Subdiffusion_CTRW::DrITSubdiffusionState::DrITSubdiffusionState()
        : G4ITTransportationState() {
    fPathLengthWasCorrected = false;
    fTimeStepReachedLimit = false;
    fComputeLastPosition = false;
    fRandomNumber = -1;
}

void DrTransportation_Subdiffusion_CTRW::StartTracking(G4Track *track) {
    fpState.reset(new DrITSubdiffusionState());
    SetInstantiateProcessState(false);
    G4ITTransportation::StartTracking(track);
}

void DrTransportation_Subdiffusion_CTRW::BuildPhysicsTable(
        const G4ParticleDefinition &particle) {
    if (verboseLevel > 0) {
        G4cout << G4endl << GetProcessName() << ":   for  " << setw(24)
               << particle.GetParticleName() << "\tSubType= " << GetProcessSubType()
               << G4endl;
    }
    // Initialize water density pointer
    fpWaterDensity = G4DNAMolecularMaterial::Instance()->GetDensityTableFor(G4Material::GetMaterial("G4_WATER"));

    fpSafetyHelper->InitialiseHelper();
    G4ITTransportation::BuildPhysicsTable(particle);
}

void DrTransportation_Subdiffusion_CTRW::ComputeStep(const G4Track &track, const G4Step &step, /*const*/ double timeStep, double &spaceStep) {
    /* If this method is called, this step
     * cannot be geometry limited.
     * In case the step is limited by the geometry,
     * this method should not be called.
     * The fTransportEndPosition calculated in
     * the method AlongStepIL should be taken
     * into account.
     * In order to do so, the flag IsLeadingStep
     * is on. Meaning : this track has the minimum
     * interaction length over all others.
     */

    if (GetIT(track)->GetTrackingInfo()->IsLeadingStep()) {
        const G4VITProcess *ITProc = ((const G4VITProcess *)step.GetPostStepPoint()->GetProcessDefinedStep());
        bool makeException = true;

        if (ITProc && ITProc->ProposesTimeStep()) makeException = false;

        if (makeException) {
            G4ExceptionDescription exceptionDescription;
            exceptionDescription << "ComputeStep is called while the track has"
                                    "the minimum interaction time";
            exceptionDescription << " so it should not recompute a timeStep ";
            G4Exception("DrTransportation_Subdiffusion_CTRW::ComputeStep",
                        "DrTransportation_Subdiffusion_CTRW001", FatalErrorInArgument,
                        exceptionDescription);
        }
    }

    State(fGeometryLimitedStep) = false;

    G4Molecule *molecule = GetMolecule(track);

    if (timeStep > 0) {
        spaceStep = DBL_MAX;

        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        G4bool isWaiting = DrBreakTable::Instance()->GetBreakMolecule(track, "Sub")->fIsWaiting;
        G4double diffCoeff = (isWaiting) ? fTrappedDiffCoef : fJumpDiffCoef;
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        static double sqrt_2 = sqrt(2.);
        double sqrt_Dt = sqrt(diffCoeff * timeStep);
        double sqrt_2Dt = sqrt_2 * sqrt_Dt;
        if (State(fTimeStepReachedLimit) == true) {
            //========================================================================
            State(fGeometryLimitedStep) = true; // important
            //========================================================================
            spaceStep = State(fEndPointDistance);
        } else {
            double x = G4RandGauss::shoot(0, sqrt_2Dt);
            double y = G4RandGauss::shoot(0, sqrt_2Dt);
            double z = G4RandGauss::shoot(0, sqrt_2Dt);
            spaceStep = sqrt(x * x + y * y + z * z);

            if (spaceStep >= State(fEndPointDistance)) {
                //======================================================================
                State(fGeometryLimitedStep) = true; // important
                //======================================================================
                if (fUseSchedulerMinTimeSteps == false) { // jump over barrier NOT used
#ifdef G4VERBOSE
                    if (fVerboseLevel > 1) {
                        G4cout << GREEN_ON_BLUE
                               << "DrTransportation_Subdiffusion_CTRW::ComputeStep() : "
                               << "Step was limited to boundary" << RESET_COLOR << G4endl;
                    }
#endif
                    if (State(fRandomNumber) >= 0) { // CDF is used
                        //==================================================================
                        // BE AWARE THAT THE TECHNIQUE USED BELOW IS A 1D APPROXIMATION
                        // Cumulative density function for the 3D case is not yet
                        // implemented
                        //==================================================================
                        double value = State(fRandomNumber) +
                                       (1 - State(fRandomNumber)) * G4UniformRand();
                        double invErfc = InvErfc(value);
                        spaceStep = invErfc * 2 * sqrt_Dt;
                        if (State(fTimeStepReachedLimit) == false) {
                            //================================================================
                            State(fGeometryLimitedStep) = false; // important
                            //================================================================
                        }
                    }
                    else if(fUseMaximumTimeBeforeReachingBoundary == false){//CDF is used
                        double min_randomNumber =
                                Erfc(State(fEndPointDistance) / 2 * sqrt_Dt);
                        double value =
                                min_randomNumber + (1 - min_randomNumber) * G4UniformRand();
                        double invErfc = InvErfc(value);
                        spaceStep = invErfc * 2 * sqrt_Dt;
                        if (spaceStep >= State(fEndPointDistance)) {
                            //================================================================
                            State(fGeometryLimitedStep) = true; // important
                            //================================================================
                        }
                        else if (State(fTimeStepReachedLimit) == false) {
                            //================================================================
                            State(fGeometryLimitedStep) = false; // important
                            //================================================================
                        }
                    }
                    else { // CDF is NOT used
                        //==================================================================
                        State(fGeometryLimitedStep) = true; // important
                        //==================================================================
                        spaceStep = State(fEndPointDistance);
                    }
                }

                State(fTransportEndPosition) =
                        spaceStep * track.GetMomentumDirection() + track.GetPosition();
            }
            else {
                //======================================================================
                State(fGeometryLimitedStep) = false; // important
                //======================================================================
                State(fTransportEndPosition) =
                        spaceStep * step.GetPostStepPoint()->GetMomentumDirection() +
                        track.GetPosition();
            }
        }
    } else {
        spaceStep = 0.;
        State(fTransportEndPosition) = track.GetPosition();
        State(fGeometryLimitedStep) = false;
    }

    State(fCandidateEndGlobalTime) =
            step.GetPreStepPoint()->GetGlobalTime() + timeStep;
    State(fEndGlobalTimeComputed) = true;

#ifdef G4VERBOSE
    if (fVerboseLevel > 1) {
        G4cout << GREEN_ON_BLUE
               << "DrTransportation_Subdiffusion_CTRW::ComputeStep() : "
               << " trackID : " << track.GetTrackID()
               << " : Molecule name: " << molecule->GetName() << G4endl
               << "Initial position:" << G4BestUnit(track.GetPosition(), "Length")
               << G4endl << "Initial direction:" << track.GetMomentumDirection()
               << G4endl << "Final position:"
               << G4BestUnit(State(fTransportEndPosition), "Length") << G4endl
               << "Initial magnitude:"
               << G4BestUnit(track.GetPosition().mag(), "Length") << G4endl
               << "Final magnitude:"
               << G4BestUnit(State(fTransportEndPosition).mag(), "Length") << G4endl
               << "Diffusion length : " << G4BestUnit(spaceStep, "Length")
               << " within time step : " << G4BestUnit(timeStep, "Time") << G4endl
               << "State(fTimeStepReachedLimit)= " << State(fTimeStepReachedLimit)
               << G4endl
               << "State(fGeometryLimitedStep)=" << State(fGeometryLimitedStep)
               << G4endl << "End point distance was: "
               << G4BestUnit(State(fEndPointDistance), "Length") << G4endl
               << RESET_COLOR << G4endl << G4endl;
    }
#endif

#ifdef DEBUG_DAMARIS
    if(spaceStep != 0.){
        DrBreakTable::Instance()->spaceStepStore.push_back(spaceStep);
    }
#endif
}

G4VParticleChange *
DrTransportation_Subdiffusion_CTRW::PostStepDoIt(const G4Track &track, const G4Step &step) {

    G4ITTransportation::PostStepDoIt(track, step);

    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    //@@@@ If the molecule is on the boundry of the nucleus it will
    //@@@@ get set to be destroyed by PostStepDoIt(...), before this
    //@@@@ gets done we create a new molecule on the boundry and
    //@@@@ copy over the relevant information. Later in
    //@@@@ Diffusion(...) the momentum direction of the new particle
    //@@@@ will be set towards the centre of the nucleus +/- 30 deg.
    if(fParticleChange.GetTrackStatus() == 2 && State(fGeometryLimitedStep)){
        //-------------------------Propose Changes----------------------------------
        G4Molecule* mother = GetMolecule(track);
        G4MoleculeDefinition *motherDef = const_cast<G4MoleculeDefinition *>(mother->GetDefinition());
        //@@@@ Create the daughter
        fParticleChange.SetNumberOfSecondaries(1);
        G4Molecule* theDaughterMolecule = new G4Molecule(motherDef);
        G4double trackTime = G4Scheduler::Instance()->GetGlobalTime();
        G4Track* theDaughterTrack = theDaughterMolecule->BuildTrack(trackTime, track.GetPosition());
        theDaughterTrack->SetTrackStatus(fAlive);
        fParticleChange.G4VParticleChange::AddSecondary(theDaughterTrack);
        //@@@@ Carries over the break information, the molecule  is changed but
        //@@@@ the break is still the same.
        G4int auxIndex = DrBreakTable::Instance()->fBreakMolAuxIndex;
        DrBreakMolecule* daughterBreak = (DrBreakMolecule*)track.GetAuxiliaryTrackInformation(auxIndex);
        track.SetAuxiliaryTrackInformation(auxIndex,nullptr);
        theDaughterTrack->SetAuxiliaryTrackInformation(auxIndex,daughterBreak);
        //--------------------------------------------------------------------------
    }
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#ifdef G4VERBOSE
    if (fVerboseLevel > 1) {
        G4cout << GREEN_ON_BLUE
               << "DrTransportation_Subdiffusion_CTRW::PostStepDoIt() :"
               << " trackID : " << track.GetTrackID()
               << " Molecule name: " << GetMolecule(track)->GetName() << G4endl;
        G4cout << "Diffusion length : "
               << G4BestUnit(step.GetStepLength(), "Length")
               << " within time step : " << G4BestUnit(step.GetDeltaTime(), "Time")
               << "\t Current global time : "
               << G4BestUnit(track.GetGlobalTime(), "Time") << RESET_COLOR << G4endl
               << G4endl;
    }
#endif

    return &fParticleChange;
}

void DrTransportation_Subdiffusion_CTRW::Diffusion(const G4Track &track) {
#ifdef DEBUG_MEM
    MemStat mem_first = MemoryUsage();
#endif

#ifdef G4VERBOSE
    if (fVerboseLevel > 1) {
        G4cout << GREEN_ON_BLUE << setw(18)
               << "DrTransportation_Subdiffusion_CTRW::Diffusion :" << setw(8)
               << GetIT(track)->GetName() << "\t trackID:" << track.GetTrackID()
               << "\t"
               << " Global Time = " << G4BestUnit(track.GetGlobalTime(), "Time")
               << RESET_COLOR << G4endl << G4endl;
    }
#endif

    G4Material *material = track.GetMaterial();

    G4double waterDensity = (*fpWaterDensity)[material->GetIndex()];

    if (waterDensity == 0.0) {
        if (fpSubDiffusiveAction) {
            // Let the user Brownian action class decide what to do
            fpSubDiffusiveAction->Transport(track, fParticleChange);
            return;
        } else {
#ifdef G4VERBOSE
            if (fVerboseLevel) {
                G4cout << "A track is outside water material : trackID = "
                       << track.GetTrackID() << " (" << GetMolecule(track)->GetName()
                       << ")" << G4endl;
                G4cout << "Local Time : " << G4BestUnit(track.GetGlobalTime(), "Time")
                       << G4endl;
                G4cout << "Step Number :" << track.GetCurrentStepNumber() << G4endl;
            }
#endif
            fParticleChange.ProposeEnergy(0.);
            fParticleChange.ProposeTrackStatus(fStopAndKill);
            return; // &fParticleChange is the final returned object
        }
    }

#ifdef DEBUG_MEM
    MemStat mem_intermediaire = MemoryUsage();
  mem_diff = mem_intermediaire - mem_first;
  G4cout
      << "\t\t\t >> || MEM || In DrTransportation_Subdiffusion_CTRW::Diffusion "
         "after dealing with waterDensity for "
      << track.GetTrackID() << ", diff is : " << mem_diff << G4endl;
#endif

    fParticleChange.ProposeMomentumDirection(G4RandomDirection());

    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    //@@@@ This code reflects the DSB molecule within +/- 30 degrees
    //@@@@ towards the nucleus centre if it is on the nucleus boundry
    if(track.GetPosition().mag() == DrDefinitions::Instance()->GetBoundingCellOrNucleusRadius()){
        G4ThreeVector unit = (-track.GetPosition()).unit();
        // first random angle between +30 and -30 degrees
        G4double angle1 = ((2 * G4UniformRand() - 1) * (CLHEP::pi / 6));
        // second random angle between +30 and -30 degrees
        G4double angle2 = ((2 * G4UniformRand() - 1) * (CLHEP::pi / 6));
        G4double uPx = unit.x() * std::cos(angle1) - unit.y() * std::sin(angle1);
        G4double uPy = unit.x() * std::sin(angle1) + unit.y() * std::cos(angle1);
        G4double uPPy = uPy * std::cos(angle2) - unit.z() * std::sin(angle2);
        G4double uPz = uPy * std::sin(angle2) + unit.z() * std::cos(angle2);
        fParticleChange.ProposeMomentumDirection(
                G4ThreeVector(uPx, uPPy, uPz).unit());
    }
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    State(fMomentumChanged) = true;
    fParticleChange.SetMomentumChanged(true);

#ifdef DEBUG_MEM
    mem_intermediaire = MemoryUsage();
  mem_diff = mem_intermediaire - mem_first;
  G4cout << "\t\t\t >> || MEM || In DrTransportation_Subdiffusion_CTRW::"
            "After proposing new direction to fParticleChange for "
         << track.GetTrackID() << ", diff is : " << mem_diff << G4endl;
#endif

    return; // &fParticleChange is the final returned object
}

// NOT USED
G4double DrTransportation_Subdiffusion_CTRW::ComputeGeomLimit(const G4Track &track, G4double &presafety, G4double limit) {
    G4double res = DBL_MAX;
    if (track.GetVolume() != fpSafetyHelper->GetWorldVolume()) {
        G4TrackStateManager &trackStateMan =
                GetIT(track)->GetTrackingInfo()->GetTrackStateManager();
        fpSafetyHelper->LoadTrackState(trackStateMan);
        res = fpSafetyHelper->CheckNextStep(
                track.GetStep()->GetPreStepPoint()->GetPosition(),
                track.GetMomentumDirection(), limit, presafety);
        fpSafetyHelper->ResetTrackState();
    }
    return res;
}

G4double DrTransportation_Subdiffusion_CTRW::AlongStepGetPhysicalInteractionLength(const G4Track &track,
                                                                                    G4double previousStepSize,
                                                                                    G4double currentMinimumStep,
                                                                                    G4double &currentSafety,
                                                                                    G4GPILSelection *selection) {
#ifdef G4VERBOSE
    if (fVerboseLevel) {
        G4cout << " DrTransportation_Subdiffusion_CTRW::"
                  "AlongStepGetPhysicalInteractionLength - track ID: "
               << track.GetTrackID() << G4endl;
        G4cout << "In volume : " << track.GetVolume()->GetName()
               << " position : " << G4BestUnit(track.GetPosition(), "Length")
               << G4endl;
    }
#endif

    G4double geometryStepLength = G4ITTransportation::AlongStepGetPhysicalInteractionLength(track,
            previousStepSize,
            currentMinimumStep,
            currentSafety,
            selection);

    if (geometryStepLength == 0) {
        if (State(fGeometryLimitedStep)) {
            G4TouchableHandle newTouchable = new G4TouchableHistory;

            newTouchable->UpdateYourself(State(fCurrentTouchableHandle)->GetVolume(), State(fCurrentTouchableHandle)->GetHistory());

            fLinearNavigator->SetGeometricallyLimitedStep();
            fLinearNavigator->LocateGlobalPointAndUpdateTouchableHandle( track.GetPosition(), track.GetMomentumDirection(), newTouchable, true);

            if (newTouchable->GetVolume() == 0) return 0;

            State(fCurrentTouchableHandle) = newTouchable;

            //=======================================
            // Longer but safer ...
            geometryStepLength = G4ITTransportation::AlongStepGetPhysicalInteractionLength(track,
                    previousStepSize,
                    currentMinimumStep,
                    currentSafety,
                    selection);
        }
    }

    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    G4bool isWaiting = DrBreakTable::Instance()->GetBreakMolecule(track, "Sub")->fIsWaiting;
    G4double diffusionCoefficient = (isWaiting)? fTrappedDiffCoef : fJumpDiffCoef;
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    State(fComputeLastPosition) = false;
    State(fTimeStepReachedLimit) = false;
    if (State(fGeometryLimitedStep)) {
        // 95 % of the space step distribution is lower than
        // d_95 = 2 * sqrt(2*D*t)
        // where t is the corresponding time step
        // so by inversion :
        if (fUseMaximumTimeBeforeReachingBoundary) {
            if (fSpeedMeUp) {
                State(theInteractionTimeLeft) = (geometryStepLength * geometryStepLength) / (diffusionCoefficient); // d_50 - use straight line
            }
            else {
                State(theInteractionTimeLeft) = (currentSafety * currentSafety) / (diffusionCoefficient); // d_50 - use safety
                //=====================================================================
                // State(theInteractionTimeLeft) = (currentSafety * currentSafety)
                //          / (8 * diffusionCoefficient); // d_95
                //=====================================================================
            }
            State(fComputeLastPosition) = true;
        }
        else {
            // Will use a random time - this is precise but long to compute in certain
            // circumstances (many particles - small volumes)
            State(fRandomNumber) = G4UniformRand();
            State(theInteractionTimeLeft) = 1 / (4 * diffusionCoefficient) * pow(geometryStepLength / InvErfc(State(fRandomNumber)), 2);
            State(fTransportEndPosition) = geometryStepLength * track.GetMomentumDirection() + track.GetPosition();
        }

        if (fUseSchedulerMinTimeSteps) {
            double minTimeStepAllowed = G4VScheduler::Instance()->GetLimitingTimeStep();

            if (State(theInteractionTimeLeft) < minTimeStepAllowed) {
                State(theInteractionTimeLeft) = minTimeStepAllowed;
                State(fTimeStepReachedLimit) = true;
                State(fComputeLastPosition) = true;
            }
        } else if (State(theInteractionTimeLeft) < fInternalMinTimeStep) {
            State(fTimeStepReachedLimit) = true;
            State(theInteractionTimeLeft) = fInternalMinTimeStep;
            if (fUseMaximumTimeBeforeReachingBoundary) {
                State(fComputeLastPosition) = true;
            }
        }

        State(fCandidateEndGlobalTime) = track.GetGlobalTime() + State(theInteractionTimeLeft);
        State(fEndGlobalTimeComputed) = true; // MK: ADDED ON 20/11/2014
        State(fPathLengthWasCorrected) = false;
    }
    else {
        // Transform geometrical step
        geometryStepLength = 2 * sqrt(diffusionCoefficient * State(theInteractionTimeLeft)) * InvErf(G4UniformRand());
        State(fPathLengthWasCorrected) = true;
        State(fTransportEndPosition) = geometryStepLength * track.GetMomentumDirection() + track.GetPosition();
    }

#ifdef G4VERBOSE
    if (fVerboseLevel > 1) {
        G4cout << GREEN_ON_BLUE
               << "DrTransportation_Subdiffusion_CTRW::"
                  "AlongStepGetPhysicalInteractionLength = "
               << G4BestUnit(geometryStepLength, "Length")
               << " | trackID = " << track.GetTrackID() << RESET_COLOR << G4endl;
    }
#endif

    return geometryStepLength;
}

//////////////////////////////////////////////////////////////////////////
//
//   Initialize ParticleChange  (by setting all its members equal
//                               to corresponding members in G4Track)
G4VParticleChange * DrTransportation_Subdiffusion_CTRW::AlongStepDoIt(const G4Track &track, const G4Step &step) {
#ifdef DEBUG_MEM
    MemStat mem_first, mem_second, mem_diff;
#endif

#ifdef DEBUG_MEM
    mem_first = MemoryUsage();
#endif

    if (GetIT(track)->GetTrackingInfo()->IsLeadingStep() && State(fGeometryLimitedStep)) {

        double spaceStep = DBL_MAX;

        if (State(theInteractionTimeLeft) <= fInternalMinTimeStep) {
            spaceStep = State(fEndPointDistance);
            State(fGeometryLimitedStep) = true;
        }
        else {
            //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            G4bool isWaiting = DrBreakTable::Instance()->GetBreakMolecule(track, "Sub")->fIsWaiting;
            G4double diffusionCoefficient = (isWaiting) ? fTrappedDiffCoef : fJumpDiffCoef;
            //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

            double sqrt_2Dt = sqrt(2 * diffusionCoefficient * State(theInteractionTimeLeft));
            double x = G4RandGauss::shoot(0, sqrt_2Dt);
            double y = G4RandGauss::shoot(0, sqrt_2Dt);
            double z = G4RandGauss::shoot(0, sqrt_2Dt);
            spaceStep = sqrt(x * x + y * y + z * z);

            if (spaceStep >= State(fEndPointDistance)) {
                State(fGeometryLimitedStep) = true;
                if (fUseSchedulerMinTimeSteps == false && spaceStep >= State(fEndPointDistance)) {
                    spaceStep = State(fEndPointDistance);
                }
            }
            else {
                State(fGeometryLimitedStep) = false;
            }
        }
        // Calculate final position
        //
        State(fTransportEndPosition) = track.GetPosition() + spaceStep * track.GetMomentumDirection();
    }

    if (fVerboseLevel) {
        G4cout << GREEN_ON_BLUE
               << "DrTransportation_Subdiffusion_CTRW::AlongStepDoIt: "
                  "GeometryLimitedStep = "
               << State(fGeometryLimitedStep) << RESET_COLOR << G4endl;
    }

    G4ITTransportation::AlongStepDoIt(track, step);

#ifdef DEBUG_MEM
    MemStat mem_intermediaire = MemoryUsage();
  mem_diff = mem_intermediaire - mem_first;
  G4cout << "\t\t\t >> || MEM || After calling G4ITTransportation::"
            "AlongStepDoIt for "
         << track.GetTrackID() << ", diff is : " << mem_diff << G4endl;
#endif

    if (track.GetStepLength() != 0) {
        Diffusion(track);
    }

#ifdef DEBUG_MEM
    mem_intermediaire = MemoryUsage();
  mem_diff = mem_intermediaire - mem_first;
  G4cout << "\t\t\t >> || MEM || After calling "
            "DrTransportation_Subdiffusion_CTRW::"
            "Diffusion for "
         << track.GetTrackID() << ", diff is : " << mem_diff << G4endl;
#endif

    return &fParticleChange;
}
