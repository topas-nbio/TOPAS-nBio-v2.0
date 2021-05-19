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
// Created by John-William Warmenhoven on 29/10/18.
// DaMaRiS is developed at the University of Manchester.
// See README for references.
//

#include "DrProteinKinetics_Generic.hh"
#include "DrPrecompiler.hh"
#include "DrBreakTable.hh"
#include <G4Molecule.hh>
#include <G4SystemOfUnits.hh>
#include <G4Scheduler.hh>
#include <G4UnitsTable.hh>

#ifndef State
#define State(theXInfo) (GetState<PKGenericState>()->theXInfo)
#endif

static G4int sProcessSubType{99};

DrProteinKinetics_Generic::DrProteinKinetics_Generic(G4String procName, G4MoleculeDefinition* fromMol,
        std::vector<G4MoleculeDefinition*> toMol, G4double procTime, G4bool checkAddLesions,G4int cleanAdditional):
G4VITDiscreteProcess(procName, fUserDefined), fVerbose(0){
    // meaning this class contains a class inheriting from G4ProcessState
    G4VITProcess::SetInstantiateProcessState(false);
    pParticleChange = &aParticleChange;
    enableAlongStepDoIt = false;
    enableAtRestDoIt = false;
    enablePostStepDoIt = true;
    SetProcessSubType(sProcessSubType);
    sProcessSubType--;
    fProposesTimeStep = true;
    fFromMolecule = fromMol;
    fToMolecule = toMol;
    fProcessTime = procTime;
    fNumberSecondaries = toMol.size();
    fProcessName = procName;
    fCheckAddLesions = checkAddLesions;
    fCleanAddLesions = cleanAdditional;

    if(fVerbose > 0){
        G4cout << "Contructing process " << procName << " for paticle "
               << fromMol->GetName() << " channging to ";
        for (auto read: toMol){
            G4cout << read->GetName() <<", ";
        }
        G4cout << "with time constant " << G4BestUnit(procTime,"Time") << G4endl;
    }
}

DrProteinKinetics_Generic::~DrProteinKinetics_Generic(){}

DrProteinKinetics_Generic::DrProteinKinetics_Generic(const DrProteinKinetics_Generic & right):
G4VITDiscreteProcess(fProcessName, fUserDefined){
    fVerbose = right.fVerbose;
}

DrProteinKinetics_Generic::PKGenericState::PKGenericState(): G4ProcessState(){
    fPreviousTimeAtPreStepPoint = -1;
}

void DrProteinKinetics_Generic::StartTracking(G4Track* track){
    G4VProcess::StartTracking(track);
    G4VITProcess::fpState.reset(new PKGenericState());
    G4VITProcess::StartTracking(track);
}

G4bool DrProteinKinetics_Generic::IsApplicable(const G4ParticleDefinition& aParticleType){
    const G4String& name = aParticleType.GetParticleName();
    if(fFromMolecule->GetName() == name) return true;
    else return false;
}

G4double DrProteinKinetics_Generic::GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*){
    return fProcessTime;
}

G4double DrProteinKinetics_Generic::PostStepGetPhysicalInteractionLength(const G4Track& track, G4double, G4ForceCondition* condition){
    *condition = NotForced;
    //@@@@ This is to get the length of the previous time step
    G4double previousTimeStep(-1.);
    if(State(fPreviousTimeAtPreStepPoint) != -1){
        previousTimeStep = track.GetGlobalTime()-State(fPreviousTimeAtPreStepPoint);
    }
    State(fPreviousTimeAtPreStepPoint) = track.GetGlobalTime();
    if((previousTimeStep < 0.0) || (fpState->theNumberOfInteractionLengthLeft <= 0.0)){
        ResetNumberOfInteractionLengthLeft();
#ifdef DEBUG_DAMARIS
        //@@@@ WARNING: This is really slow!
        G4String name = fFromMolecule->GetName()+"_to_"+fToMolecule[0]->GetName();
        G4double value = fpState->theNumberOfInteractionLengthLeft * GetMeanFreePath(track,previousTimeStep,condition);
        if(DrBreakTable::Instance()->debugProcMap.find(name) != DrBreakTable::Instance()->debugProcMap.end()){
            std::vector<G4double> tempVec = DrBreakTable::Instance()->debugProcMap[name];
            tempVec.push_back(value/s);
            DrBreakTable::Instance()->debugProcMap[name] = tempVec;
        }
        else{
            std::vector<G4double> tempVec = {value/s};
            DrBreakTable::Instance()->debugProcMap[name] = tempVec;
        }
#endif /*DEBUG_DAMARIS*/

    }
    else if(previousTimeStep > 0.0){
        SubtractNumberOfInteractionLengthLeft(previousTimeStep);
        if(fpState->theNumberOfInteractionLengthLeft <= perMillion){
            fpState->theNumberOfInteractionLengthLeft=0.0;
        }
    }
    else if(fVerbose > 2){
        G4cout<<"DrProteinKinetics_Generic returning 0 time step"<<G4endl;
    }
    fpState->currentInteractionLength = GetMeanFreePath(track,previousTimeStep,condition);
    G4double value = fpState->theNumberOfInteractionLengthLeft * fpState->currentInteractionLength;

    if(value < 0){
        G4cerr <<"ERROR: Post step interaction length returned by"
               <<"DrProteinKinetics_Generic cannot be negative!"
               <<G4endl;
        assert(value >= 0);
    }
    //@@@@ negative lets the stepper know we are returning a time
    return value*-1;
}

G4VParticleChange* DrProteinKinetics_Generic::PKGenericAction(const G4Track& track){
    //@@@@ Safety to make sure everything is cleared and fresh
    aParticleChange.Clear();
    aParticleChange.Initialize(track);

    //--------------------Setting Up Involved Molecules-------------------------
    DrBreakMolecule* motherBreak = DrBreakTable::Instance()->GetBreakMolecule(track, "Generic");
    //--------------------------------------------------------------------------

    //--------------------Cleaning Additional Lesions---------------------------
    if(fCleanAddLesions){
        if(fCleanAddLesions == 1){
            G4int lesion1 = motherBreak->sBreakEndA->fDamageTypes[1];
            G4int lesion2 = motherBreak->sBreakEndB->fDamageTypes[1];
            G4int total = lesion1+lesion2;
            if(total != 0) {
                //Random probability between 0 and 1 to check which action to take
                G4double castDie = G4UniformRand();
                if (castDie <= (G4double) lesion1 / (G4double) total) motherBreak->sBreakEndA->fDamageTypes[1]--;
                else motherBreak->sBreakEndB->fDamageTypes[1]--;
            }
        }
        else if(fCleanAddLesions == 2){
            G4int lesion1 = motherBreak->sBreakEndA->fDamageTypes[0];
            G4int lesion2 = motherBreak->sBreakEndB->fDamageTypes[0];
            G4int total = lesion1+lesion2;
            if(total != 0){
                //Random probability between 0 and 1 to check which action to take
                G4double castDie = G4UniformRand();
                if(castDie <= (G4double)lesion1/(G4double) total) motherBreak->sBreakEndA->fDamageTypes[0]--;
                else motherBreak->sBreakEndB->fDamageTypes[0]--;
            }
        }
    }
    //--------------------------------------------------------------------------

    //-----------------------Checking for Complexities--------------------------
    //@@@@ Checks all base and backbone breaks are fixed before allowing the
    //@@@@ molecule to progress to the final fixed state
    if(fCheckAddLesions){
        G4int Back1 = motherBreak->sBreakEndA->fDamageTypes[1];
        G4int Back2 = motherBreak->sBreakEndB->fDamageTypes[1];
        G4int Base1  = motherBreak->sBreakEndA->fDamageTypes[0];
        G4int Base2  = motherBreak->sBreakEndB->fDamageTypes[0];

        if(Back1 != 0 || Back2 != 0 || Base1 != 0 || Base2 != 0){
            //@@@@ If the molecule still has backbone or base breaks associated
            //@@@@ with it, this will kick the process back without allowing it
            //@@@@ to change the molecule.
            aParticleChange.ProposeTrackStatus(fAlive);
            return &aParticleChange;
        }
    }
    //--------------------------------------------------------------------------

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

    //-------------------------Propose Changes----------------------------------
    //@@@@ Create the daughter
    aParticleChange.SetNumberOfSecondaries(fNumberSecondaries);

    std::vector<G4Molecule*> theDaughterMolecules;
    theDaughterMolecules.reserve(fNumberSecondaries);
    std::vector<G4Track*> theDaughterTracks;
    theDaughterTracks.reserve(fNumberSecondaries);
    for(auto daughterDef: fToMolecule){
        theDaughterMolecules.emplace_back(new G4Molecule(daughterDef));
    }
    for(auto daughterMol: theDaughterMolecules){
        theDaughterTracks.emplace_back(daughterMol->BuildTrack(time,track.GetPosition()));
    }
    for(auto daughterTrack: theDaughterTracks){
        daughterTrack->SetTrackStatus(fAlive);
        aParticleChange.G4VParticleChange::AddSecondary(daughterTrack);
    }

    //@@@@ Update Molecules Present In Simulation
    if(fNumberSecondaries == 1){
        //@@@@ Carries over the break information, the molecule type is changed but
        //@@@@ the break is still the same.
        G4int auxIndex = DrBreakTable::Instance()->fBreakMolAuxIndex;
        auto daughterBreak = (DrBreakMolecule*)track.GetAuxiliaryTrackInformation(auxIndex);
        track.SetAuxiliaryTrackInformation(auxIndex,nullptr);
        theDaughterTracks[0]->SetAuxiliaryTrackInformation(auxIndex,daughterBreak);
    }
    else if(fNumberSecondaries == 2){
        //@@@@ Carries over the break information, the synaptic complex is split and
        //@@@@ each break end keeps it's information.
        DrBreakTable::Instance()->SplitBreakMolecule(track, theDaughterTracks[0], theDaughterTracks[1]);
    }
    //@@@@ Remove mother from the simulation
    aParticleChange.ProposeTrackStatus(fStopAndKill);
    //--------------------------------------------------------------------------
    return &aParticleChange;
}

G4VParticleChange* DrProteinKinetics_Generic::PostStepDoIt(const G4Track& track, const G4Step&){
    if ((track.GetTrackStatus() == fStopButAlive ) || (track.GetTrackStatus() == fStopAndKill )){
        aParticleChange.Initialize(track);
        return &aParticleChange;
    }
    else{
        //@@@@ Clean the number of interaction length
        //@@@@ left when the process triggers
        ClearNumberOfInteractionLengthLeft();
        return PKGenericAction(track);
    }
}
