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
#include "DrBiologyTimeStepAction.hh"
#include <G4Scheduler.hh>

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@JWW_DaMaRiS
#include "DrBreakTable.hh"
#include "DrPrecompiler.hh"
#include <G4Molecule.hh>
#include <G4SystemOfUnits.hh>
#include <TsParameterManager.hh>
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

DrBiologyTimeStepAction::DrBiologyTimeStepAction(TsParameterManager* pM, G4String name)
:G4UserTimeStepAction(), fPm(pM)
{
    G4String parmName = "Ch/" + name + "/";
    if ( fPm->ParameterExists(parmName + "AddTimeStepHighEdge") &&
         fPm->ParameterExists(parmName + "AddTimeStepResolution")) {
        G4int nSteps = fPm->GetVectorLength(parmName + "AddTimeStepHighEdge");
        if ( nSteps != fPm->GetVectorLength(parmName + "AddTimeStepResolution")) {
            G4cerr << "TOPAS is exiting due to an error in parameters definition" << G4endl;
            G4cerr << "Length of double vector parameter: " << parmName + "AddTimeStepHighEdge" <<
                   " does not match with legth of double vector parameter: " <<
                   parmName + "AddTimeStepResolution" << G4endl;
        }
        G4double* highEdges = fPm->GetDoubleVector(parmName + "AddTimeStepHighEdge","Time");
        G4double* resolution = fPm->GetDoubleVector(parmName + "AddTimeStepResolution","Time");
        for ( int i = 0; i < nSteps; i++ ) AddTimeStep(highEdges[i], resolution[i]);
    }
}

DrBiologyTimeStepAction::~DrBiologyTimeStepAction() {}

DrBiologyTimeStepAction &DrBiologyTimeStepAction::
operator=(const DrBiologyTimeStepAction &rhs) {
    if (this == &rhs)
        return *this;
    return *this;
}

void DrBiologyTimeStepAction::UserPreTimeStepAction(){
    G4double globalTime = G4Scheduler::Instance()->GetGlobalTime();
    if(globalTime == 0.0) DrBreakTable::Instance()->SetUpTypeTrackingTable();
}

void DrBiologyTimeStepAction::UserPostTimeStepAction() {

    //@@@@ Get track list for checks below
    DrBreakTable *bTable = DrBreakTable::Instance();
    G4TrackManyList *mainList = G4ITTrackHolder::Instance()->GetMainList();
    G4TrackManyList::iterator it = mainList->begin();
    G4TrackManyList::iterator end = mainList->end();

#ifdef DEBUG_DAMARIS
    std::ofstream file("debug_positions.out",std::ios_base::app);
    it = mainList->begin();
    end = mainList->end();
    for (; it != end; ++it) {
        G4Track *track = *it;
        if(!std::isnan(track->GetPosition().mag())){
            file<<track->GetGlobalTime()<<" "
                <<track->GetPosition().mag()/nm<<" "
                <<track->GetPosition().x()/nm<<" "
                <<track->GetPosition().y()/nm<<" "
                <<track->GetPosition().z()/nm
                <<G4endl;
        }
    }
#endif /*DEBUG_DAMARIS*/

#ifdef PLOT_MOTION
     {
        it = mainList->begin();
        end = mainList->end();
        for (; it != end; ++it) {
            G4Track *track = *it;
            G4String molName = GetMolecule(track)->GetName();

            if( molName.substr(0,5) != "Clock"){

                G4int IDA = bTable->GetBreakMolecule(*track)
                        ->sBreakEndA->fOriginalBreakMoleculeID;
                G4int IDACorr = bTable->GetBreakMolecule(*track)
                        ->sBreakEndA->fCorrectPartnerBreakMoleculeID;
                G4int IDB = bTable->GetBreakMolecule(*track)
                        ->sBreakEndB->fOriginalBreakMoleculeID;

                G4String fileName;

                if(IDB != -1){
                    if(IDB != IDACorr){
                        G4String ID;
                        if(IDA > IDB) ID = std::to_string(IDB)+std::to_string(IDA);
                        else ID = std::to_string(IDA)+std::to_string(IDB);
                        fileName = "DEBUG_MOTION"+ID+"_Mis.out";
                    }
                    else{
                        G4String ID;
                        if(IDA > IDB) ID = std::to_string(IDB)+std::to_string(IDA);
                        else ID = std::to_string(IDA)+std::to_string(IDB);
                        fileName = "DEBUG_MOTION"+ID+"_Corr.out";
                    }
                }
                else {
                    fileName = "DEBUG_MOTION"
                               +std::to_string(IDA)
                               +".out";
                }

                std::ofstream debFile(fileName,std::ios_base::app);
                debFile << G4Scheduler::Instance()->GetGlobalTime()/s << " "
                        << track->GetPosition().x()/nm << " "
                        << track->GetPosition().y()/nm << " "
                        << track->GetPosition().z()/nm << " "
                        << molName << " "
                        << G4endl;
            }
        }
    }
#endif /*PLOT_MOTION*/

    if (G4Scheduler::Instance()->GetGlobalTime() ==
        G4Scheduler::Instance()->GetEndTime()) {

        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        //@@@@ This calculates the displacement of every particle left in
        //@@@@ the simulation at the end and pushes it onto a static list
        //@@@@ contained in DrBreakTable to be binned and plotted at
        //@@@@ the end of the whole run.
        it = mainList->begin();
        end = mainList->end();

        for (; it != end; ++it) {
            G4Track *track = *it; //^^^^^^^
            G4Molecule *molecule = GetMolecule(track);
            const G4MoleculeDefinition *moleculeDefinition =
                    molecule->GetDefinition();

            if (moleculeDefinition->GetName().substr(0, 3) == "DSB") {
                DrBreakMolecule *breakMolecule =
                        bTable->GetBreakMolecule(*track, "UserPTSA");

                G4ThreeVector startPosition1;
                G4ThreeVector startPosition2;
                G4double displacement1;
                G4double displacement2;

                if (moleculeDefinition->GetName().substr(0, 6) == "DSBEnd") {
                    startPosition1 = breakMolecule->sBreakEndA->fOriginalPosition;
                    displacement1 = (startPosition1 - track->GetPosition()).mag();
                    bTable->fDisplacementTrackingStore.push_back(displacement1);
                    bTable->fResidualEndDisplacementStore.push_back(displacement1);
                } else if (moleculeDefinition->GetName().substr(0, 12) == "DSB_Fixed_HR"){
                    startPosition1 = breakMolecule->sBreakEndA->fOriginalPosition;
                    displacement1 = (startPosition1 - track->GetPosition()).mag();
                    bTable->fDisplacementTrackingStore.push_back(displacement1);
                } else if (moleculeDefinition->GetName().substr(0, 6) == "DSBSyn") {
                    startPosition1 = breakMolecule->sBreakEndA->fOriginalPosition;
                    startPosition2 = breakMolecule->sBreakEndB->fOriginalPosition;
                    displacement1 = (startPosition1 - track->GetPosition()).mag();
                    displacement2 = (startPosition2 - track->GetPosition()).mag();
                    bTable->fDisplacementTrackingStore.push_back(displacement1);
                    bTable->fDisplacementTrackingStore.push_back(displacement2);
                    bTable->fResidualSynDisplacementStore.push_back(displacement1);
                    bTable->fResidualSynDisplacementStore.push_back(displacement2);
                } else if (moleculeDefinition->GetName().substr(0, 6) == "DSB_Fi") {
                    startPosition1 = breakMolecule->sBreakEndA->fOriginalPosition;
                    startPosition2 = breakMolecule->sBreakEndB->fOriginalPosition;
                    displacement1 = (startPosition1 - track->GetPosition()).mag();
                    displacement2 = (startPosition2 - track->GetPosition()).mag();
                    bTable->fDisplacementTrackingStore.push_back(displacement1);
                    bTable->fDisplacementTrackingStore.push_back(displacement2);
                }
            }
        }
    }
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
}

void DrBiologyTimeStepAction::UserReactionAction(
        const G4Track &trackA, const G4Track &trackB,
        const std::vector<G4Track *> *productsVector) {

#ifdef DEBUG_DAMARIS
  G4ThreeVector tAPosition = trackA.GetPosition();
  G4ThreeVector tBPosition = trackB.GetPosition();
  G4double seperation = (tAPosition-tBPosition).mag();
  DrBreakTable::Instance()->actualReactionRangeStore.push_back(seperation);
#endif /*DEBUG_DAMARIS*/

    G4Track *secondary = (*productsVector)[0];
    if (GetMolecule(trackA)->GetDefinition()->GetName().substr(0, 3) == "DSB") {
        DrBreakTable::Instance()->JoinBreakMolecule(trackA, trackB, secondary);
    }
}
