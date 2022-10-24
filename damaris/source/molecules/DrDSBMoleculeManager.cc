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
// Created by john-william on 24/06/21.
//

#include "DrDSBMoleculeManager.hh"
#include "DrDefinitions.hh"
#include "TsParameterManager.hh"
#include <G4Scheduler.hh>
#include <G4Molecule.hh>

G4int DrDSBMoleculeManager::fNextBreakID = 0;

DrDSBMoleculeManager::DrDSBMoleculeManager() {}
DrDSBMoleculeManager::~DrDSBMoleculeManager() {}

void DrDSBMoleculeManager::NewBreakID(DrBreakMolecule* breakMolecule){
    breakMolecule->sBreakEndA->fOriginalBreakMoleculeID = fNextBreakID;
    fNextBreakID++;
}

DrBreakMolecule* DrDSBMoleculeManager::LinkNewBreakMolecule(G4Track* track){
    auto* breakMolecule = new DrBreakMolecule();
    track->SetAuxiliaryTrackInformation(G4PhysicsModelCatalog::GetIndex("DrBreakMolecule"),breakMolecule);
    return breakMolecule;
}

DrBreakMolecule* DrDSBMoleculeManager::CopyBreakInfo(DrBreakMolecule* mother){
    auto* daughter = new DrBreakMolecule();
    daughter->operator =(*mother);
    return daughter;
}

void DrDSBMoleculeManager::LinkThisBreakMolecule (G4Track* track, DrBreakMolecule* breakMolecule){
    DrBreakMolecule* newBreakMolecule = CopyBreakInfo(breakMolecule);
    track->SetAuxiliaryTrackInformation(G4PhysicsModelCatalog::GetIndex("DrBreakMolecule"),newBreakMolecule);
}


void DrDSBMoleculeManager::LinkThisBreakMolecule(const G4Track& trackOld, const G4Track& trackNew){
    DrBreakMolecule* breakMolecule = (DrBreakMolecule*)trackOld.GetAuxiliaryTrackInformation(G4PhysicsModelCatalog::GetIndex("DrBreakMolecule"));
    trackOld.SetAuxiliaryTrackInformation(G4PhysicsModelCatalog::GetIndex("DrBreakMolecule"),nullptr);
    trackNew.SetAuxiliaryTrackInformation(G4PhysicsModelCatalog::GetIndex("DrBreakMolecule"),breakMolecule);
}

void DrDSBMoleculeManager::JoinBreakMolecule(const G4Track& trackReactant1, const G4Track& trackReactant2, G4Track* trackProduct){

    DrBreakMolecule* break1 = (DrBreakMolecule*)(trackReactant1.GetAuxiliaryTrackInformation(G4PhysicsModelCatalog::GetIndex("DrBreakMolecule")));
    DrBreakMolecule* break2 = (DrBreakMolecule*)(trackReactant2.GetAuxiliaryTrackInformation(G4PhysicsModelCatalog::GetIndex("DrBreakMolecule")));

    //@@@@ copying BreakEndA info from Break2 into BreakEndB of Break1
    //@@@@ Break2 will be destroyed after this and Break1 (which now
    //@@@@ contains all the data from both BreakEnds) will be associated
    //@@@@ to the product of the reaction.
    break1->sBreakEndB = break2->sBreakEndA;

    if(break1->sBreakEndA->fCorrectPartnerBreakMoleculeID != break1->sBreakEndB->fOriginalBreakMoleculeID){
        DrDefinitions::Instance()->fMisrepairTimeStore.push_back(trackProduct->GetGlobalTime());
    }

    break1->sBreakEndB->fPKcsIsBleached = break2->sBreakEndA->fPKcsIsBleached;

    LinkThisBreakMolecule(trackProduct, break1);
}

void DrDSBMoleculeManager::SplitBreakMolecule(const G4Track& trackMother, G4Track* trackDaughter1, G4Track* trackDaughter2){

    DrBreakMolecule* mother = (DrBreakMolecule*)(trackMother.GetAuxiliaryTrackInformation(G4PhysicsModelCatalog::GetIndex("DrBreakMolecule")));
    DrBreakMolecule* daughter1 = LinkNewBreakMolecule(trackDaughter1);
    DrBreakMolecule* daughter2 = LinkNewBreakMolecule(trackDaughter2);

    daughter1->sBreakEndA = mother->sBreakEndA;
    daughter1->sBreakEndA->fPKcsIsBleached = false;

    daughter2->sBreakEndA = mother->sBreakEndB;
    daughter2->sBreakEndA->fPKcsIsBleached = false;
}

void DrDSBMoleculeManager::UpdateTypeTracking(std::map<const G4MoleculeDefinition*,G4int> updateList, G4double time){

    if(DrDefinitions::Instance()->fMoleculesRecord.size() < 2 ) DrDefinitions::Instance()->fMoleculesRecord.resize(2);
    else if(DrDefinitions::Instance()->fMoleculesRecord.size() > 2){
        G4cout << "ERROR: MoleculesRecord storage not of expected size"
               << G4endl
               << "Expected 0,1, or 2 but got: "
               << DrDefinitions::Instance()->fMoleculesRecord.size()
               << G4endl;
        DrDefinitions::Instance()->GetParameterManager()->AbortSession(1);
    }
    for(auto read: updateList){
        // [vector index][TIME][MOL-DEF][NTS, err, repeat, NTS, err]
        DrDefinitions::Instance()->fMoleculesRecord[1][time][read.first].resize(5);
        DrDefinitions::Instance()->fMoleculesRecord[1][time][read.first] = {(G4double)read.second, 0.,1., (G4double)read.second, 0.};
    }
}

void DrDSBMoleculeManager::SetUpTypeTrackingTable(){

    std::map<const G4MoleculeDefinition*, G4int > moleculeCount; // Molecule Name: time: count)

    for(auto read: DrDefinitions::Instance()->GetNameMap()){
        moleculeCount[read.second] = 0;
    }

    G4TrackManyList *mainList = G4ITTrackHolder::Instance()->GetMainList();
    G4TrackManyList::iterator it = mainList->begin();
    G4TrackManyList::iterator end = mainList->end();

    for (; it != end; ++it) {
        G4Track *track = *it;
        const G4MoleculeDefinition *molDef = GetMolecule(track)->GetDefinition();

        if (moleculeCount.find(molDef) == moleculeCount.end()) {
            // not found molecule, it is not one defined in pathway.in
        }
        else {
            // found molecule, add one to count
            moleculeCount[molDef]++;
        }
    }
    //----------------------------------------------------------------------
    // Add the recorded moleculeCount above into TypeTable
    //----------------------------------------------------------------------
    UpdateTypeTracking(moleculeCount,0);

    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    //set up bleached molecules too
    DrDefinitions::Instance()->fBleachedStoreRun[0.0] = {0.0};
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
}
