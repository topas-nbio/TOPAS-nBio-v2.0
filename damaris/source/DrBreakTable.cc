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
#include "DrBreakTable.hh"
#include <G4Scheduler.hh>
#include "DrDefinitions.hh"
#include <G4Molecule.hh>

using namespace std;

DrBreakTable* DrBreakTable::fBreakTableInstance(nullptr);

DrBreakTable::DrBreakTable(): fNextBreakID(0), fCurrentBiologyRepeatNumber(0), fInitialBreakNumber(0), fBreakMolAuxIndex(0) {
    fBreakMolAuxIndex = G4PhysicsModelCatalog::GetIndex("DrBreakMolecule");
}
DrBreakTable::~DrBreakTable(){}

DrBreakTable* DrBreakTable::Instance(){
    if (!fBreakTableInstance) fBreakTableInstance = new DrBreakTable;
    return fBreakTableInstance;
}

void DrBreakTable::NewBreakID(DrBreakMolecule* breakMolecule){
    breakMolecule->sBreakEndA->fOriginalBreakMoleculeID = fNextBreakID;
    AddToTable(fNextBreakID,breakMolecule);
    fNextBreakID++;
}

DrBreakMolecule* DrBreakTable::LinkNewBreakMolecule(G4Track* track){
    auto* breakMolecule = new DrBreakMolecule();
    track->SetAuxiliaryTrackInformation
            (G4PhysicsModelCatalog::GetIndex("DrBreakMolecule"),breakMolecule);
    return breakMolecule;
}

DrBreakMolecule* DrBreakTable::CopyBreakInfo(DrBreakMolecule* mother){
    auto* daughter = new DrBreakMolecule();

    daughter->operator =(*mother);

    return daughter;
}

void DrBreakTable::LinkThisBreakMolecule (G4Track* track,
                                          DrBreakMolecule* breakMolecule){

    DrBreakMolecule* newBreakMolecule = CopyBreakInfo(breakMolecule);
    RemoveFromTable(breakMolecule->sBreakEndA->fOriginalBreakMoleculeID);
    AddToTable(newBreakMolecule->sBreakEndA->fOriginalBreakMoleculeID,newBreakMolecule);

    G4int newBreakID = newBreakMolecule->sBreakEndA->fOriginalBreakMoleculeID;

    auto itFindBreak = fBreakTable.find(newBreakID);
    if(itFindBreak != fBreakTable.end()){
        track->SetAuxiliaryTrackInformation
                (G4PhysicsModelCatalog::GetIndex("DrBreakMolecule"),newBreakMolecule);
    }
    else{
        G4cerr <<"WARNING: Cannot find break molecule with ID "
               <<newBreakID <<" in break table (LinkThisBreak)"
               <<G4endl;
    }
}


void DrBreakTable::LinkThisBreakMolecule(const G4Track& trackOld,
                                         const G4Track& trackNew){

    DrBreakMolecule* breakMolecule = (DrBreakMolecule*)trackOld.GetAuxiliaryTrackInformation(fBreakMolAuxIndex);
    trackOld.SetAuxiliaryTrackInformation(fBreakMolAuxIndex,nullptr);
    trackNew.SetAuxiliaryTrackInformation(fBreakMolAuxIndex,breakMolecule);
}

void DrBreakTable::RemoveFromTable(G4int breakID){
    auto itFindBreak = fBreakTable.find(breakID);
    if (itFindBreak != fBreakTable.end()) fBreakTable.erase(breakID);
    else{
        G4cerr <<"WARNING: Track has no associated DSB molecule (Remove)"
               <<G4endl;
    }
}

void DrBreakTable::AddToTable(G4int breakID,DrBreakMolecule* breakMolecule){
    DrBreakMolecule* storageCopyOfBreakMol = CopyBreakInfo(breakMolecule);
    fBreakTable [breakID] = storageCopyOfBreakMol;
}

void DrBreakTable::JoinBreakMolecule
        (const G4Track& trackReactant1, const G4Track& trackReactant2, G4Track* trackProduct){
    DrBreakMolecule* break1 = GetBreakMolecule(trackReactant1, "JoinBreak");
    DrBreakMolecule* break2 = GetBreakMolecule(trackReactant2, "JoinBreak");
    G4int break1ID = break1->sBreakEndA->fOriginalBreakMoleculeID;
    G4int break2ID = break2->sBreakEndA->fOriginalBreakMoleculeID;

    auto itFindBreak1 = fBreakTable.find(break1ID);
    auto itFindBreak2 = fBreakTable.find(break2ID);
    if(itFindBreak1 != fBreakTable.end() && itFindBreak2 != fBreakTable.end()){

        //@@@@ copying BreakEndA info from Break2 into BreakEndB of Break1
        //@@@@ Break2 will be destroyed after this and Break1 (which now
        //@@@@ contains all the data from both BreakEnds) will be associated
        //@@@@ to the product of the reaction.
        break1->sBreakEndB = break2->sBreakEndA;

        if(break1->sBreakEndA->fCorrectPartnerBreakMoleculeID != break1->sBreakEndB->fOriginalBreakMoleculeID){
            DrBreakTable::Instance()->fMisrepairTimeStore.push_back(trackProduct->GetGlobalTime());
        }

        break1->sBreakEndB->fPKcsIsBleached = break2->sBreakEndA->fPKcsIsBleached;

        LinkThisBreakMolecule(trackProduct, break1);
        RemoveFromTable(break2ID);
    }else{
        G4cerr <<"WARNING: One or both track(s) has no associated DSB molecule (Join)"
               <<G4endl;
        exit(EXIT_FAILURE);
    }
}

void DrBreakTable::SplitBreakMolecule
        (const G4Track& trackMother, G4Track* trackDaughter1, G4Track* trackDaughter2){
    DrBreakMolecule* mother = GetBreakMolecule(trackMother, "SplitBreak");
    G4int motherID = mother->sBreakEndA->fOriginalBreakMoleculeID;

    auto itFindMother = fBreakTable.find(motherID);
    if(itFindMother != fBreakTable.end()){
        DrBreakMolecule* daughter1 = LinkNewBreakMolecule(trackDaughter1);
        DrBreakMolecule* daughter2 = LinkNewBreakMolecule(trackDaughter2);

        daughter1->sBreakEndA = mother->sBreakEndA;
        daughter1->sBreakEndA->fPKcsIsBleached = false;

        daughter2->sBreakEndA = mother->sBreakEndB;
        daughter2->sBreakEndA->fPKcsIsBleached = false;

        RemoveFromTable(motherID);
        AddToTable(daughter1->sBreakEndA->fOriginalBreakMoleculeID,daughter1);
        AddToTable(daughter2->sBreakEndA->fOriginalBreakMoleculeID,daughter2);
    }else{
        G4cerr <<"WARNING: Mother track has no associated DSB molecule (Split)"
               <<G4endl;
    }
}

DrBreakMolecule* DrBreakTable::GetBreakMolecule(const G4Track& track){
    DrBreakMolecule* breakMolecule = (DrBreakMolecule*)
            (track.GetAuxiliaryTrackInformation(fBreakMolAuxIndex));
    G4int breakID = breakMolecule->sBreakEndA->fOriginalBreakMoleculeID;

    auto itFindBreak = fBreakTable.find(breakID);
    if (itFindBreak != fBreakTable.end()) return breakMolecule;
    else{
        G4cerr <<"WARNING: Break " <<breakID
               <<" cannot be found in the break table (Get)"
               <<G4endl;
        return 0;
    }
}

DrBreakMolecule* DrBreakTable::GetBreakMolecule(const G4Track& track, G4String procName){
    DrBreakMolecule* breakMolecule = (DrBreakMolecule*)
            (track.GetAuxiliaryTrackInformation(fBreakMolAuxIndex));
    G4int breakID = breakMolecule->sBreakEndA->fOriginalBreakMoleculeID;

    auto itFindBreak = fBreakTable.find(breakID);
    if (itFindBreak != fBreakTable.end()) return breakMolecule;
    else{
        G4cerr <<"WARNING (Get from " <<procName <<"): Break " <<breakID
               <<" associated with track" <<track.GetTrackID()
               <<" cannot be found in the break table" <<G4endl
               <<"The molecule type is "
               <<GetMolecule(track)->GetDefinition()->GetName() <<G4endl
               <<"With attached AuxInfo (#"
               <<track.GetAuxiliaryTrackInformationMap()->size()<<")"
               <<G4endl;
        return 0;
    }
}

void DrBreakTable::UpdateTypeTracking(std::map<const G4MoleculeDefinition*,G4int> updateList, G4double time){

    if(fMoleculesRecord.size() < 2 ) fMoleculesRecord.resize(2);
    else if(fMoleculesRecord.size() > 2){
        G4cout << "ERROR: MoleculesRecord storage not of expected size"
                << G4endl
                << "Expected 0,1, or 2 but got: "
                << fMoleculesRecord.size()
                << G4endl;
        exit(EXIT_FAILURE);
    }
    for(auto read: updateList){
        // [vector index][TIME][MOL-DEF][NTS, err, repeat, NTS, err]
        fMoleculesRecord[1][time][read.first].resize(5);
        fMoleculesRecord[1][time][read.first] = {(G4double)read.second, 0.,1., (G4double)read.second, 0.};
    }
}

void DrBreakTable::SetUpTypeTrackingTable(){

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
    fBleachedStoreRun[0.0] = {0.0};
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
}

void DrBreakTable::ClearStaticMembers(){

    fBreakTable.clear();
    fInitialBreakNumber = 0;
    fNextBreakID = 0;
}

