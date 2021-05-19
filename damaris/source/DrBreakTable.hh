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
#pragma once

#include "DrBreakMolecule.hh"
#include "DrBreakIterator.hh"
#include "DrPrecompiler.hh"
#include <G4Track.hh>
#include <G4MoleculeDefinition.hh>

typedef DrBreakIterator<DrBreakMolecule> DrBreakMoleculeIterator;

class G4String;

struct checkBreakStore{
    G4int fNumberBreakMoleculeLeft{-1};
    G4int fNumberUnrepaired{-1};
    G4int fMisrepairCounter{-1};
    G4int fInterChromosomeAberration{-1};
    G4int fInterChromatidAberration{-1};
    G4int fInterArmAberration{-1};
    G4int fIntraArmAberration{-1};
};

class DrBreakTable {
    friend class DrUtils;
//Functions
public:
    ~DrBreakTable();
    static DrBreakTable* Instance();
    void NewBreakID(DrBreakMolecule*);
    DrBreakMolecule* LinkNewBreakMolecule(G4Track*);
    void LinkThisBreakMolecule(G4Track*, DrBreakMolecule*);
    void LinkThisBreakMolecule(const G4Track&, const G4Track&);
    void JoinBreakMolecule(const G4Track&, const G4Track&, G4Track*);
    void SplitBreakMolecule(const G4Track&, G4Track*, G4Track*);
    void RemoveFromTable(G4int);
    void AddToTable(G4int,DrBreakMolecule*);
    DrBreakMolecule* GetBreakMolecule(const G4Track&);
    DrBreakMolecule* GetBreakMolecule(const G4Track&, G4String);
    DrBreakMolecule* CopyBreakInfo(DrBreakMolecule*);
    DrBreakMoleculeIterator GetBreakIterator() { return DrBreakMoleculeIterator(this->fBreakTable); }
    void SetUpTypeTrackingTable();
    void UpdateTypeTracking(std::map<const G4MoleculeDefinition*,G4int>,G4double);
    void ClearStaticMembers();
    void PrintSegment(G4String,G4MoleculeDefinition*,std::ofstream&);

protected:
    DrBreakTable();

// Parameters
public:
    static DrBreakTable* fBreakTableInstance;
    typedef std::map<G4int, DrBreakMolecule*> breakTableTypeDef;
    std::vector<G4double> fDisplacementTrackingStore;
    std::map<G4double,std::vector<G4double> > fMSDTrackingStore;
    std::vector<G4double> fResidualEndDisplacementStore;
    std::vector<G4double> fResidualSynDisplacementStore;
    std::map<G4double,std::vector<G4double> > fBleachedStoreRun;
    std::vector<G4double> fMisrepairSeparationsStore;
    std::vector<G4double> fMisrepairTimeStore;
    std::vector<G4int> fMisrepairNumberStore;
    std::map<G4int,checkBreakStore> fCheckBreakStore;
    breakTableTypeDef fBreakTable;
    std::vector< std::map <G4double, std::map <const G4MoleculeDefinition*, std::vector <G4double> > > > fMoleculesRecord;

    G4int fCurrentBiologyRepeatNumber;
    G4int fInitialBreakNumber;
    G4int fBreakMolAuxIndex;
private:
    G4int fNextBreakID;

#ifdef DEBUG_DAMARIS
    public:
    std::map<G4String, std::vector<G4double> > debugProcMap;
    std::vector<G4double> spaceStepStore;
    std::vector<G4double> actualReactionRangeStore;
    std::vector<G4double> suggestedReactionRangeStore;
    std::vector<G4double> diffTimeStore;
#endif /*DEBUG_DAMARIS*/
};
