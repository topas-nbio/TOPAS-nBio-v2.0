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

#pragma once

#include "DrBreakMolecule.hh"
#include <G4MoleculeDefinition.hh>
#include <G4Track.hh>


class DrDSBMoleculeManager {
public:
    DrDSBMoleculeManager();
    ~DrDSBMoleculeManager();

    void NewBreakID(DrBreakMolecule*);
    DrBreakMolecule* LinkNewBreakMolecule(G4Track*);
    void LinkThisBreakMolecule(G4Track*, DrBreakMolecule*);
    void LinkThisBreakMolecule(const G4Track&, const G4Track&);
    void JoinBreakMolecule(const G4Track&, const G4Track&, G4Track*);
    void SplitBreakMolecule(const G4Track&, G4Track*, G4Track*);
    DrBreakMolecule* CopyBreakInfo(DrBreakMolecule*);

    //@@@@@@@@@@@@@@@@@@@
    // TypeTracking
    //
    void SetUpTypeTrackingTable();
    void UpdateTypeTracking(std::map<const G4MoleculeDefinition*,G4int>,G4double);

private:
    static G4int fNextBreakID;
};