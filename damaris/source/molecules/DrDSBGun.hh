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

#include "DrMoleculeGun.hh"
#include "DrDamageEvent.hh"

struct DoubleStrandBreak {
    G4int fBackboneBreak1{-1};
    G4int fBaseBreak1{-1};
    G4int fBackboneBreak2{-1};
    G4int fBaseBreak2{-1};
};

using breakListTypeDef = std::vector<std::pair<G4ThreeVector, DoubleStrandBreak>>;
using storeTypeDef = std::vector<std::vector<std::vector< DrDamageEvent* > > >;

class DrDSBGun {
public:
    DrDSBGun();
    ~DrDSBGun();

    void DefineAllTracks();
    G4int PlaceBreaks(std::vector<std::vector<DrDamageEvent*> >);
    DoubleStrandBreak ConstructBreak(std::vector<std::vector<G4int>>);
    void PlaceClockMolecule();
    void AddMoleculeToGun(const G4ThreeVector& position, const G4String& name, const G4double& time);
    G4int BuildDSBEnd_Offset(G4String name, G4double offset);
    G4int BuildDSBEnd_Origin(G4int number, G4String name);
    G4int BuildDSB_Column(G4int number, G4double width_nm);
    G4int BuildDSB_SepSpaceAndTime(G4double seperation, G4int backbone, G4int base, G4int DSB_or_DSBEnd, G4double timeDelay);
    G4int BuildDSB_Pattern(G4int backbone, G4int base, std::vector<G4int> dsbPattern);
    G4int PlaceFromSTDFile(const G4String&);
    G4int PlaceFromSTDInput(std::vector<std::vector<DrDamageEvent*> >);
    DoubleStrandBreak ParseSTDInput(DrDamageEvent*);
    G4int SelectRandomFromSTDInput(storeTypeDef);
    void SetBreakMolPropertiesFromEvent(DrDamageEvent*,DrBreakMolecule*);

private:
    DrMoleculeGun *fMoleculeGun;
    G4double fBoundingRadius;
};
