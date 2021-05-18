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
// Created by John-William Warmenhoven on 18/10/18.
// Please cite: https://doi.org/10.1667/RR15209.1
//

#pragma once

#include <vector>
#include <G4ThreeVector.hh>

class DrDamageEvent {
public:
    DrDamageEvent();
    ~DrDamageEvent();
    //@@@@ This function will print to screen the contents of the damage event.
    //@@@@ All fields are printed regardless of if they exist or not and are
    //@@@@ seperated by semi-colons.
    void PrintEvent();
    //@@@@ This function will print to the specified file the contents of the
    //@@@@ damage event. All fields are printed regardless of if they exist or
    //@@@@ not and are seperated by semi-colons.
    void PrintToFileEvent(G4String);

    //@@@@ 0) Same event 1) Start of damages associated with a new particle
    //@@@@ 2) Start of damages associated with a new dose
    std::vector<G4int> fNewEvent;
    //@@@@ Geometric position in nucleus
    std::vector<G4ThreeVector> fPosition;
    //@@@@ ([DNA Structure];[Chromosome #];[Chromatid #];[Chromosome Arm])
    //@@@@ ---- DNA Structure:
    //@@@@ 0) unspecified, 1) heterochromatin, 2) euchromatin 3) free DNA fragment
    //@@@@ 4) mitocondrial/bacteria/viral DNA
    //@@@@ ---- Chromosome Arm:
    //@@@@ 0) short, 1) long
    std::vector<G4int> fChromasomeID;
    G4double fChromasomePosition{-1.0};
    //@@@@ ([Damage Type],[Direct #],[Indirect #])
    //@@@@ ---- Damage type
    //@@@@ 0) Direct, 1) Indirect, 2) Combination of 1&2, 3) Charge Migration
    std::vector<G4int> fCause;
    //@@@@ ([Base Damage #],[Single Strand Break #],[Double Strand Break #])
    std::vector<G4int> fDamageTypes;
    G4String fFull_Break_Spec{""};
    G4String fDNA_Sequence{""};
    std::vector<std::vector<G4int> > fFullBreakStructure;
    //@@@@ 0) Missing, 1) A, 2) C, 3) T, 4) G
    std::vector<std::vector<G4int> > fDNASequence;
    //@@@@ Single Number Only: Time for whole damage simulate
    //@@@@ Otherwise: time for each individual lesion in fFullBreakStructure order
    std::vector<G4double> fLesionTime;

    std::vector<G4int> fParticleTypes;
    std::vector<G4double> fEnergies;
    std::vector<G4ThreeVector> fTranslation;
    std::vector<G4ThreeVector> fDirection;
    std::vector<G4double> fParticleTime;
};
