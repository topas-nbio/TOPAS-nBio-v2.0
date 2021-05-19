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

#include <G4VAuxiliaryTrackInformation.hh>
#include <vector>
#include <G4ThreeVector.hh>

struct breakInformation{
    //@@@@ Incorporating FRAP
    G4bool fKuIsBleached{false};
    G4bool fPKcsIsBleached{false};
    G4bool fXRCC4IsBleached{false};

    //@@@@ Protein Attachment
    G4int fKuAttached{0};
    G4int fPKcsAttached{0};
    G4int fXRCC4Attached{0};

    //@@@@ More positional data
    G4bool fTowardsCentromere;

    //@@@@ Double strand break variables
    // store original IDs of DSB in synapsis
    G4int fOriginalBreakMoleculeID = -1;
    // ID of correct partner DSB
    G4int fCorrectPartnerBreakMoleculeID = -1;
    // used to record and asses the starting points of DSB molecules
    G4ThreeVector fOriginalPosition;

   //@@@@ ([DNA Structure];[Chromosome #];[Chromatid #];[Chromosome Arm])
    //@@@@ ---- DNA Structure:
    //@@@@ 0) unspecified, 1) heterochromatin, 2) euchromatin 3) free DNA fragment
    //@@@@ 4) mitocondrial/bacteria/viral DNA
    //@@@@ ---- Chromosome Arm:
    //@@@@ 0) short, 1) long
    std::vector<G4int> fChromasomeID{-1,-1,-1,-1};
    G4double fChromasomePosition{-1.0};
    //@@@@ ([Damage Type],[Direct #],[Indirect #])
    //@@@@ ---- Damage type
    //@@@@ 0) Direct, 1) Indirect, 2) Combination of 1&2, 3) Charge Migration
    std::vector<G4int> fCause{-1,-1,-1};
    //@@@@ ([Base Damage #],[Single Strand Break #],[Double Strand Break #])
    std::vector<G4int> fDamageTypes{-1,-1,-1};
    std::vector<std::vector<G4int> > fFullBreakStructure;
    //@@@@ 0) Missing, 1) A, 2) C, 3) T, 4) G
    std::vector<std::vector<G4int> > fDNASequence;
    //@@@@ Single Number Only: Time for whole damage simulate
    //@@@@ Otherwise: time for each individual lesion in fFullBreakStructure order
    std::vector<G4double> fLesionTime;
};

class DrBreakMolecule : public G4VAuxiliaryTrackInformation
{
public:
    DrBreakMolecule();
    ~DrBreakMolecule();
    DrBreakMolecule& operator=(const DrBreakMolecule& other);
    void PrintBreakDetails();
    void PrintPartBreakDetails(breakInformation* info, std::ofstream&);

    breakInformation* sBreakEndA;
    breakInformation* sBreakEndB;

    //@@@@ Sub-diffusion Variables
    // used to switch between waiting and diffusing states
    G4bool fIsWaiting;
    // waiting time used for anomolous diffusion
    G4double fWaitingTime;
    // diffusion time used for anomolous diffusion
    G4double fDiffusionTime;

    //@@@@ Protein Attachment
    // Are some proteins only attached to the whole site?

    //@@@@ Might be interesting to report the final
    //@@@@ ligated structure?
    std::vector<std::vector<G4int> > fFullBreakStructure;
    // 0) Missing, 1) A, 2) C, 3) T, 4) G
    std::vector<std::vector<G4int> > fDNASequence;
};