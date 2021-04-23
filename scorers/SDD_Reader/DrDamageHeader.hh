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

#include "DrDamageEvent.hh"

class DrDamageHeader {
public:
    DrDamageHeader();
    ~DrDamageHeader();
    //@@@@ This function will print to screen the contents of the damage event.
    //@@@@ All fields are printed regardless of if they exist or not and are
    //@@@@ seperated by semi-colons.
    void PrintHeader();
    //@@@@ This function will print to the specified file the contents of the
    //@@@@ damage event. All fields are printed regardless of if they exist or
    //@@@@ not and are seperated by semi-colons.
    void PrintToFileHeader(G4String);

    G4String fSDD_version{"unspecified"};
    G4String fSoftware{"unspecified"};
    G4String fAuthor{"unspecified"};
    G4String fSimulation_details{"unspecified"};
    G4String fSource{"unspecified"};
    G4int fSource_type{-1};
    std::vector<G4int> fIncident_particles;
    std::vector<G4double> fMean_particle_energy;
    std::vector<G4double> fParticle_fraction;
    G4double fDNA_density{-1.0};
    G4int fIn_vitro_in_vivo{-1};
    G4double fTime{-1.0};
    G4String fAdditional_Information{"unspecified"};
    std::map<G4String,G4int> fData_entries;

    G4String fEnergy_distribution{"unspecified"};
    G4String fDose_or_fluence{"unspecified"};
    //std::pair<G4int,std::pair<G4double,G4double> > fDose_or_fluence;
    G4String fVolumes{"unspecified"};
    G4String fDamage_definition{"unspecified"};
    G4String fDose_rate{"unspecified"};
    //std::vector<G4double> fDose_rate;
    G4String fIrradiation_target{"unspecified"};
    G4String fChromosome_sizes{"unspecified"};
    //std::pair<G4int,std::vector<G4double> > fChromosome_sizes;
    G4String fCell_cycle_phase{"unspecified"};
    //std::pair<G4int,G4double> fCell_cycle_phase;
    G4String fDNA_structure{"unspecified"};
    //std::pair<G4int,G4int> fDNA_structure;
    G4String fProliferation_status{"unspecified"};
    //std::pair<G4int,G4String> fProliferation_status;
    G4String fMicroenvironment{"unspecified"};
    //std::pair<G4double,G4double> fMicroenvironment;
    G4String fDamage_and_primary_count{"unspecified"};
    //std::pair<G4int,G4int> fDamage_and_primary_count;

    //@@@@ This is here in case the user wants to concretely associate the list
    //@@@@ of damage events to the header and access it only though this
    //@@@@ structure.
    std::vector<DrDamageEvent*> fpListOfDamages;
};

