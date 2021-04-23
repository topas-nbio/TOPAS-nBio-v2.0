// Extra Class for use by TOPAS
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

#include "DrDamageHeader.hh"
#include <fstream>

DrDamageHeader::DrDamageHeader() {}

DrDamageHeader::~DrDamageHeader() {}

void DrDamageHeader::PrintHeader() {
    G4cout<<"SDD version: "<<fSDD_version<<G4endl
          <<"Software: "<<fSoftware<<G4endl
          <<"Author: "<<fAuthor<<G4endl
          <<"Simulation Details: "<<fSimulation_details<<G4endl
          <<"Source: "<<fSource<<G4endl
          <<"Source Type: "<<fSource_type<<G4endl
          <<"Incident Particles:";
    for(auto read: fIncident_particles){
        G4cout<<" "<<read;
    }
    G4cout<<G4endl;
    G4cout<<"Mean Particle Energies:";
    for(auto read: fMean_particle_energy){
        G4cout<<" "<<read;
    }
    G4cout<<G4endl;
    G4cout<<"Energy Distrobutions: "<<fEnergy_distribution<<G4endl;
    G4cout<<"Particle Fractions:";
    for(auto read: fParticle_fraction){
        G4cout<<" "<<read;
    }
    G4cout<<G4endl
          <<"Dose/Fluence: "<<fDose_or_fluence<<G4endl
          <<"Dose Rate: "<<fDose_rate<<G4endl
          <<"Irradiation Target: "<<fIrradiation_target<<G4endl
          <<"Volumes: "<<fVolumes<<G4endl
          <<"Chromasome Sizes: "<<fChromosome_sizes<<G4endl
          <<"DNA Density: "<<fDNA_density<<G4endl
          <<"Cell Cycle Phase: "<<fCell_cycle_phase<<G4endl
          <<"DNA Structure: "<<fDNA_structure<<G4endl
          <<"In vitro/vivo: "<<fIn_vitro_in_vivo<<G4endl
          <<"Proliferation State: "<<fProliferation_status<<G4endl
          <<"Microenvironment: "<<fMicroenvironment<<G4endl
          <<"Damage Definition: "<<fDamage_definition<<G4endl
          <<"Time: "<<fTime<<G4endl
          <<"Damage and Primary Count: "<<fDamage_and_primary_count<<G4endl
          <<"Additional Information: "<<fAdditional_Information<<G4endl<<G4endl
          <<"Data Entries:";
    for(auto read: fData_entries){
        G4cout<<" "<<read.second;
    }
    G4cout<<G4endl<<G4endl;
}

void DrDamageHeader::PrintToFileHeader(G4String fileName) {
    std::ofstream file;
    file.open(fileName, std::ios_base::app);

    file<<"SDD version:"<<fSDD_version<<";"<<G4endl
        <<"Software:"<<fSoftware<<";"<<G4endl
        <<"Author:"<<fAuthor<<";"<<G4endl
        <<"Simulation Details:"<<fSimulation_details<<";"<<G4endl
        <<"Source:"<<fSource<<";"<<G4endl
        <<"Source type:"<<fSource_type<<";"<<G4endl
        <<"Incident particles:";
    G4bool multi{false};
    for(auto read: fIncident_particles){
        if(multi)file<<",";
        file<<read;
        multi=true;
    }
    file<<";"<<G4endl;
    file<<"Mean particle energy:";
    multi=false;
    for(auto read: fMean_particle_energy){
        if(multi)file<<",";
        file<<read;
        multi=true;
    }
    file<<";"<<G4endl;
    file<<"Energy Distribution:"<<fEnergy_distribution<<";"<<G4endl;
    file<<"Particle fraction:";
    multi=false;
    for(auto read: fParticle_fraction){
        if(multi)file<<",";
        file<<read;
        multi=true;
    }
    file<<";"<<G4endl
        <<"Dose or fluence:"<<fDose_or_fluence<<";"<<G4endl
        <<"Dose Rate:"<<fDose_rate<<";"<<G4endl
        <<"Irradiation target:"<<fIrradiation_target<<";"<<G4endl
        <<"Volumes:"<<fVolumes<<";"<<G4endl
        <<"Chromosome sizes:"<<fChromosome_sizes<<";"<<G4endl
        <<"DNA Density:"<<fDNA_density<<";"<<G4endl
        <<"Cell Cycle Phase:"<<fCell_cycle_phase<<";"<<G4endl
        <<"DNA Structure:"<<fDNA_structure<<";"<<G4endl
        <<"In vitro / in vivo:"<<fIn_vitro_in_vivo<<";"<<G4endl
        <<"Proliferation status:"<<fProliferation_status<<";"<<G4endl
        <<"Microenvironment:"<<fMicroenvironment<<";"<<G4endl
        <<"Damage definition:"<<fDamage_definition<<";"<<G4endl
        <<"Time:"<<fTime<<";"<<G4endl
        <<"Damage and primary count:"<<fDamage_and_primary_count<<";"<<G4endl
        <<"Data entries:";
    multi=false;
    for(auto read: fData_entries){
        if(multi) file<<",";
        file<<read.second;
        multi=true;
    }
    file<<";"<<G4endl<<"Additional information:"<<fAdditional_Information<<";"<<G4endl
        <<"***EndOfHeader***;"<<G4endl;
}
