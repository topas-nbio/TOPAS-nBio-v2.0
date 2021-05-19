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

#include "DrDamageEvent.hh"
#include <fstream>

DrDamageEvent::DrDamageEvent(){}
DrDamageEvent::~DrDamageEvent(){}
void DrDamageEvent::PrintEvent(){
    G4cout<<G4endl
          <<"Classification: ";
    for(auto read: fNewEvent){
        G4cout<<read<<" ";
    }
    G4cout<<G4endl;

    G4cout<<"Position: "<<G4endl;
    for(auto read: fPosition){
        G4cout<<read<<G4endl;
    }

    G4cout<<"Chromosome ID: "<<G4endl;
    for(auto read: fChromasomeID){
        G4cout<<read<<" ";
    }
    G4cout<<G4endl;

    G4cout<<"Chromosome Position: "
          <<G4endl
          <<fChromasomePosition
          <<G4endl;

    G4cout<<"Cause: "<<G4endl;
    for(auto read: fCause){
        G4cout<<read<<" ";
    }
    G4cout<<G4endl;

    G4cout<<"Damage Types: "<<G4endl;
    for(auto read: fDamageTypes){
        G4cout<<read<<" ";
    }
    G4cout<<G4endl;

    G4cout<<"Full Break Structure: "<<G4endl;
    for(auto read: fFullBreakStructure){
        G4cout<<read[0]<<" "<<read[1]<<" "<<read[2]<<" / ";
    }
    G4cout<<G4endl;

    G4cout<<"DNA Sequence: "<<G4endl;
    for(auto read: fDNASequence){
        G4cout<<read[0]<<" "<<read[1]<<" "<<read[2]<<" / ";
    }
    G4cout<<G4endl;

    G4cout<<"Lesion Time: "<<G4endl;
    for(auto read: fLesionTime){
        G4cout<<read<<" ";
    }
    G4cout<<G4endl;

    G4cout<<"Particle Types: "<<G4endl;
    for(auto read: fParticleTypes){
        G4cout<<read<<" ";
    }
    G4cout<<G4endl;

    G4cout<<"Energies: "<<G4endl;
    for(auto read: fEnergies){
        G4cout<<read<<" ";
    }
    G4cout<<G4endl;

    G4cout<<"Translation: "<<G4endl;
    for(auto read: fTranslation){
        G4cout<<read<<" ";
    }
    G4cout<<G4endl;

    G4cout<<"Direction: "<<G4endl;
    for(auto read: fDirection){
        G4cout<<read<<" ";
    }
    G4cout<<G4endl;

    G4cout<<"Particle Time: "<<G4endl;
    for(auto read: fParticleTime){
        G4cout<<read<<" ";
    }
    G4cout<<G4endl<<G4endl;
}

void DrDamageEvent::PrintToFileEvent(G4String fileName){

    std::ofstream file;
    file.open(fileName, std::ios::app);
    G4bool multi{false};

    file<<G4endl;
    for(auto read: fNewEvent){
        if(multi) file<<",";
        file<<read;
        multi=true;
    }
    file<<";";

    multi=false;
    for(auto read: fPosition){
        if(multi) file<<"/";
        file<<read.x()<<","<<read.y()<<","<<read.z();
        multi=true;
    }
    file<<";";

    multi=false;
    for(auto read: fChromasomeID){
        if(multi) file<<",";
        file<<read;
        multi=true;
    }
    file<<";";

    file<<fChromasomePosition<<";";

    multi=false;
    for(auto read: fCause){
        if(multi) file<<",";
        file<<read;
        multi=true;
    }
    file<<";";

    multi=false;
    for(auto read: fDamageTypes){
        if(multi) file<<",";
        file<<read;
        multi=true;
    }
    file<<";";

    for(auto read: fFullBreakStructure){
        file<<read[0]<<","<<read[1]<<","<<read[2]<<"/";
    }
    file<<";";

    for(auto read: fDNASequence){
        file<<read[0]<<","<<read[1]<<","<<read[2]<<"/";
    }
    file<<";";

    multi=false;
    for(auto read: fLesionTime){
        if(multi) file<<"/";
        file<<read;
        multi=true;
    }
    file<<";";

    multi=false;
    for(auto read: fParticleTypes){
        if(multi) file<<",";
        file<<read;
        multi=true;
    }
    file<<";";

    multi=false;
    for(auto read: fEnergies){
        if(multi) file<<",";
        file<<read;
        multi=true;
    }
    file<<";";

    multi=false;
    for(auto read: fTranslation){
        if(multi) file<<",";
        file<<read.x()<<"/"<<read.y()<<"/"<<read.z();
        multi=true;
    }
    file<<";";

    multi=false;
    for(auto read: fDirection){
        if(multi) file<<",";
        file<<read.x()<<"/"<<read.y()<<"/"<<read.z();
        multi=true;
    }
    file<<";";

    multi=false;
    for(auto read: fParticleTime){
        if(multi) file<<",";
        file<<read;
        multi=true;
    }
    file<<";"<<G4endl;
}
