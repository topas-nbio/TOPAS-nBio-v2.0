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
#include "DrBreakMolecule.hh"
#include "DrBreakTable.hh"
#include <fstream>

DrBreakMolecule::DrBreakMolecule():G4VAuxiliaryTrackInformation(), fIsWaiting(true), fWaitingTime(-1.0), fDiffusionTime(-1.0){
	sBreakEndA = new breakInformation();
	sBreakEndB = new breakInformation();
}
DrBreakMolecule::~DrBreakMolecule(){}

DrBreakMolecule& DrBreakMolecule::operator =(const DrBreakMolecule& rhs){

    this->sBreakEndA = rhs.sBreakEndA;
    this->sBreakEndB = rhs.sBreakEndB;

	//@@@@ Sub-diffusion Variables
	this->fIsWaiting = rhs.fIsWaiting;
	this->fWaitingTime = rhs.fWaitingTime;
	this->fDiffusionTime = rhs.fDiffusionTime;

	this->fFullBreakStructure = rhs.fFullBreakStructure;
	this->fDNASequence = rhs.fDNASequence;

	return *this;
}

void DrBreakMolecule::PrintBreakDetails() {
    G4String ID = DrBreakTable::Instance()->fCurrentBiologyRepeatNumber;
    std::ofstream file("ReportBreakStructure.out", std::ios_base::app);
    file << "===Start of run " << ID <<"===" << G4endl;
    file << "---Start of BreakEndA---" << G4endl;
    PrintPartBreakDetails(sBreakEndA, file);
    file << "---Start of BreakEndB---" << G4endl;
    PrintPartBreakDetails(sBreakEndB, file);
    file << "===End of run " << ID <<"===" << G4endl;
    file << G4endl << G4endl;
}

void DrBreakMolecule::PrintPartBreakDetails(breakInformation* info, std::ofstream& file) {

    file
            << info->fKuIsBleached << G4endl
            << info->fPKcsIsBleached << G4endl
            << info->fXRCC4IsBleached << G4endl
            << info->fKuAttached << G4endl
            << info->fPKcsAttached << G4endl
            << info->fXRCC4Attached << G4endl
            << info->fTowardsCentromere << G4endl
            << info->fOriginalBreakMoleculeID << G4endl
            << info->fCorrectPartnerBreakMoleculeID << G4endl
            << info->fOriginalPosition.x() << ", "
            << info->fOriginalPosition.y() << ", "
            << info->fOriginalPosition.z() << G4endl
            << info->fChromasomeID[0] << ", "
            << info->fChromasomeID[1] << ", "
            << info->fChromasomeID[2] << ", "
            << info->fChromasomeID[4] << G4endl
            << info->fChromasomePosition << G4endl
            << info->fCause[0] << ", "
            << info->fCause[1] << ", "
            << info->fCause[2] << G4endl
            << info->fDamageTypes[0] << ", "
            << info->fDamageTypes[1] << ", "
            << info->fDamageTypes[2] << G4endl;
    for (const auto& read: info->fFullBreakStructure){
        for (const auto& read2: read){
            file << read2 << ", ";
        }
        file << G4endl;
    }
    for (const auto& read: info->fDNASequence){
        for (const auto& read2: read){
            file << read2 << ", ";
        }
        file << G4endl;
    }
    for (const auto& read: info->fLesionTime){
        file << read << ", ";
    }
    file << G4endl;
}