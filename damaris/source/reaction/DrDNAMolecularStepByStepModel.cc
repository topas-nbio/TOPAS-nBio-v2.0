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
#include "DrDNAMolecularStepByStepModel.hh"
#include "DrDNAMolecularReaction.hh"
#include "DrDNAMoleculeEncounterStepper.hh"
#include "DrDNASmoluchowskiReactionModel.hh"
#include <G4DNAMolecularReactionTable.hh>


DrDNAMolecularStepByStepModel::DrDNAMolecularStepByStepModel(const G4String& name)
        : G4DNAMolecularStepByStepModel(name,
                                        std::unique_ptr<DrDNAMoleculeEncounterStepper>(new DrDNAMoleculeEncounterStepper()),
                                        std::unique_ptr<DrDNAMolecularReaction>(new DrDNAMolecularReaction()))
{}

void DrDNAMolecularStepByStepModel::Initialize(){

    if(fpReactionTable == nullptr){
        SetReactionTable(G4DNAMolecularReactionTable::GetReactionTable());
    }

    if(!fpReactionModel) fpReactionModel.reset(new DrDNASmoluchowskiReactionModel());

    fpReactionModel->SetReactionTable((const G4DNAMolecularReactionTable*) fpReactionTable);

    ((DrDNAMolecularReaction*) fpReactionProcess.get())->SetReactionModel(fpReactionModel.get());
    ((DrDNAMoleculeEncounterStepper*) fpTimeStepper.get())->SetReactionModel(fpReactionModel.get());

    G4VITStepModel::Initialize();
}