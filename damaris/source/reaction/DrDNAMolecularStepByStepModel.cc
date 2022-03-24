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
#include <G4Molecule.hh>

DrDNAMolecularStepByStepModel::DrDNAMolecularStepByStepModel(const G4String& name)
        : DrDNAMolecularStepByStepModel(name,
                                        std::unique_ptr<DrDNAMoleculeEncounterStepper>(new DrDNAMoleculeEncounterStepper()),
                                        std::unique_ptr<DrDNAMolecularReaction>(new DrDNAMolecularReaction()))
{}

DrDNAMolecularStepByStepModel::DrDNAMolecularStepByStepModel(const G4String& name,
                                                             std::unique_ptr<G4VITTimeStepComputer> pTimeStepper,
                                                             std::unique_ptr<G4VITReactionProcess> pReactionProcess)
        : G4VITStepModel(std::move(pTimeStepper),
                         std::move(pReactionProcess),
                         name)
        , fMolecularReactionTable(reference_cast<const G4DNAMolecularReactionTable*>(fpReactionTable))
{
    //@@@@
    fpTimeStepper.reset(new DrDNAMoleculeEncounterStepper());
    fpReactionProcess.reset(new DrDNAMolecularReaction());
    //@@@@
    fType1 = G4Molecule::ITType();
    fType2 = G4Molecule::ITType();
}

DrDNAMolecularStepByStepModel::~DrDNAMolecularStepByStepModel() = default;

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

void DrDNAMolecularStepByStepModel::PrintInfo()
{
#ifdef G4VERBOSE
    G4cout << fName << " will be used" << G4endl;
#endif
}

void DrDNAMolecularStepByStepModel::SetReactionModel(G4VDNAReactionModel* pReactionModel)
{
    fpReactionModel.reset(pReactionModel);
}

G4VDNAReactionModel* DrDNAMolecularStepByStepModel::GetReactionModel()
{
    return fpReactionModel.get();
}