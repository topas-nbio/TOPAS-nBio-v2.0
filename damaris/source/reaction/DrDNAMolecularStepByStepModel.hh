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

#include <G4VITStepModel.hh>

class G4DNAMolecularReactionTable;
class G4VDNAReactionModel;


/**
  * G4DNAMolecularStepByStepModel :
  *  - TimeStepper : DrDNAMolecularEncounterStepper
  *  - ReactionProcess : DrDNAMolecularReaction
  * Before each step, the next minimum encounter time is calculated for each
  * pair of molecule. The minimum time step is selected. All the molecules are stepped
  * within this time step. Then, only the relevant pair of molecules are checked for
  * reaction.
  */

class DrDNAMolecularStepByStepModel : public G4VITStepModel{
public:
    DrDNAMolecularStepByStepModel(const G4String& name = "DrDNAMolecularStepByStepModel");
    void Initialize() override ;

    DrDNAMolecularStepByStepModel(const G4String& name,
                                  std::unique_ptr<G4VITTimeStepComputer> pTimeStepper,
                                  std::unique_ptr<G4VITReactionProcess> pReactionProcess);

    DrDNAMolecularStepByStepModel& operator=(const DrDNAMolecularStepByStepModel&) = delete;
    DrDNAMolecularStepByStepModel(const DrDNAMolecularStepByStepModel&) = delete;
    ~DrDNAMolecularStepByStepModel() override;

    void PrintInfo() override;

    void SetReactionModel(G4VDNAReactionModel*);
    G4VDNAReactionModel* GetReactionModel();

protected:
    const G4DNAMolecularReactionTable*& fMolecularReactionTable;
    std::unique_ptr<G4VDNAReactionModel> fpReactionModel;

};