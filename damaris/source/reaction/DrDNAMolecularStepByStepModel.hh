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

#include <G4DNAMolecularStepByStepModel.hh>

class G4DNAMolecularReactionTable;

/**
  * G4DNAMolecularStepByStepModel :
  *  - TimeStepper : DrDNAMolecularEncounterStepper
  *  - ReactionProcess : DrDNAMolecularReaction
  * Before each step, the next minimum encounter time is calculated for each
  * pair of molecule. The minimum time step is selected. All the molecules are stepped
  * within this time step. Then, only the relevant pair of molecules are checked for
  * reaction.
  */

class DrDNAMolecularStepByStepModel : public G4DNAMolecularStepByStepModel{
public:
    DrDNAMolecularStepByStepModel(const G4String& name = "DrDNAMolecularStepByStepModel");
   virtual void Initialize();
};