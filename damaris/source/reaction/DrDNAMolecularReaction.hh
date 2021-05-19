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

#include <G4DNAMolecularReaction.hh>

class G4DNAMolecularReactionTable;
class G4VDNAReactionModel;

/**
  * G4DNAMolecularReaction is the reaction process
  * used in G4DNAMolecularStepByStepModel between
  * two molecules.
  * After the global track steps, it test whether
  * the molecules can react. If so, the reaction is made.
  */

class DrDNAMolecularReaction : public G4DNAMolecularReaction
{
    public:
        std::unique_ptr<G4ITReactionChange> MakeReaction(const G4Track&, const G4Track&) override;
};
