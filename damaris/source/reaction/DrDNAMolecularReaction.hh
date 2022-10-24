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

#include <G4VITReactionProcess.hh>

class G4DNAMolecularReactionTable;
class G4VDNAReactionModel;

/**
  * G4DNAMolecularReaction is the reaction process
  * used in G4DNAMolecularStepByStepModel between
  * two molecules.
  * After the global track steps, it test whether
  * the molecules can react. If so, the reaction is made.
  */

class DrDNAMolecularReaction : public G4VITReactionProcess
{
    public:
    DrDNAMolecularReaction();
    explicit DrDNAMolecularReaction(G4VDNAReactionModel*);
    ~DrDNAMolecularReaction() override = default;
    DrDNAMolecularReaction(const DrDNAMolecularReaction& other) = delete;
    DrDNAMolecularReaction& operator=(const DrDNAMolecularReaction& other) = delete;

    G4bool TestReactibility(const G4Track&,
                            const G4Track&,
                            double currentStepTime,
                            bool userStepTimeLimit) override;

    std::vector<std::unique_ptr<G4ITReactionChange>> FindReaction(G4ITReactionSet*, const double, const double, const bool) override;
    std::unique_ptr<G4ITReactionChange> MakeReaction(const G4Track&, const G4Track&) override;

    void SetReactionModel(G4VDNAReactionModel*);

protected:
    const G4DNAMolecularReactionTable*& fMolReactionTable;
    G4VDNAReactionModel* fpReactionModel;
};