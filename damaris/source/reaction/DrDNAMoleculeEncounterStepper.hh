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

#include <G4DNAMoleculeEncounterStepper.hh>

class G4VDNAReactionModel;
class G4DNAMolecularReactionTable;
class G4MolecularConfiguration;

class G4Molecule;

/**
 * Given a molecule G4DNAMoleculeEncounterStepper will calculate for its possible reactants
 * what will be the minimum encounter time and the associated molecules.*
 *
 * This model includes dynamical time steps as explained in
 * "Computer-Aided Stochastic Modeling of the Radiolysis of Liquid Water",
 * V. Michalik, M. Begusová, E. A. Bigildeev,
 * Radiation Research, Vol. 149, No. 3 (Mar., 1998), pp. 224-236
 *
 */

class DrDNAMoleculeEncounterStepper : public G4VITTimeStepComputer
{
public:
    DrDNAMoleculeEncounterStepper();
    ~DrDNAMoleculeEncounterStepper() override;
    DrDNAMoleculeEncounterStepper(const DrDNAMoleculeEncounterStepper&) = delete;
    DrDNAMoleculeEncounterStepper& operator=(const DrDNAMoleculeEncounterStepper&) = delete;

    void Prepare() override ;
    G4double CalculateStep(const G4Track&, const G4double&) override;
    G4double CalculateMinTimeStep(G4double, G4double) override;

    void SetReactionModel(G4VDNAReactionModel*);
    G4VDNAReactionModel* GetReactionModel();

    void SetVerbose(int);
    // Final time returned when reaction is available in the reaction table = 1
    // All details = 2

private:
    void InitializeForNewTrack();

    class Utils;
    void CheckAndRecordResults(const Utils&,
#ifdef G4VERBOSE
                               const G4double reactionRange,
#endif
                               G4KDTreeResultHandle&);

  G4bool fHasAlreadyReachedNullTime;
  const G4DNAMolecularReactionTable*& fMolecularReactionTable;
  G4VDNAReactionModel* fReactionModel;
    G4ITTrackHolder* fpTrackContainer;
    G4ITReactionSet* fReactionSet;
  G4int fVerbose;


    class Utils
    {
    public:
        Utils(const G4Track& tA, const G4MolecularConfiguration* mB);
        ~Utils() = default;

        G4double GetConstant() const { return fConstant; }

        const G4Track& fpTrackA;
        const G4MolecularConfiguration* fpMoleculeB;
        const G4Molecule* fpMoleculeA;
        G4double fDiffusionCoefficientMoleculeA;
        G4double fDiffusionCoefficientMoleculeB;
        G4double fConstant;
    };
};
