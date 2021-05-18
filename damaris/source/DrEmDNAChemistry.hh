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

#include "G4VPhysicsConstructor.hh"
#include "G4VUserChemistryList.hh"

class TsParameterManager;

class G4DNAMolecularReactionTable;

class DrEmDNAChemistry
        : public G4VUserChemistryList,
          public G4VPhysicsConstructor
{

public:

    explicit DrEmDNAChemistry();
    DrEmDNAChemistry(TsParameterManager* pM);

    virtual ~DrEmDNAChemistry();

    virtual void ConstructParticle()
    {
      ConstructMolecule();
    }

    virtual void ConstructMolecule();
    virtual void ConstructProcess();

    virtual void ConstructDissociationChannels();
    virtual void ConstructReactionTable(G4DNAMolecularReactionTable* reactionTable);
    virtual void ConstructTimeStepModel(G4DNAMolecularReactionTable* reactionTable);

    G4bool IsWaterMolecularConfigurationActivated();

private:
    void DefineParameters();
    G4bool MoleculeExists(G4String);
    G4String GetFullParmName(G4String);
    void Quit(G4String, G4String);

private:
    TsParameterManager* fPm;
    std::map<G4String, G4String> fExistingMolecules;
    G4String fName;

    std::map<G4String, G4double> fDiffusionCoefficients;
    std::vector< std::vector<G4String> > fReactionSpecies;
    std::vector< std::vector<G4String> > fReactionProducts;
    std::vector<G4double> fReactionRates;

    G4bool fSetWaterConfiguration;

    G4double fIonizationStates;
    G4double fA1B1DissociativeDecay;
    G4double fA1B1Relaxation;
    G4double fB1A1AutoIonization;
    G4double fB1A1DissociativeDecay;
    G4double fB1A1Relaxation;
    G4double fRydDiffAutoIonization;
    G4double fRydDiffRelaxation;
    G4double fDissociativeAttachment;
};
