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

#ifndef TsEmDNAChemistryExtended_hh
#define TsEmDNAChemistryExtended_hh 1

#include "G4VPhysicsConstructor.hh"
#include "G4VUserChemistryList.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "globals.hh"
#include <map>

class TsParameterManager;

class G4DNAMolecularReactionTable;

class TsEmDNAChemistryExtended
: public G4VUserChemistryList,
  public G4VPhysicsConstructor
{
    
public:
    
    explicit TsEmDNAChemistryExtended();
    TsEmDNAChemistryExtended(TsParameterManager* pM);
    
    virtual ~TsEmDNAChemistryExtended();
    
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
    std::map<G4String, G4double> fRadius;
    std::vector< std::vector<G4String> > fReactionSpecies;
    std::vector< std::vector<G4String> > fReactionProducts;
    std::vector<G4double> fReactionRates;
    std::vector<G4String> fReactionTypes;

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


    std::vector<G4String> fScavengedMolecules;
    std::vector<G4double> fScavengedCapacities;
    std::vector<G4bool> fScavengerHasProducts;
    std::vector<G4String> fRemoveInMaterialTheseMolecules;
    G4String fRemoveInMaterial;
    
    struct TsExtendedMolecularReactionData {
        G4int reactionType;
        G4double kobs;
        G4double probabilityForReaction;
        G4double reactionRadius;
        G4double effectiveReactionRadius;
        G4double effectiveTildeReactionRadius;
        G4double alpha;
    };
    
};

#endif
