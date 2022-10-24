//

#ifndef TsDNAFirstOrderReaction_HH
#define TsDNAFirstOrderReaction_HH

#include "G4VITProcess.hh"

class G4MolecularConfiguration;

class TsDNAFirstOrderReaction : public G4VITProcess
{
public:
    
    TsDNAFirstOrderReaction(const G4String& aName =  "DNAFirstOrderReaction",
                            G4ProcessType type = fDecay);
    
    virtual ~TsDNAFirstOrderReaction();
    
    G4IT_ADD_CLONE(G4VITProcess,TsDNAFirstOrderReaction)
    
    TsDNAFirstOrderReaction(const TsDNAFirstOrderReaction&);
    
    TsDNAFirstOrderReaction& operator=(const TsDNAFirstOrderReaction&);
    
    void StartTracking(G4Track*);
    
    void SetReaction(const G4MolecularConfiguration*, double /*scavengingCapacity*/);
    
    void SetReaction(const G4MolecularConfiguration*, std::vector<G4MolecularConfiguration*>, double /*scavengingCapacity*/);

    virtual void BuildPhysicsTable(const G4ParticleDefinition&);
    
    virtual G4double PostStepGetPhysicalInteractionLength(const G4Track& track,
                                                          G4double   previousStepSize,
                                                          G4ForceCondition* condition);
    
    virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);
    
    virtual G4double AtRestGetPhysicalInteractionLength(const G4Track&,
                                                        G4ForceCondition*){ return -1.0; }
    
    virtual G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&){ return 0; }
    
    virtual G4double AlongStepGetPhysicalInteractionLength(const G4Track&,
                                                           G4double  ,
                                                           G4double  ,
                                                           G4double& ,
                                                           G4GPILSelection*){ return -1.0; }
    
    virtual G4VParticleChange* AlongStepDoIt(const G4Track& ,
                                             const G4Step&) {return 0;}
    
protected:
    struct FirstOrderReactionState : public G4ProcessState
    {
        FirstOrderReactionState();
        virtual ~FirstOrderReactionState(){;}
        G4double fPreviousTimeAtPreStepPoint;
        G4bool fIsInGoodMaterial;
    };
    
    G4bool fIsInitialized;
    G4double fReturnedValue;
    const std::vector<double>* fpMoleculeDensity;
    G4double fReactionRate;
    G4double fConcentration;
    G4double fMolarMassOfMaterial;
    G4ParticleChange fParticleChange;
    
    const G4MolecularConfiguration* fpMolecularConfiguration;
    const G4Material* fpMaterial;
    std::vector<G4MolecularConfiguration*> fpProductsMolecularConfiguration;
    G4bool fHasProducts;
private :
    void Create();
    G4double fScavengingCapacity;
};

#endif
