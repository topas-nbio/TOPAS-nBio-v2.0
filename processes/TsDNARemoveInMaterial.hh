//

#ifndef TsDNARemoveInMaterial_HH
#define TsDNARemoveInMaterial_HH

#include "G4VITProcess.hh"

class G4MolecularConfiguration;

class TsDNARemoveInMaterial : public G4VITProcess
{
public:
    
    TsDNARemoveInMaterial(const G4String& aName =  "DNARemoveInMaterial",
                            G4ProcessType type = fDecay);
    
    virtual ~TsDNARemoveInMaterial();
    
    G4IT_ADD_CLONE(G4VITProcess,TsDNARemoveInMaterial)
    
    TsDNARemoveInMaterial(const TsDNARemoveInMaterial&);
    
    TsDNARemoveInMaterial& operator=(const TsDNARemoveInMaterial&);
    
    void StartTracking(G4Track*);
    
    void SetReaction(const G4MolecularConfiguration*, const G4Material*); 

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
    struct RemoveInMaterialState : public G4ProcessState
    {
        RemoveInMaterialState();
        virtual ~RemoveInMaterialState(){;}
        G4double fPreviousTimeAtPreStepPoint;
        G4bool fIsInGoodMaterial;
    };
    
    G4bool fIsInitialized;
    G4double fReturnedValue;
    G4ParticleChange fParticleChange;
    
    const G4Material* fpMaterial;
    const G4MolecularConfiguration* fpMolecularConfiguration;

private :
    void Create();

};

#endif
