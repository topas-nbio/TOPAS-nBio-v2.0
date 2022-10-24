// Extra Class for TsEmDNAChemistry

#include "TsDNARemoveInMaterial.hh"

#include <G4VScheduler.hh>
#include "G4SystemOfUnits.hh"
#include "G4Molecule.hh"
#include "G4DNAMolecularMaterial.hh"
#include "G4MolecularConfiguration.hh"
#include "G4MoleculeFinder.hh"
#include "G4VMoleculeCounter.hh"
#include "G4UnitsTable.hh"
#include "G4TrackingInformation.hh"
#include "G4DNADamage.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"

#ifndef State
#define State(theXInfo) (GetState<RemoveInMaterialState>()->theXInfo)
#endif

void TsDNARemoveInMaterial::Create() {
    pParticleChange = &fParticleChange;
    enableAtRestDoIt    = false;
    enableAlongStepDoIt = false;
    enablePostStepDoIt  = true;
    
    SetProcessSubType(60);
    
    G4VITProcess::SetInstantiateProcessState(false);
    
    fIsInitialized = false;
	fpMolecularConfiguration = 0;
    fpMaterial = 0;
    
    fProposesTimeStep = true;
    fReturnedValue = DBL_MAX;
    
    verboseLevel = 0;
}


TsDNARemoveInMaterial::TsDNARemoveInMaterial(const G4String &aName, G4ProcessType type)
:G4VITProcess(aName,type)
{
    Create();
}


TsDNARemoveInMaterial::TsDNARemoveInMaterial(const TsDNARemoveInMaterial& rhs)
:G4VITProcess(rhs)
{
    Create();
}

TsDNARemoveInMaterial::~TsDNARemoveInMaterial()
{;}


TsDNARemoveInMaterial& TsDNARemoveInMaterial::operator=(const TsDNARemoveInMaterial& rhs)
{
    if (this == &rhs) return *this;
    return *this;
}


TsDNARemoveInMaterial::RemoveInMaterialState::RemoveInMaterialState() : G4ProcessState()
{
    fPreviousTimeAtPreStepPoint = -1;
    fIsInGoodMaterial = false;
}


void TsDNARemoveInMaterial::BuildPhysicsTable(const G4ParticleDefinition&)
{
    fIsInitialized = true;
}


void TsDNARemoveInMaterial::StartTracking(G4Track* track)
{
    G4VProcess::StartTracking(track);
    G4VITProcess::fpState.reset(new RemoveInMaterialState());
    G4VITProcess::StartTracking(track);
}


void TsDNARemoveInMaterial::SetReaction(const G4MolecularConfiguration* molConf,
                                        const G4Material* material) 
{
    if(fIsInitialized) {
        G4ExceptionDescription exceptionDescription ;
        exceptionDescription << "TsDNARemoveInMaterial was already initialised. ";
        exceptionDescription << "You cannot set a reaction after initialisation.";
        G4Exception("TsDNARemoveInMaterial::SetReaction","TsDNARemoveInMaterial001",
                    FatalErrorInArgument,exceptionDescription);
    }
   
	
	fpMolecularConfiguration = molConf;

    fpMaterial = material; 
}


G4double TsDNARemoveInMaterial::PostStepGetPhysicalInteractionLength(const G4Track& track,
                                                                       G4double,
                                                                       G4ForceCondition* pForceCond)
{
    const G4Material* material = track.GetMaterial();
    if (material != fpMaterial) {
        return DBL_MAX;
    }

    G4Molecule* mol = GetMolecule(track);
    if(!mol) return DBL_MAX;
	
	if(mol->GetMolecularConfiguration() != fpMolecularConfiguration)
		return DBL_MAX;
	
    *pForceCond = Forced;
	
    return 0.0 * ps;
}


G4VParticleChange* TsDNARemoveInMaterial::PostStepDoIt(const G4Track& track,const G4Step&)
{
    fParticleChange.Initialize(track);
    fParticleChange.ProposeTrackStatus(fStopAndKill);
    State(fPreviousTimeAtPreStepPoint) = -1;
    return &fParticleChange;
}

