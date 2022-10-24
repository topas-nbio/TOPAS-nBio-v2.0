// Extra Class for TsEmDNAChemistry

#include "TsDNAFirstOrderReaction.hh"

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
#define State(theXInfo) (GetState<FirstOrderReactionState>()->theXInfo)
#endif

void TsDNAFirstOrderReaction::Create() {
    pParticleChange = &fParticleChange;
    enableAtRestDoIt    = false;
    enableAlongStepDoIt = false;
    enablePostStepDoIt  = true;
    
    SetProcessSubType(60);
    
    G4VITProcess::SetInstantiateProcessState(false);
    
    fIsInitialized = false;
    fpMolecularConfiguration = 0;
    
    fProposesTimeStep = true;
    fReturnedValue = DBL_MAX;
    
    verboseLevel = 0;
}


TsDNAFirstOrderReaction::TsDNAFirstOrderReaction(const G4String &aName, G4ProcessType type)
:G4VITProcess(aName,type)
{
    Create();
}


TsDNAFirstOrderReaction::TsDNAFirstOrderReaction(const TsDNAFirstOrderReaction& rhs)
:G4VITProcess(rhs)
{
    Create();
}

TsDNAFirstOrderReaction::~TsDNAFirstOrderReaction()
{;}


TsDNAFirstOrderReaction& TsDNAFirstOrderReaction::operator=(const TsDNAFirstOrderReaction& rhs)
{
    if (this == &rhs) return *this;
    return *this;
}


TsDNAFirstOrderReaction::FirstOrderReactionState::FirstOrderReactionState() : G4ProcessState()
{
    fPreviousTimeAtPreStepPoint = -1;
    fIsInGoodMaterial = false;
}


void TsDNAFirstOrderReaction::BuildPhysicsTable(const G4ParticleDefinition&)
{
    fIsInitialized = true;
}


void TsDNAFirstOrderReaction::StartTracking(G4Track* track)
{
    G4VProcess::StartTracking(track);
    G4VITProcess::fpState.reset(new FirstOrderReactionState());
    G4VITProcess::StartTracking(track);
}


void TsDNAFirstOrderReaction::SetReaction(const G4MolecularConfiguration* molConf,
                                          double scavengingCapacity)
{
    if(fIsInitialized) {
        G4ExceptionDescription exceptionDescription ;
        exceptionDescription << "TsDNAFirstOrderReaction was already initialised. ";
        exceptionDescription << "You cannot set a reaction after initialisation.";
        G4Exception("TsDNAFirstOrderReaction::SetReaction","TsDNAFirstOrderReaction001",
                    FatalErrorInArgument,exceptionDescription);
    }
    
    fpMolecularConfiguration = molConf;
    fScavengingCapacity = scavengingCapacity;
    fHasProducts = false;
}


void TsDNAFirstOrderReaction::SetReaction(const G4MolecularConfiguration* molConf,
                                          std::vector<G4MolecularConfiguration*> products,
                                          double scavengingCapacity)
{
    if(fIsInitialized) {
        G4ExceptionDescription exceptionDescription ;
        exceptionDescription << "TsDNAFirstOrderReaction was already initialised. ";
        exceptionDescription << "You cannot set a reaction after initialisation.";
        G4Exception("TsDNAFirstOrderReaction::SetReaction","TsDNAFirstOrderReaction001",
                    FatalErrorInArgument,exceptionDescription);
    }
    
    fpMolecularConfiguration = molConf;
    fScavengingCapacity = scavengingCapacity;
    for ( size_t u = 0; u < products.size(); u++ ) {
        fpProductsMolecularConfiguration.push_back(products[u]);
    }
    fHasProducts = true;
}


G4double TsDNAFirstOrderReaction::PostStepGetPhysicalInteractionLength(const G4Track& track,
                                                                       G4double,
                                                                       G4ForceCondition* pForceCond)
{
    G4Molecule* mol = GetMolecule(track);
    if(!mol) return DBL_MAX;
    
    if(mol->GetMolecularConfiguration() != fpMolecularConfiguration) {
        return DBL_MAX;
    }
    
    State(fIsInGoodMaterial) = true;
    
    G4double previousTimeStep(-1.);
    
    if(State(fPreviousTimeAtPreStepPoint) != -1) {
        previousTimeStep = track.GetGlobalTime() -
        State(fPreviousTimeAtPreStepPoint) ;
    }
    
    State(fPreviousTimeAtPreStepPoint) = track.GetGlobalTime();
    
    // condition is set to "Not Forced"
    *pForceCond = NotForced;
    
    if ((previousTimeStep < 0.0) || (fpState->theNumberOfInteractionLengthLeft<=0.0)) {
        ResetNumberOfInteractionLengthLeft();
    } else if ( previousTimeStep  > 0.0) {
        SubtractNumberOfInteractionLengthLeft( previousTimeStep );
    }
    
    // input variable from parameter manager
    fpState->currentInteractionLength = 1./fScavengingCapacity;
    
    G4double value;
    if (fpState->currentInteractionLength <DBL_MAX) {
        value = fpState->theNumberOfInteractionLengthLeft * (fpState->currentInteractionLength);
    } else {
        value = DBL_MAX;
    }
    
    if(value < fReturnedValue)
        fReturnedValue  = value;
    
    return value*-1;
}


G4VParticleChange* TsDNAFirstOrderReaction::PostStepDoIt(const G4Track& track,const G4Step&)
{
    fReturnedValue  = DBL_MAX;
    fParticleChange.Initialize(track);
    if ( fHasProducts ) {
        for ( size_t u = 0; u < fpProductsMolecularConfiguration.size(); u++ ) {
            G4Molecule* product = new G4Molecule(fpProductsMolecularConfiguration[u]);
            G4DNADamage::Instance()->AddIndirectDamage("NULL",product,track.GetPosition(),track.GetGlobalTime());
	    
            //G4Track* productTrack = product->BuildTrack(track.GetGlobalTime(),track.GetPosition());
	    
            //productTrack->SetTrackStatus(fAlive);
           
	    //G4MoleculeFinder::Instance()->Push(productTrack);

            //std::cout << " Created product from scavenging " << std::endl;
            //if(G4VMoleculeCounter::InUse())
            //    G4VMoleculeCounter::Instance()->AddAMoleculeAtTime(
            //                                                   GetMolecule(productTrack)->GetMolecularConfiguration(),
            //                                                  track.GetGlobalTime(),&(track.GetPosition()));
        }
    }
    fParticleChange.ProposeTrackStatus(fStopAndKill);
    State(fPreviousTimeAtPreStepPoint) = -1;
    return &fParticleChange;
}

