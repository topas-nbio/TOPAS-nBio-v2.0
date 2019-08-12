// Extra Class for TsEmDNAChemistry
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

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAWaterExcitationStructure.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4NistManager.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAMolecularMaterial.hh"
#include "G4TransportationManager.hh"
#include "G4ITNavigator.hh"
#include "G4Navigator.hh"

//#define MODEL_VERBOSE


template<typename MODEL>
TsTDNAOneStepThermalizationModelh<MODEL>::
TsTDNAOneStepThermalizationModelh(const G4ParticleDefinition*,
                                const G4String& nam) :
G4VEmModel(nam), fIsInitialised(false)
{
  fVerboseLevel = 0;
  SetLowEnergyLimit(0.);
  G4DNAWaterExcitationStructure exStructure;
  SetHighEnergyLimit(exStructure.ExcitationEnergy(0));
  fParticleChangeForGamma = 0;
  fpWaterDensity = 0;
  fNavigator = 0;
}

//------------------------------------------------------------------------------

template<typename MODEL>
TsTDNAOneStepThermalizationModelh<MODEL>::~TsTDNAOneStepThermalizationModelh()
{
  if(fNavigator)
  {
    delete fNavigator;
  }
}

//------------------------------------------------------------------------------
template<typename MODEL>
void TsTDNAOneStepThermalizationModelh<MODEL>::
Initialise(const G4ParticleDefinition* particleDefinition,
           const G4DataVector&)
{
#ifdef MODEL_VERBOSE
  if(fVerboseLevel)
    G4cout << "Calling TsDNAOneStepThermalizationModel::Initialise()"
           << G4endl;
#endif
  if (particleDefinition->GetParticleName() != "e-")
  {
    G4ExceptionDescription errMsg;
    errMsg << "TsDNAOneStepThermalizationModel can only be applied "
    "to electrons";
    G4Exception("TsDNAOneStepThermalizationModel::CrossSectionPerVolume",
                "TsDNAOneStepThermalizationModel001",
                FatalErrorInArgument,errMsg);
    return;
  }
  
  if(!fIsInitialised)
  {
    fIsInitialised = true;
    fParticleChangeForGamma = GetParticleChangeForGamma();
  }
  
  G4Navigator* navigator =
  G4TransportationManager::GetTransportationManager()->
  GetNavigatorForTracking();
  
  fNavigator = new G4Navigator();
  
  if(navigator){ // add these checks for testing mode
    auto world=navigator->GetWorldVolume();
    if(world){
      fNavigator->SetWorldVolume(world);
      //fNavigator->NewNavigatorState();
    }
  }
  
  fpWaterDensity =
  G4DNAMolecularMaterial::Instance()->
  GetNumMolPerVolTableFor(G4Material::GetMaterial("G4_WATER"));
}

//------------------------------------------------------------------------------
template<typename MODEL>
G4double TsTDNAOneStepThermalizationModelh<MODEL>::
CrossSectionPerVolume(const G4Material* material,
                      const G4ParticleDefinition*,
                      G4double ekin,
                      G4double,
                      G4double)
{
#ifdef MODEL_VERBOSE
  if(fVerboseLevel > 1)
    G4cout << "Calling CrossSectionPerVolume() of TsDNAOneStepThermalizationModel"
    << G4endl;
#endif
  
  if(ekin > HighEnergyLimit()){
    return 0.0;
  }
  
  G4double waterDensity = (*fpWaterDensity)[material->GetIndex()];
  
  if(waterDensity!= 0.0){
    return DBL_MAX;
  }
  return 0.;
}

//------------------------------------------------------------------------------
template<typename MODEL>
double TsTDNAOneStepThermalizationModelh<MODEL>::GetRmean(double k){
  return MODEL::GetRmean(k);
}


//------------------------------------------------------------------------------

template<typename MODEL>
void TsTDNAOneStepThermalizationModelh<MODEL>::
GetPenetration(G4double k, G4ThreeVector& displacement)
{
  return MODEL::GetPenetration(k, displacement);
}

//------------------------------------------------------------------------------
template<typename MODEL>
void TsTDNAOneStepThermalizationModelh<MODEL>::
SampleSecondaries(std::vector<G4DynamicParticle*>*,
                  const G4MaterialCutsCouple*,
                  const G4DynamicParticle* particle,
                  G4double,
                  G4double)
{
#ifdef MODEL_VERBOSE
  if(fVerboseLevel)
    G4cout << "Calling SampleSecondaries() of TsDNAOneStepThermalizationModel"
    << G4endl;
#endif
  
  G4double k = particle->GetKineticEnergy();
  
  if (k <= HighEnergyLimit())
  {
    fParticleChangeForGamma->ProposeTrackStatus(fStopAndKill);
    fParticleChangeForGamma->ProposeLocalEnergyDeposit(k);
    
    if(G4DNAChemistryManager::IsActivated())
    {
      G4ThreeVector displacement(0,0,0);
      GetPenetration(k, displacement);
      
      //______________________________________________________________
      const G4Track * theIncomingTrack =
      fParticleChangeForGamma->GetCurrentTrack();
      G4ThreeVector finalPosition(theIncomingTrack->GetPosition()+displacement);
      
      fNavigator->SetWorldVolume(theIncomingTrack->GetTouchable()->
                                 GetVolume(theIncomingTrack->GetTouchable()->
                                           GetHistoryDepth()));
      
      double displacementMag = displacement.mag();
      double safety = DBL_MAX;
      G4ThreeVector direction = displacement/displacementMag;
      
      double mag_displacement = displacement.mag();
      G4ThreeVector displacement_direction = displacement/mag_displacement;
      
      fNavigator->ResetHierarchyAndLocate(theIncomingTrack->GetPosition(),
                                          direction,
                                          *((G4TouchableHistory*)
                                            theIncomingTrack->GetTouchable()));
      
      fNavigator->ComputeStep(theIncomingTrack->GetPosition(),
                              displacement/displacementMag,
                              displacementMag,
                              safety);
      
      if(safety <= displacementMag)
      {
        finalPosition = theIncomingTrack->GetPosition()
        + (displacement/displacementMag)*safety*0.80;
      }
      
      G4DNAChemistryManager::Instance()->CreateSolvatedElectron(theIncomingTrack,
                                                                &finalPosition);
      
      fParticleChangeForGamma->SetProposedKineticEnergy(25.e-3*eV);
    }
  }
}
