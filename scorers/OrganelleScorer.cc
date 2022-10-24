// Scorer for OrganelleScorer
//
// ********************************************************************
// *                                                                  *
// * This file is based on code from the TOPAS-nBio extensions to the *
// * TOPAS Simulation Toolkit.                                        *
// * The TOPAS-nBio extensions are freely available under the license *
// * agreement set forth at: https://topas-nbio.readthedocs.io/       *
// *                                                                  *
// ********************************************************************
// *  Modifications by Marc B. Hahn (2020)                            *
// *  Please report bugs to hahn@physiks.fu-berlin.de                 *
// *  or on https://github.com/BAMresearch/TOPAS-CellModels           *
// ********************************************************************
//
// Tuple scorer for scoring energy deposited in the cell organelles

#include "OrganelleScorer.hh"
#include "G4SystemOfUnits.hh"

OrganelleScorer::OrganelleScorer(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                                 G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
: TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{
    //SetScorer();
    fNtuple->RegisterColumnF(&fEnergyDep, "Energy Deposited", "MeV");
    fNtuple->RegisterColumnS(&fVolName, "Volume Name");

}


OrganelleScorer::~OrganelleScorer() {;}


G4bool OrganelleScorer::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
    ResolveSolid(aStep);
    
    G4double flagEnergyDep      = aStep->GetTotalEnergyDeposit();
        
    //Find volume name
    G4Track* aTrack = aStep->GetTrack();
    G4String volumeName = aTrack->GetVolume()->GetName();
    
    
    //Score events that deposit energy in the cell and it's organelles:
    if ((flagEnergyDep > 0)){
        
        //Get Edep
        fEnergyDep = flagEnergyDep;

        //Get volume Name
        fVolName = volumeName;

        fNtuple->Fill();

        return true;
    }
    return false;   
}
