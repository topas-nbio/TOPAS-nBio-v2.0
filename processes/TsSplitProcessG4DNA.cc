// Extra Class for TsPhysicsListNBIO
// ********************************************************************
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * TOPAS collaboration.                                             *
// * Use or redistribution of this code is not permitted without the  *
// * explicit approval of the TOPAS collaboration.                    *
// * Contact: Joseph Perl, perl@slac.stanford.edu                     *
// *                                                                  *
// ********************************************************************
//

#include "TsSplitProcessG4DNA.hh"
#include "TsTrackInformation.hh"

#include "G4Track.hh"
#include "G4VParticleChange.hh"
#include "G4LogicalVolume.hh"
#include "G4TouchableHandle.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"

#include <vector>

TsSplitProcessG4DNA::TsSplitProcessG4DNA(G4String regName, G4int nsplit)
:fRegionName(regName), fNSplit(nsplit)
{
    fRegion = G4RegionStore::GetInstance()->FindOrCreateRegion(fRegionName);
}


TsSplitProcessG4DNA::~TsSplitProcessG4DNA() {
}


G4VParticleChange* TsSplitProcessG4DNA::PostStepDoIt(const G4Track& track, const G4Step& step) {
    G4VParticleChange* particleChange(0);
    
    if ( fRegion != step.GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetRegion() ) {
        particleChange = pRegProcess->PostStepDoIt(track, step);
        assert ( 0 != particleChange);
        return particleChange;
    }
    
    if (  1 == fNSplit ) {
        particleChange = pRegProcess->PostStepDoIt(track, step);
        assert ( 0 != particleChange);
        return particleChange;
    }
    
    TsTrackInformation* parentInformation = (TsTrackInformation*)(step.GetTrack()->GetUserInformation());
    G4int initialSplitTrackID = parentInformation->GetSplitTrackID();
    
    if ( 1 < initialSplitTrackID ) {
        particleChange = pRegProcess->PostStepDoIt(track, step);
        assert ( 0 != particleChange);
        return particleChange;
    }
    
    G4double weight = track.GetWeight()/fNSplit;
    G4int splitTrackID = 3;
    
    std::vector<G4Track*> secondaries;
    std::vector<G4int> vSplitTrack;
    
    for ( int i = 0; i < fNSplit; i++ ) {
        particleChange = pRegProcess->PostStepDoIt(track, step);
        assert( 0 !=  particleChange);
        particleChange->SetVerboseLevel(0);
        G4Track* newTrack = new G4Track(*(particleChange->GetSecondary(0)));
        
        secondaries.push_back( newTrack );
        vSplitTrack.push_back( splitTrackID );
        
        splitTrackID++;
    }
    
    parentInformation->SetSplitTrackID(2);
    
    particleChange->SetNumberOfSecondaries(secondaries.size());
    particleChange->SetSecondaryWeightByProcess(true);
    
    std::vector<G4Track*>::iterator iter = secondaries.begin();
    G4int i = 0;
    while( iter != secondaries.end() ) {
        G4Track* newTrack = *iter;
        newTrack->SetWeight(weight);
        
        TsTrackInformation* secondaryInformation = new TsTrackInformation();
        secondaryInformation->SetSplitTrackID(vSplitTrack[i]);
        newTrack->SetUserInformation(secondaryInformation);
        
        particleChange->AddSecondary(newTrack);
        
        iter++;
        i++;
    }
    
    return particleChange;
}
