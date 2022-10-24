//
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

#ifndef TsSplitProcessG4DNA_hh
#define TsSplitProcessG4DNA_hh 1

#include "G4WrapperProcess.hh"

class TsSplitProcessG4DNA : public G4WrapperProcess {
    
public:
    TsSplitProcessG4DNA(G4String regName, G4int split);
    
    virtual ~TsSplitProcessG4DNA();
    
    G4VParticleChange* PostStepDoIt(const G4Track& track, const G4Step& step);
    
private:
    G4Region* fRegion;
    G4String fRegionName;
    G4int fNSplit;
};

#endif
