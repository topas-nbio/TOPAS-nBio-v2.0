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
