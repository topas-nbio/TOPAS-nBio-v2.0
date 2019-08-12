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

#ifndef TsCharltonDNA_hh
#define TsCharltonDNA_hh

#include "TsVGeometryComponent.hh"

class TsCharltonDNA : public TsVGeometryComponent
{
public:
    TsCharltonDNA(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
                  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
    ~TsCharltonDNA();
    
    G4VPhysicalVolume* Construct();
    
    void ResolveParameters();
    
private:
    
    G4double BoxHLX, BoxHLY, BoxHLZ;
    
    G4int fNumberOfBasePairs;
    G4double BasePairLength;
    G4double z0;
    
};

#endif
