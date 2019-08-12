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

#ifndef TsLinearDNA_hh
#define TsLinearDNA_hh

#include "TsVGeometryComponent.hh"


class TsLinearDNA : public TsVGeometryComponent
{    
public:
	TsLinearDNA(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
				  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	~TsLinearDNA();
	
	G4VPhysicalVolume* Construct();
    
    void ResolveParameters();
    
private:
    G4int fNumberOfBasePairs;
    G4double HL;
  
};

#endif
