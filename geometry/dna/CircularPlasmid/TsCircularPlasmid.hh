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


#ifndef TsCircularPlasmid_hh
#define TsCircularPlasmid_hh

#include "TsVGeometryComponent.hh"

class TsCircularPlasmid : public TsVGeometryComponent
{    
public:
	TsCircularPlasmid(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
				  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	~TsCircularPlasmid();
	
	G4VPhysicalVolume* Construct();
	
	void ResolveParameters();
	
private:
	G4int fNumberOfBasePairs;
	G4double fRMin;
	G4double fRMax;
};

#endif
