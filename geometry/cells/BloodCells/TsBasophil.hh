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
#ifndef TsBasophil_hh
#define TsBasophil_hh

#include "TsVGeometryComponent.hh"


class TsBasophil : public TsVGeometryComponent
{    
public:
	TsBasophil(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
				  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	~TsBasophil();
	
	G4VPhysicalVolume* Construct();
    
    void ResolveParameters();
    
private:
    G4double BasophilRadius;
    G4double GranuleRadius;
};

#endif
