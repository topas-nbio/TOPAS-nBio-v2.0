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

#ifndef TsCylindericalNucleosome_hh
#define TsCylindericalNucleosome_hh

#include "TsVGeometryComponent.hh"

class TsCylindericalNucleosome : public TsVGeometryComponent
{    
public:
	TsCylindericalNucleosome(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
				  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	~TsCylindericalNucleosome();
	
	G4VPhysicalVolume* Construct();
};

#endif
