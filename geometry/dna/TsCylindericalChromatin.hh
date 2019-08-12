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

#ifndef TsCylindericalChromatin_hh
#define TsCylindericalChromatin_hh

#include "TsVGeometryComponent.hh"

class TsCylindericalChromatin : public TsVGeometryComponent
{    
public:
	TsCylindericalChromatin(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
				  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	~TsCylindericalChromatin();
	
	G4VPhysicalVolume* Construct();
};

#endif
