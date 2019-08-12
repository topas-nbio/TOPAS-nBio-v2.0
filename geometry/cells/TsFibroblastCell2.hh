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

#ifndef TsFibroblastCell2_hh
#define TsFibroblastCell2_hh

#include "TsVGeometryComponent.hh"

class TsFibroblastCell2 : public TsVGeometryComponent
{    
public:
	TsFibroblastCell2(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
				  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	~TsFibroblastCell2();
	
	G4VPhysicalVolume* Construct();
};

#endif
