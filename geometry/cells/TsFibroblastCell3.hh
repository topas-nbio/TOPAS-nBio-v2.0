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

#ifndef TsFibroblastCell3_hh
#define TsFibroblastCell3_hh

#include "TsVGeometryComponent.hh"

class TsFibroblastCell3 : public TsVGeometryComponent
{    
public:
	TsFibroblastCell3(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
				  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	~TsFibroblastCell3();
	
	G4VPhysicalVolume* Construct();
};

#endif
