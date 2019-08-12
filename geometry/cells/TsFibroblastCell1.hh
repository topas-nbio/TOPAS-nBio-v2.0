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

#ifndef TsFibroblastCell1_hh
#define TsFibroblastCell1_hh

#include "TsVGeometryComponent.hh"

class TsFibroblastCell1 : public TsVGeometryComponent
{    
public:
	TsFibroblastCell1(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
				  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	~TsFibroblastCell1();
	
	G4VPhysicalVolume* Construct();
};

#endif
