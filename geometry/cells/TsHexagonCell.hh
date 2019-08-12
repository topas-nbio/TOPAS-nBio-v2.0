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

#ifndef TsHexagonCell_hh
#define TsHexagonCell_hh

#include "TsVGeometryComponent.hh"

class TsHexagonCell : public TsVGeometryComponent
{    
public:
	TsHexagonCell(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
				  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	~TsHexagonCell();
	
	G4VPhysicalVolume* Construct();
};

#endif
