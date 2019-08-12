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

#ifndef TsOsteoclast_hh
#define TsOsteoclast_hh

#include "TsVGeometryComponent.hh"


class TsOsteoclast : public TsVGeometryComponent
{    
public:
	TsOsteoclast(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
				  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	~TsOsteoclast();
	
	G4VPhysicalVolume* Construct();
    
    void ResolveParameters();
    
private:
    G4double cell_radius;
};

#endif
