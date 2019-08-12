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
#ifndef TsRedBloodCell_hh
#define TsRedBloodCell_hh

#include "TsVGeometryComponent.hh"


class TsRedBloodCell : public TsVGeometryComponent
{    
public:
	TsRedBloodCell(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
				  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	~TsRedBloodCell();
	
	G4VPhysicalVolume* Construct();
    
    void ResolveParameters();
    
private:
    G4double RedBloodCellRadius;
    G4double TorusWidth;
};

#endif
