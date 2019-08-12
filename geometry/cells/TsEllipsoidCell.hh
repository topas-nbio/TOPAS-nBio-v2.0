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

#ifndef TsEllipsoidCell_hh
#define TsEllipsoidCell_hh

#include "TsVGeometryComponent.hh"


class TsEllipsoidCell : public TsVGeometryComponent
{    
public:
	TsEllipsoidCell(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
				  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	~TsEllipsoidCell();
	
	G4VPhysicalVolume* Construct();

    void ResolveParameters();
    
private:
    G4double xSA;
    G4double ySA;
    G4double zSA;
};

#endif
