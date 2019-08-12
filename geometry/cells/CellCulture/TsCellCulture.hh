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

#ifndef TsCellCulture_hh
#define TsCellCulture_hh

#include "TsVGeometryComponent.hh"


class TsCellCulture : public TsVGeometryComponent
{    
public:
	TsCellCulture(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
				  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	~TsCellCulture();
	
	G4VPhysicalVolume* Construct();
    
    void ResolveParameters();
    
private:
    G4double HLX;
    G4double HLY;
    G4double HLZ;
    
    G4double CellRadius;
    G4double NuclRadius;
    
    G4int NbOfCells;
};

#endif
