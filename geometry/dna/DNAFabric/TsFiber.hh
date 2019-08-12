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

#ifndef TsFiber_hh
#define TsFiber_hh

#include "TsVGeometryComponent.hh"
#include "GeoManager.hh"


class TsFiber : public TsVGeometryComponent
{    
public:
	TsFiber(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
				  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	~TsFiber();
    
    GeoManager* fGeoManager;
	G4VPhysicalVolume* Construct();
    
    //void SetMaterial (G4String materialChoice);
    
    void UpdateGeometry();
    
private:
    
    G4Material *OtherMaterial(G4String materialName);
    G4double fWrapperRadius;
    G4double fWrapperHeight;

};

#endif
