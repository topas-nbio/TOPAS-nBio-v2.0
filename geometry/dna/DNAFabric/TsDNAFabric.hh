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

#ifndef TsDNAFabric_hh
#define TsDNAFabric_hh

#include "TsVGeometryComponent.hh"
#include "PhysGeoImport.hh"


class TsDNAFabric : public TsVGeometryComponent
{    
public:
	TsDNAFabric(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
				  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	~TsDNAFabric();
    
    //GeoManager* fGeoManager;
	G4VPhysicalVolume* Construct();
    

    
    //void UpdateGeometry();
    
private:
    
    G4Material *DefineMaterial(G4String materialName);
    G4double fFactor; ///< Scale factor: should be set to 1
    G4bool fbVisu;  ///< Enable or disable the Geant4 visualization. Slow.
    G4String fInputFile; ///< Geometry input file
};

#endif
