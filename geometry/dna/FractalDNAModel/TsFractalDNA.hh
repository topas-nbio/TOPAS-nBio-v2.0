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

#ifndef TsFractalDNA_hh
#define TsFractalDNA_hh

#include "TsVGeometryComponent.hh"
#include "GeoManager.hh"


class TsFractalDNA : public TsVGeometryComponent
{    
public:
	TsFractalDNA(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
				  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	~TsFractalDNA();
    
    GeoManager* fGeoManager;
	
    void ResolveParameters();
	G4VPhysicalVolume* Construct();
    
    std::vector<G4ThreeVector> GetSugar1Info();
    
private:
    
    std::vector<G4VPhysicalVolume*> pLoop;
    
    std::vector<G4double> fx;
    std::vector<G4double> fy;
    std::vector<G4double> fz;
    
    std::vector<G4LogicalVolume*> lLoop;
    
    G4bool fBuildBases;
    

    

};

#endif
