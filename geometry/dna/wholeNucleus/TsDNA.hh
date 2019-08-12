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

#ifndef TsDNA_hh
#define TsDNA_hh

#include "TsVGeometryComponent.hh"

class TsDNA : public TsVGeometryComponent
{    
public:
	TsDNA(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
				  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	~TsDNA();
	
	G4VPhysicalVolume* Construct();
    
    void ResolveParameters();
    
private:
    
    G4bool fBuildChromatinFiber;
    G4bool fBuildBases;
    
    //Methods used by parametrisation
    void LoadChromosome(const char* filename,
                        G4VPhysicalVolume* chromoBox,
                        G4LogicalVolume* lFlower);
    
    G4LogicalVolume* lEnv;
    
    std::vector<G4ThreeVector*> fFlowPositions;
    std::vector<G4RotationMatrix*> fFlowRotations;
};

#endif
