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

#ifndef TsNeuroMorpho_hh
#define TsNeuroMorpho_hh

#include "TsVGeometryComponent.hh"

class TsNeuroMorpho : public TsVGeometryComponent
{    
public:
	TsNeuroMorpho(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
				  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	~TsNeuroMorpho();
	
	G4VPhysicalVolume* Construct();
    void ResolveParameters();
    
private:
    
    G4int TotalParts;
    G4int SomaNumber, AxonNumber, BasalDendriteNumber, ApicalDendriteNumber;
    
    // NeuroMorpho filename
    G4String fNeuroMorphoFileName;
    G4double HLX, HLY, HLZ;
    
    G4double mxx, mxy, mxz, myx, myy, myz, mzx, mzy, mzz;
    G4double rx, ry, rz, r;
    G4double bx, by, bz;
    G4double mx, my, mz;
    
    std::vector<G4ThreeVector*> fPositions;
    std::vector<G4int> fPartID;
    std::vector<G4double> fRadius;
    std::vector<G4int> fID1;
    std::vector<G4int> fID2;
    std::vector<G4int> fSomaRadii;
    


};

#endif
