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

#ifndef TsMembraneRing_hh
#define TsMembraneRing_hh

#include "TsVGeometryComponent.hh"

class TsMembraneRing : public TsVGeometryComponent
{    
public:
	TsMembraneRing(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
				  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	~TsMembraneRing();
	
    void ResolveParameters();
    
	G4VPhysicalVolume* Construct();
private:
   // std::vector<G4VPhysicalVolume*> pLipidwrap;
   // std::vector<G4LogicalVolume*> lLipidwrap;
    
    G4int fNumberOfLipids;
    G4double fLipidHeadRadius;
    G4double fLipidTailLength;
    G4double fLipidTailRadius;
    
    G4double fLipidLength;
    G4double fRMin;
    G4double fRMax;
};

#endif
