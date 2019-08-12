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
#ifndef TsEosinophil_hh
#define TsEosinophil_hh

#include "TsVGeometryComponent.hh"


class TsEosinophil : public TsVGeometryComponent
{    
public:
	TsEosinophil(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
				  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	~TsEosinophil();
	
	G4VPhysicalVolume* Construct();
    
    void ResolveParameters();
    
private:
    G4double EosinophilRadius;
    G4double GranuleRadius;
};

#endif
