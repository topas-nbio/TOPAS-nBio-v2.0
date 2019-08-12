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

#ifndef TsChannel_hh
#define TsChannel_hh

#include "TsVGeometryComponent.hh"

class TsChannel : public TsVGeometryComponent
{    
public:
	TsChannel(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
				  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	~TsChannel();
	
	G4VPhysicalVolume* Construct();
    
    void ResolveParameters();
    
private:
    G4double fLipidHeadRadius;
    G4double fLipidTailLength;
    G4double fLipidTailRadius;
    
    G4double fLipidLength;
    G4int fNumberOfRows;
    G4int fNumberOfCols;
    
    G4double fChannelRadius;
    G4double fChannelInnerRadius;
    G4double fChannelLength;
};

#endif
