// Component for TsFiber
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
// Code that reads in the DNA fiber geometry from the DNAFabric software package tool
// More details about the fiber geometry can be found in:
// Meylan et al. (2016) Commputer Physics Comunications, 204, 159
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:


#include "TsFiber.hh"
#include "GeoManager.hh"

#include "TsParameterManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include "G4UserLimits.hh"
#include "G4NistManager.hh"
#include "G4RunManager.hh"
#include "G4Region.hh"

TsFiber::TsFiber(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
			 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
    
    fGeoManager = new GeoManager(0, 1.);
    
}


TsFiber::~TsFiber()
{
     delete fGeoManager;
}


G4VPhysicalVolume* TsFiber::Construct()
{
	BeginConstruction();
    
    fWrapperRadius = 22*nm;
    fWrapperHeight = 85*nm;
    
    //Wrapping component for whole fiber
    
    G4Tubs* sWrapper = new G4Tubs("solid_wrapper", 0., fWrapperRadius, fWrapperHeight, 0, 360);
    fEnvelopeLog = CreateLogicalVolume(sWrapper);
    fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);
    
    G4cout << "Building the Fiber" << G4endl;
    
    fGeoManager->Initialize();
    
    G4LogicalVolume* lFiber = fGeoManager->BuildLogicFiber(true);
    
    CreatePhysicalVolume("Fiber", lFiber, fEnvelopePhys);
    
	return fEnvelopePhys;
}

G4Material * TsFiber::OtherMaterial(G4String materialName)
{
    G4Material * material(0);
    
    if(materialName == "G4_WATER"){
        // Water is defined from NIST material database
        G4NistManager * man = G4NistManager::Instance();
        material = man->FindOrBuildMaterial("G4_WATER");

    }
    
    return material;
}




