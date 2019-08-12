// Component for TsDNAFabric
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


#include "TsDNAFabric.hh"
#include "GeoManager.hh"

#include "TsParameterManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include "G4UserLimits.hh"
#include "G4NistManager.hh"
#include "G4RunManager.hh"
#include "G4Region.hh"

#include "G4UnionSolid.hh"

#include "G4SubtractionSolid.hh"

TsDNAFabric::TsDNAFabric(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
			 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
    
    fInputFile = "";
    
}


TsDNAFabric::~TsDNAFabric()
{;}


G4VPhysicalVolume* TsDNAFabric::Construct()
{
	BeginConstruction();

    fFactor = 100.0;
    fbVisu = true;
    
    //Read in user defined parameters
    
    G4String pname = GetFullParmName("DNAVis");
    if (fPm->ParameterExists(pname)) {
        fbVisu = fPm->GetBooleanParameter(pname);
    }
    
    //Default value is set to 100.
   pname = GetFullParmName("fFactor");
    if (fPm->ParameterExists(pname)) {
        fFactor = fPm->GetIntegerParameter(pname);
    }

    
    // Create the import class
    PhysGeoImport geo(fbVisu);
    // Set the scale factor (default setting = 1)
    geo.SetFactor(fFactor);
    
    // Wrapper volume for all contained geometry
    G4Box* sWrapper = new G4Box("solidWorld",12000.*fFactor*nm, 12000.*fFactor*nm, 12000.*fFactor*nm);
    fEnvelopeLog = CreateLogicalVolume(sWrapper);
    fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);
    
    // Load the fiber
    geo.CreateLogicVolumeDNA("fiber.dnafab");
    
	return fEnvelopePhys;
}

G4Material * TsDNAFabric::DefineMaterial(G4String materialName)
{
    G4Material * material(0);
    
    if(materialName == "G4_WATER"){
        // Water is defined from NIST material database
        G4NistManager * man = G4NistManager::Instance();
        material = man->FindOrBuildMaterial("G4_WATER");
        
    }
    
    return material;
}



