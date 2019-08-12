// Component for TsCylindericalDNA
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

#include "TsCylindericalDNA.hh"

#include "TsParameterManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4Tubs.hh"

TsCylindericalDNA::TsCylindericalDNA(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
			 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{;}


TsCylindericalDNA::~TsCylindericalDNA()
{;}


G4VPhysicalVolume* TsCylindericalDNA::Construct()
{
	BeginConstruction();

    G4String name = GetFullParmName("DNARadius");
    G4double DNARadius;
    if (!fPm->ParameterExists(name)) {
        DNARadius = 1*nm;
    }
    else{
       DNARadius = fPm->GetDoubleParameter(name, "Length");
    }
    
    name = GetFullParmName("DNAHalfLength");
    G4double DNALength;
    if (!fPm->ParameterExists(name)) {
        DNALength = 1*nm;
    }
    else{
        DNALength = fPm->GetDoubleParameter(name, "Length");
    }
    
	G4Tubs* envelopeSolid = new G4Tubs(fName, 0, DNARadius, DNALength,0.0*deg, 360.0*deg);
	fEnvelopeLog = CreateLogicalVolume(envelopeSolid);
	fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);
    
	InstantiateChildren(fEnvelopePhys);
	
	return fEnvelopePhys;
}
