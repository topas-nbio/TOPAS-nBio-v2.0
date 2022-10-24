// Component for TsCharltonDNA
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
// DNA segment consisting of a central cylinder (base pair) and the sugarphosphate backbone
// The DNA backbone consists of two half cylinders wrapped around the base pair.
// Basepair length is 0.34 nm
// The sugar phosphate cylinders are rotated by 36 degrees to imitate the double helix structure of DNA
// Based on model in Charlton, Nikjoo & Humm (1989) Int J Radiat Biol 56(1), 1-19

#include "TsCharltonDNA.hh"

#include "TsParameterManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4TwoVector.hh"

#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"

#include "G4Box.hh"

TsCharltonDNA::TsCharltonDNA(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
                             TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
    ResolveParameters();
}


TsCharltonDNA::~TsCharltonDNA()
{;}

void TsCharltonDNA::ResolveParameters() {
    
    fNumberOfBasePairs = fPm->GetIntegerParameter(GetFullParmName("NumberOfBasePairs"));
    BasePairLength = fNumberOfBasePairs*0.17*nm;
    z0 = 0 - fNumberOfBasePairs*0.17*nm + 0.17*nm;
}


G4VPhysicalVolume* TsCharltonDNA::Construct()
{
    BeginConstruction();
    
    //****************************************************************************
    //              Envelope geometry: Box containing entire strand
    //****************************************************************************
    
    BoxHLX = 1.15*nm*2;
    BoxHLY = 1.15*nm*2;
    BoxHLZ = fNumberOfBasePairs*0.17*nm*2;
    
    G4Box* envelopeSolid = new G4Box("box", BoxHLX, BoxHLY, BoxHLZ);
    fEnvelopeLog = CreateLogicalVolume(envelopeSolid);
    fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);
    
    //****************************************************************************
    //              Subcomponent 1: DNA base-pair
    //****************************************************************************
    
    //Dimensions are taken from Charlton et al. (1989)
    
    G4String subComponentName1 = "BasePair";
    
    G4Tubs* gBasePair = new G4Tubs(subComponentName1, 0, 0.5*nm, BasePairLength, 0.0*deg, 360.0*deg);
    G4LogicalVolume* lBasePair = CreateLogicalVolume(subComponentName1, gBasePair);
    CreatePhysicalVolume(subComponentName1, lBasePair, fEnvelopePhys);
    
    //****************************************************************************
    //              Subcomponent 2: DNA sugar-phosphate backbone strands
    //****************************************************************************
    
    //Strand 1
    
    G4String subComponentName2 = "Backbone1";
    G4Tubs* gSugar1 = new G4Tubs(subComponentName2, 0.5*nm, 1.15*nm, 0.17*nm, 0*deg, 180*deg);
    G4LogicalVolume* lSugar1 = CreateLogicalVolume(subComponentName2, gSugar1);
    
    
    //Strand 2
    
    G4String subComponentName3 = "Backbone2";
    G4Tubs* gSugar2 = new G4Tubs(subComponentName3, 0.5*nm, 1.15*nm, 0.17*nm, 180*deg, 180*deg);
    G4LogicalVolume* lSugar2 = CreateLogicalVolume(subComponentName3, gSugar2);
    
    
    //****************************************************************************
    //              Rotation of strands around the base pair
    //****************************************************************************
    
    G4double x = 0.0;
    G4double y = 0.0;
    
    for (int j = 0; j < fNumberOfBasePairs; j++){
        
        G4double theta = 36*deg*j;
        G4double z = z0 + j*0.34*nm;
        
        G4ThreeVector* position = new G4ThreeVector(x, y, z);
        G4RotationMatrix* rot = new G4RotationMatrix();
        rot ->rotateZ(theta);
        
        CreatePhysicalVolume(subComponentName2, j, true, lSugar1, rot, position, fEnvelopePhys);
        CreatePhysicalVolume(subComponentName3, j, true, lSugar2, rot, position, fEnvelopePhys);
    }
    
    
    InstantiateChildren(fEnvelopePhys);
    
    return fEnvelopePhys;
}
