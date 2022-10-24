// Component for TsEosinophil
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
// White blood cell : Eosinophil.
// Nucleus is bilobed.

#include "TsEosinophil.hh"

#include "TsParameterManager.hh"
#include "G4VPhysicalVolume.hh"

#include "G4Orb.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "G4UnionSolid.hh"

TsEosinophil::TsEosinophil(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
			 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
    ResolveParameters();
}


TsEosinophil::~TsEosinophil()
{;}

void TsEosinophil::ResolveParameters() {
    EosinophilRadius = fPm->GetDoubleParameter(GetFullParmName("EosinophilRadius"), "Length");
    
    G4String name = GetFullParmName("Granule/GranuleRadius");
    if (fPm->ParameterExists(name)) {
        GranuleRadius = fPm->GetDoubleParameter(GetFullParmName("Granule/GranuleRadius"), "Length");
    }
    else {
        GranuleRadius = 0.25*micrometer;
    }
    
}


G4VPhysicalVolume* TsEosinophil::Construct()
{
	BeginConstruction();
        
    //***********************************************************************
    //              Envelope Geometry : spherical cell
    //***********************************************************************
    
    G4Orb* gEosinophil = new G4Orb(fName, EosinophilRadius);
    fEnvelopeLog = CreateLogicalVolume(gEosinophil);
    fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);
    
    
    //***********************************************************************
    //Nucleus is bilobed
    //***********************************************************************
    
    G4String subComponentName1 = "Nucleus";
    G4double Sphere1Radius = 2.0*um;
    G4double Sphere2Radius = 2.0*um;
    
    G4double NuclRadius = Sphere1Radius + Sphere2Radius;
    
    G4Orb* solidSphere1 = new G4Orb("Sphere1", Sphere1Radius);
    G4Orb* solidSphere2 = new G4Orb("Sphere2", Sphere2Radius);
    
    G4RotationMatrix* nuclRotation = new G4RotationMatrix();
    nuclRotation->rotateX(0.*deg);
    nuclRotation->rotateY(0.*deg);
    nuclRotation->rotateZ(0.*rad);
    
    G4ThreeVector* trans = new G4ThreeVector(-2.6*um,0,0);
    
    G4UnionSolid* gNucleus = new G4UnionSolid(subComponentName1, solidSphere1, solidSphere2, nuclRotation, G4ThreeVector(4*um,0,0));
    G4LogicalVolume* lNucleus = CreateLogicalVolume(subComponentName1, gNucleus);
    G4VPhysicalVolume* pNucleus = CreatePhysicalVolume(subComponentName1, lNucleus, nuclRotation, trans, fEnvelopePhys);
    
    G4bool OverlapCheck = pNucleus->CheckOverlaps();
    
    if (OverlapCheck == true){
        G4cerr << "Topas is exiting due to a serious error in geometry setup." << G4endl;
        G4cerr << "Nucleus overlaps with the cell." << G4endl;
        fPm->AbortSession(1);
    }
    

    //*****************************************
    // Subcomponent: Granules
    // Randomly distributed in cytoplasm
    //*****************************************
    
    G4String subComponentName2 = "Granule";
    
    G4String name = GetFullParmName("Granule/NumberOfGranules");
    
    if (fPm->ParameterExists(name)) {
        //number of granules
        G4int NumberOfGranules  = fPm->GetIntegerParameter( GetFullParmName("Granule/NumberOfGranules") );
    
    
        G4Orb* gGranule = new G4Orb(subComponentName2, GranuleRadius);
        G4LogicalVolume* lGranule = CreateLogicalVolume(subComponentName2, gGranule);
    
        //Randomly distribute granules throughout cytoplasm volume
        for (int j = 0; j < NumberOfGranules; j++){
        
            G4bool Overlap = true;
            while (Overlap == true){
            
                G4double u = G4UniformRand()*2*pi;
                G4double v = std::acos(2*G4UniformRand()-1);
                G4double dr = G4UniformRand()*(EosinophilRadius - NuclRadius);
            
                G4double x = 0.0;
                G4double y = 0.0;
                G4double z = 0.0;
            
                x = (NuclRadius + dr)* std::cos(u) * std::sin(v);
                y = (NuclRadius + dr)* std::sin(u) * std::sin(v);
                z = (NuclRadius + dr)* std::cos(v);
            
                G4ThreeVector* position = new G4ThreeVector(x,y,z);
            
                G4RotationMatrix* rotm = new G4RotationMatrix();
            
                rotm->rotateX(0);
                rotm->rotateY(0);
            
                G4VPhysicalVolume* pGranule = CreatePhysicalVolume(subComponentName2, j, true, lGranule, rotm, position, fEnvelopePhys);
            
                OverlapCheck = pGranule->CheckOverlaps();
            
                if (OverlapCheck == false){break;}
                if (OverlapCheck == true){
                    pGranule = NULL;
                    G4cout << "**** Finding new position for volume " << subComponentName2 << ":" << j <<  " ****" << G4endl;
                }
            }
        }
    }

    InstantiateChildren(fEnvelopePhys);
	
	return fEnvelopePhys;
}
