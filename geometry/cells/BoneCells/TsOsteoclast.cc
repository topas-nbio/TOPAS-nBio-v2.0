// Component for TsOsteoclast
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
// Bone cells: Osteoclast
// Multi-nuclei 

#include "TsOsteoclast.hh"

#include "TsParameterManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"

#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Ellipsoid.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

TsOsteoclast::TsOsteoclast(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
			 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
    ResolveParameters();
}


TsOsteoclast::~TsOsteoclast()
{;}


void TsOsteoclast::ResolveParameters() {
    cell_radius = fPm->GetDoubleParameter(GetFullParmName("CellRadius"), "Length");
}



G4VPhysicalVolume* TsOsteoclast::Construct()
{
	BeginConstruction();

    //***********************************************************************
    //              Envelope Geometry : Cell
    //***********************************************************************
    
    G4String subComponentName1 = "Osteoclast";
    
    G4Orb* gCell1 = new G4Orb("gCell1", cell_radius);
    fEnvelopeLog = CreateLogicalVolume(subComponentName1, gCell1);
    fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);
    
    //***********************************************************************
    // Optional : include multiple nuclei in the cell
    //***********************************************************************
    
    G4String name = GetFullParmName("Nucleus/NumberOfNuclei");
    
    if (fPm->ParameterExists(name)) {
        
        const G4int NbOfNucl  = fPm->GetIntegerParameter( GetFullParmName("Nucleus/NumberOfNuclei") );
        
        G4double NuclRadius = 1.0*um;
        
        name = GetFullParmName("Nucleus/NucleusRadius");
        if (fPm->ParameterExists(name)) {NuclRadius = fPm->GetDoubleParameter(name, "Length");}
        
        G4String subComponentName2 = "Nuclei";
    
        G4Orb* gNucleus = new G4Orb("gNucleus", NuclRadius);
        G4LogicalVolume* lNucleus = CreateLogicalVolume(subComponentName2, gNucleus);
        
        //Randomly distribute nucleus throughout cell volume
        for (int j = 0; j < NbOfNucl; j++){
                
            G4bool Overlap = true;
            while (Overlap == true){
                    
                G4double u = G4UniformRand()*2*pi;
                G4double v = std::acos(2*G4UniformRand()-1);
                G4double dr = G4UniformRand()*(cell_radius);
                
                G4double x = 0.0;
                G4double y = 0.0;
                G4double z = 0.0;
                    
                x = (dr)* std::cos(u) * std::sin(v);
                y = (dr)* std::sin(u) * std::sin(v);
                z = (dr)* std::cos(v);
                
                G4ThreeVector* position = new G4ThreeVector(x,y,z);
                G4RotationMatrix* rotm = new G4RotationMatrix();
                    
                rotm->rotateX(0);
                rotm->rotateY(0);
                    
                G4VPhysicalVolume* pNucleus = CreatePhysicalVolume(subComponentName2, j, true, lNucleus, rotm, position, fEnvelopePhys);
                    
                G4bool OverlapCheck = pNucleus->CheckOverlaps();
                    
                if (OverlapCheck == false){break;}
                if (OverlapCheck == true){
					G4PhysicalVolumeStore::DeRegister(pNucleus);
                    G4cout << "**** Finding new position for volume " << subComponentName2 << ":" << j <<  " ****" << G4endl;
                }
            }
        }
    }
        

    InstantiateChildren(fEnvelopePhys);
	
	return fEnvelopePhys;
}
