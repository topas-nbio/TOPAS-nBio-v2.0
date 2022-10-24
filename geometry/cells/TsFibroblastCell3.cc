// Component for TsFibroblastCell3
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
// Fibroblast cell.
// User has the option of including organelles: nucleus and/or mitochondria.

#include "TsFibroblastCell3.hh"

#include "TsParameterManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"

#include "G4TwoVector.hh"
#include "G4ExtrudedSolid.hh"
#include "G4Orb.hh"
#include "G4Ellipsoid.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"


TsFibroblastCell3::TsFibroblastCell3(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
			 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{;}


TsFibroblastCell3::~TsFibroblastCell3()
{;}


G4VPhysicalVolume* TsFibroblastCell3::Construct()
{
	BeginConstruction();

    //***********************************************************************
    //              Envelope Geometry : Cell
    //***********************************************************************
    
    //Base of the cell:
    std::vector<G4TwoVector> cellpoly(111);
    cellpoly[0] = G4TwoVector(25*micrometer, 0*micrometer);
    cellpoly[1] = G4TwoVector(23.5*micrometer, 1*micrometer);
    cellpoly[2] = G4TwoVector(22*micrometer, 2*micrometer);
    cellpoly[3] = G4TwoVector(21*micrometer, 3*micrometer);
    cellpoly[4] = G4TwoVector(20.5*micrometer, 4*micrometer);
    cellpoly[5] = G4TwoVector(22*micrometer,8*micrometer);
    cellpoly[6] = G4TwoVector(23.5*micrometer,9*micrometer);
    cellpoly[7] = G4TwoVector(24*micrometer,10*micrometer);
    cellpoly[8] = G4TwoVector(22*micrometer,11*micrometer);
    cellpoly[9] = G4TwoVector(21*micrometer,11.5*micrometer);
    cellpoly[10] = G4TwoVector(19.5*micrometer,11.8*micrometer);
    cellpoly[11] = G4TwoVector(18*micrometer,12*micrometer);
    cellpoly[12] = G4TwoVector(17*micrometer,12.5*micrometer);
    cellpoly[13] = G4TwoVector(15*micrometer,13*micrometer);
    cellpoly[14] = G4TwoVector(13*micrometer,15*micrometer);
    cellpoly[15] = G4TwoVector(14*micrometer,15.5*micrometer);
    cellpoly[16] = G4TwoVector(14.5*micrometer,17*micrometer);
    cellpoly[17] = G4TwoVector(15*micrometer,18*micrometer);
    cellpoly[18] = G4TwoVector(15.5*micrometer,19*micrometer);
    cellpoly[19] = G4TwoVector(16*micrometer,20*micrometer);
    cellpoly[20] = G4TwoVector(13*micrometer,19.5*micrometer);
    cellpoly[22] = G4TwoVector(11*micrometer,19*micrometer);
    cellpoly[21] = G4TwoVector(9.5*micrometer,18.5*micrometer);
    cellpoly[22] = G4TwoVector(7.5*micrometer,19*micrometer);
    cellpoly[23] = G4TwoVector(5.5*micrometer,19*micrometer);
    cellpoly[24] = G4TwoVector(4*micrometer,20*micrometer);
    cellpoly[25] = G4TwoVector(2*micrometer,21*micrometer);
    cellpoly[26] = G4TwoVector(0*micrometer,22*micrometer);
    cellpoly[27] = G4TwoVector(-3*micrometer,20*micrometer);
    cellpoly[28] = G4TwoVector(-5*micrometer,19.5*micrometer);
    cellpoly[29] = G4TwoVector(-7*micrometer,20*micrometer);
    cellpoly[30] = G4TwoVector(-8.5*micrometer,21*micrometer);
    cellpoly[31] = G4TwoVector(-10*micrometer,21.5*micrometer);
    cellpoly[32] = G4TwoVector(-12.5*micrometer,22*micrometer);
    cellpoly[33] = G4TwoVector(-14*micrometer,24*micrometer);
    cellpoly[34] = G4TwoVector(-15*micrometer,25*micrometer);
    cellpoly[35] = G4TwoVector(-15*micrometer,22*micrometer);
    cellpoly[36] = G4TwoVector(-15.5*micrometer,19*micrometer);
    cellpoly[37] = G4TwoVector(-16*micrometer,17*micrometer);
    cellpoly[38] = G4TwoVector(-16.5*micrometer,16*micrometer);
    cellpoly[39] = G4TwoVector(-17*micrometer,14.5*micrometer);
    cellpoly[40] = G4TwoVector(-18*micrometer,14*micrometer);
    cellpoly[41] = G4TwoVector(-21*micrometer,13*micrometer);
    cellpoly[42] = G4TwoVector(-20*micrometer,12*micrometer);
    cellpoly[43] = G4TwoVector(-19*micrometer,11*micrometer);
    cellpoly[44] = G4TwoVector(-18*micrometer,9*micrometer);
    cellpoly[45] = G4TwoVector(-17.5*micrometer,8*micrometer);
    cellpoly[46] = G4TwoVector(-17.3*micrometer,7*micrometer);
    cellpoly[47] = G4TwoVector(-17.2*micrometer,6*micrometer);
    cellpoly[48] = G4TwoVector(-17*micrometer,5*micrometer);
    cellpoly[49] = G4TwoVector(-16.5*micrometer,4*micrometer);
    cellpoly[50] = G4TwoVector(-16*micrometer,0*micrometer);
    cellpoly[51] = G4TwoVector(-16.5*micrometer,-4*micrometer);
    cellpoly[52] = G4TwoVector(-17*micrometer,-6*micrometer);
    cellpoly[53] = G4TwoVector(-16.5*micrometer,-8*micrometer);
    cellpoly[54] = G4TwoVector(-19*micrometer,-11*micrometer);
    cellpoly[55] = G4TwoVector(-20*micrometer,-12*micrometer);
    cellpoly[56] = G4TwoVector(-22*micrometer,-15*micrometer);
    cellpoly[57] = G4TwoVector(-20*micrometer,-14*micrometer);
    cellpoly[58] = G4TwoVector(-18*micrometer,-12*micrometer);
    cellpoly[59] = G4TwoVector(-17*micrometer,-10*micrometer);
    cellpoly[60] = G4TwoVector(-15*micrometer,-9*micrometer);
    cellpoly[61] = G4TwoVector(-14*micrometer,-8*micrometer);
    cellpoly[62] = G4TwoVector(-13*micrometer,-10*micrometer);
    cellpoly[63] = G4TwoVector(-12*micrometer,-11*micrometer);
    cellpoly[64] = G4TwoVector(-11*micrometer,-14*micrometer);
    cellpoly[65] = G4TwoVector(-11.5*micrometer,-15.5*micrometer);
    cellpoly[66] = G4TwoVector(-11.5*micrometer,-18*micrometer);
    cellpoly[67] = G4TwoVector(-11.5*micrometer,-21*micrometer);
    cellpoly[68] = G4TwoVector(-10*micrometer,-19*micrometer);
    cellpoly[69] = G4TwoVector(-9*micrometer,-17*micrometer);
    cellpoly[70] = G4TwoVector(-8*micrometer,-14.5*micrometer);
    cellpoly[71] = G4TwoVector(-7*micrometer,-13*micrometer);
    cellpoly[72] = G4TwoVector(-5.5*micrometer,-11*micrometer);
    cellpoly[73] = G4TwoVector(-4*micrometer,-10.5*micrometer);
    cellpoly[74] = G4TwoVector(-2*micrometer,-10*micrometer);
    cellpoly[75] = G4TwoVector(0*micrometer,-10*micrometer);
    cellpoly[76] = G4TwoVector(2*micrometer,-11*micrometer);
    cellpoly[77] = G4TwoVector(3*micrometer,-12*micrometer);
    cellpoly[78] = G4TwoVector(5*micrometer,-14*micrometer);
    cellpoly[79] = G4TwoVector(7*micrometer,-15*micrometer);
    cellpoly[80] = G4TwoVector(8*micrometer,-17*micrometer);
    cellpoly[81] = G4TwoVector(8.5*micrometer,-18*micrometer);
    cellpoly[82] = G4TwoVector(7*micrometer,-21*micrometer);
    cellpoly[83] = G4TwoVector(8*micrometer,-22*micrometer);
    cellpoly[84] = G4TwoVector(7.5*micrometer,-24*micrometer);
    cellpoly[85] = G4TwoVector(7*micrometer,-26*micrometer);
    cellpoly[86] = G4TwoVector(10*micrometer,-25*micrometer);
    cellpoly[87] = G4TwoVector(11*micrometer,-24.5*micrometer);
    cellpoly[88] = G4TwoVector(12*micrometer,-22*micrometer);
    cellpoly[89] = G4TwoVector(13*micrometer,-20*micrometer);
    cellpoly[90] = G4TwoVector(14*micrometer,-18*micrometer);
    cellpoly[91] = G4TwoVector(15*micrometer,-16*micrometer);
    cellpoly[92] = G4TwoVector(15.5*micrometer,-14*micrometer);
    cellpoly[93] = G4TwoVector(16*micrometer,-13*micrometer);
    cellpoly[94] = G4TwoVector(18*micrometer,-13*micrometer);
    cellpoly[95] = G4TwoVector(19*micrometer,-14.5*micrometer);
    cellpoly[96] = G4TwoVector(20*micrometer,-16*micrometer);
    cellpoly[97] = G4TwoVector(21*micrometer,-18*micrometer);
    cellpoly[98] = G4TwoVector(21.5*micrometer,-19*micrometer);
    cellpoly[99] = G4TwoVector(22*micrometer,-18*micrometer);
    cellpoly[100] = G4TwoVector(21*micrometer,-16*micrometer);
    cellpoly[101] = G4TwoVector(20*micrometer,-14*micrometer);
    cellpoly[102] = G4TwoVector(19.4*micrometer,-11.5*micrometer);
    cellpoly[103] = G4TwoVector(19*micrometer,-11*micrometer);
    cellpoly[104] = G4TwoVector(19.4*micrometer,-9*micrometer);
    cellpoly[105] = G4TwoVector(19.8*micrometer,-6*micrometer);
    cellpoly[106] = G4TwoVector(20*micrometer,-5*micrometer);
    cellpoly[107] = G4TwoVector(21*micrometer,-4*micrometer);
    cellpoly[108] = G4TwoVector(22*micrometer,-3*micrometer);
    cellpoly[109] = G4TwoVector(23*micrometer,-2*micrometer);
    cellpoly[110] = G4TwoVector(24*micrometer,-1*micrometer);
    
    // Height of cell (z):
    G4double hz = 10*micrometer;  //half length along z

    // Used for defining the area that mitochondria are distributed in
    G4double CellRadius = 10*micrometer;
    
    G4ExtrudedSolid* gFibroCell = new G4ExtrudedSolid(fName,
                                                      cellpoly,
                                                      hz,
                                                      G4TwoVector(0,0), 0.5, G4TwoVector(0,0), 1.0); // shape along z
    
    fEnvelopeLog = CreateLogicalVolume(gFibroCell);
    fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);
    
    //***********************************************************************
    // Optional : include a nucleus and/or mitochondria in the cell
    //***********************************************************************
    
    //***************************
    // Subcomponent: Nucleus
    //***************************
    
    G4double NuclRadius = 0.0*um;
    G4String name = GetFullParmName("Nucleus/NucleusRadius");
    if (fPm->ParameterExists(name)) {
        NuclRadius = fPm->GetDoubleParameter(name, "Length");
  
        G4String subComponentName1 = "Nucleus";
        G4Orb* gNucleus = new G4Orb("gNucleus", NuclRadius);
        G4LogicalVolume* lNucleus = CreateLogicalVolume(subComponentName1, gNucleus);
        
        G4VPhysicalVolume* pNucleus = CreatePhysicalVolume(subComponentName1, lNucleus, fEnvelopePhys);
        
        G4bool OverlapCheck = pNucleus->CheckOverlaps();
        
        if (OverlapCheck == true){
            G4cerr << "Topas is exiting due to a serious error in geometry setup." << G4endl;
            G4cerr << "Nucleus overlaps with the cell." << G4endl;
            fPm->AbortSession(1);
        }
    }
    
    
    //*******************************
    // Subcomponent: Mitochondria
    //*******************************
    name = GetFullParmName("Mitochondria/NumberOfMitochondria");
    if (fPm->ParameterExists(name)) {
        
        //number of mitochondria
        const G4int NbOfMito  = fPm->GetIntegerParameter( GetFullParmName("Mitochondria/NumberOfMitochondria") );
        
        //Semi-axis lengths of the ellpsoid
        G4double EllA = 0.5*micrometer;
        G4double EllB = 0.3*micrometer;
        G4double EllC = 0.9*micrometer;
        
        name=GetFullParmName("Mitochondria/a");
        if (fPm->ParameterExists(name)){EllA = fPm->GetDoubleParameter(GetFullParmName("Mitochondria/a"), "Length" );}
        
        name=GetFullParmName("Mitochondria/b");
        if (fPm->ParameterExists(name)){EllB = fPm->GetDoubleParameter(GetFullParmName("Mitochondria/b"), "Length" );}
        
        name=GetFullParmName("Mitochondria/c");
        if (fPm->ParameterExists(name)){EllC = fPm->GetDoubleParameter(GetFullParmName("Mitochondria/c"), "Length" );}
        
        
        G4String subComponentName2 = "Mitochondria";
        G4Ellipsoid* gMito = new G4Ellipsoid("gMito", EllA, EllB, EllC);
        G4LogicalVolume* lMito = CreateLogicalVolume(subComponentName2, gMito);
    
        //Randomly distribute mitochondria throughout cell volume outside nucleus (default)
        for (int j = 0; j < NbOfMito; j++){
    
            G4bool Overlap = true;
            while (Overlap == true){
        
                G4double u = G4UniformRand()*2*pi;
                G4double v = std::acos(2*G4UniformRand()-1);
                G4double dr = G4UniformRand()*(CellRadius - NuclRadius);
                G4double phi = G4UniformRand()*2*pi;
                G4double psi = G4UniformRand()*2*pi;
                G4double x = 0.0;
                G4double y = 0.0;
                G4double z = 0.0;
    
                x = (NuclRadius + dr)* std::cos(u) * std::sin(v);
                y = (NuclRadius + dr)* std::sin(u) * std::sin(v);
                z = (NuclRadius + dr)* std::cos(v);

                G4ThreeVector* position = new G4ThreeVector(x,y,z);
                
                G4RotationMatrix* rotm = new G4RotationMatrix();
                
                rotm->rotateX(psi);
                rotm->rotateY(phi);
                
                G4VPhysicalVolume* pMito = CreatePhysicalVolume(subComponentName2, j, true, lMito, rotm, position, fEnvelopePhys);
                
                G4bool OverlapCheck = pMito->CheckOverlaps();
                
                if (OverlapCheck == false){break;}
                if (OverlapCheck == true){
                    G4PhysicalVolumeStore::DeRegister(pMito);
                    G4cout << "**** Finding new position for volume " << subComponentName2 << ":" << j <<  " ****" << G4endl;
                }
            }
        }
    }
    
	InstantiateChildren(fEnvelopePhys);
	
	return fEnvelopePhys;
}
