// Component for TsDNA
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
//Geometry is based on the extended Geant4 example: wholeNuclearDNA
//This example is provided by the Geant4-DNA collaboration. Any report or published results obtained using
// this DNA geometry shall cite the following Geant4-DNA collaboration publications:
//[1] NIM B 298 (2013) 47-54
//[2] Med. Phys. 37 (2010) 4692-4708 [3] Phys. Med. 31 (2015) 861-874

#include "TsDNA.hh"

#include "TsParameterManager.hh"
#include "G4VPhysicalVolume.hh"

#include "G4UnionSolid.hh"
#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Tubs.hh"
#include "G4Ellipsoid.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

#include "G4VPVParameterisation.hh"
#include "ChromosomeParameterisation.hh"

#include <fstream>

#define countof(x) (sizeof(x) / sizeof(x[0]))

using namespace std;
using CLHEP::mm;
using CLHEP::degree;
using CLHEP::nanometer;
using CLHEP::micrometer;

#include "globals.hh"
#include "CLHEP/Units/SystemOfUnits.h"
#include "G4LogicalVolume.hh"
#include "G4UnionSolid.hh"

#include "G4VisAttributes.hh"
#include "G4NistManager.hh"


TsDNA::TsDNA(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
             TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
    ResolveParameters();
}

TsDNA::~TsDNA()
{;}

void TsDNA::ResolveParameters(){
    fBuildChromatinFiber = fPm->GetBooleanParameter(GetFullParmName("BuildChromatinFiber"));
    fBuildBases = fPm->GetBooleanParameter(GetFullParmName("BuildBases"));
    
    if (fBuildChromatinFiber) {
        G4cout << "Topas is building Chromatin fibres." << G4endl;
        G4cout << "This can slow down simulations with visualization" << G4endl;
    }
    
    if (fBuildBases){
        G4cout << "Topas is building Bases." << G4endl;
        G4cout << "This can slow down simulations with visualization" << G4endl;
    }
}


void TsDNA::LoadChromosome(const char* filename,
                            G4VPhysicalVolume* chromBox,
                            G4LogicalVolume* lFlower)
{
    ChromosomeParameterisation* cp = new ChromosomeParameterisation(filename);
    new G4PVParameterised("box ros",
                          lFlower,
                          chromBox,
                          kUndefined,
                          cp->GetNumRosettes(),
                          cp);

    G4cout << filename << " done" << G4endl;
}


G4VPhysicalVolume* TsDNA::Construct()
{
	BeginConstruction();

    /****************************************************************************/
    //                             Box (envelope)
    /****************************************************************************/

    G4Box* gTin = new G4Box("gTin",
                            13 * micrometer,
                            10 * micrometer,
                            5 * micrometer);

	fEnvelopeLog = CreateLogicalVolume(gTin);
	fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);

    /****************************************************************************/
    //                           Subcomponent 1: Cell nucleus
    /****************************************************************************/
    G4String subComponentName1 = "Nucleus";
    G4Ellipsoid* gNucleus = new G4Ellipsoid("gNucleus",
                                            13 * micrometer,  //11.83 micrometer
                                            10 * micrometer,  // 8.52 micrometer
                                            4 * micrometer,   // 3 micrometer
                                            0,
                                            0);

    G4LogicalVolume* lNucleus = CreateLogicalVolume(subComponentName1, gNucleus);

    G4VPhysicalVolume* pNucleus = CreatePhysicalVolume(subComponentName1, lNucleus, fEnvelopePhys);

    /****************************************************************************/
    //                        Subcomponent 2 : Chromosomes territories
    /****************************************************************************/
    // NOTE: The only supported values for the rotation are
    // 0 and 90 degrees on the Y axis.
    G4double chromosomePositionSizeRotation[][7] = {
        {4.467, 2.835, 0, 1.557, 1.557, 1.557, 90},
        {-4.467, 2.835, 0, 1.557, 1.557, 1.557, 0},
        {4.423, -2.831, 0, 1.553, 1.553, 1.553, 90},
        {-4.423, -2.831, 0, 1.553, 1.553, 1.553, 0},
        {1.455, 5.63, 0, 1.455, 1.455, 1.455, 0},
        {-1.455, 5.63, 0, 1.455, 1.455, 1.455, 90},
        {1.435, 0, 1.392, 1.435, 1.435, 1.435, 0},
        {-1.435, 0, 1.392, 1.435, 1.435, 1.435, 90},
        {1.407, 0, -1.450, 1.407, 1.407, 1.407, 90}, // 5 right
        {-1.407, 0, -1.450, 1.407, 1.407, 1.407, 0}, // 5 left
        {1.380, -5.437, 0, 1.380, 1.380, 1.380, 0},
        {-1.380, -5.437, 0, 1.380, 1.380, 1.380, 90},
        {1.347, 2.782, -1.150, 1.347, 1.347, 1.347, 90},
        {-1.347, 2.782, -1.150, 1.347, 1.347, 1.347, 0},
        {1.311, -2.746, -1.220, 1.311, 1.311, 1.311, 90},
        {-1.311, -2.746, -1.220, 1.311, 1.311, 1.311, 0},
        {7.251, -2.541, 0, 1.275, 1.275, 1.275, 0},
        {-6.701, 0, -0.85, 1.275, 1.275, 1.275, 90},
        {4.148, 0, 1.278, 1.278, 1.278, 1.278, 90}, // 10 right
        {-4.148, 0, 1.278, 1.278, 1.278, 1.278, 0}, // 10 left
        {4.147, 0, -1.277, 1.277, 1.277, 1.277, 0},
        {-4.147, 0, -1.277, 1.277, 1.277, 1.277, 90},
        {8.930, 0.006, 0, 1.272, 1.272, 1.272, 90},
        {-7.296, 2.547, 0, 1.272, 1.272, 1.272, 90},
        {1.207, -2.642, 1.298, 1.207, 1.207, 1.207, 0},
        {-1.207, -2.642, 1.298, 1.207, 1.207, 1.207, 90},
        {1.176, 2.611, 1.368, 1.176, 1.176, 1.176, 0},
        {-1.176, 2.611, 1.368, 1.176, 1.176, 1.176, 90},
        {4.065, 5.547, 0, 1.155, 1.155, 1.155, 90}, // 15 right
        {-4.065, 5.547, 0, 1.155, 1.155, 1.155, 0}, // 15 left
        {6.542, 0.159, 1.116, 1.116, 1.116, 1.116, 0},
        {-9.092, 0, 0, 1.116, 1.116, 1.116, 0},
        {6.507, 0.159, -1.081, 1.081, 1.081, 1.081, 90},
        {-7.057, -2.356, 0, 1.081, 1.081, 1.081, 90},
        {3.824, -5.448, 0, 1.064, 1.064, 1.064, 90},
        {-3.824, -5.448, 0, 1.064, 1.064, 1.064, 0},
        {5.883, -5.379, 0, 0.995, 0.995, 0.995, 0},
        {-9.133, -2.111, 0, 0.995, 0.995, 0.995, 0},
        {6.215, 5.387, 0, 0.995, 0.995, 0.995, 0}, // 20 right
        {-6.971, -4.432, 0, 0.995, 0.995, 0.995, 90}, // 20 left
        {9.583, 2.177, 0, 0.899, 0.899, 0.899, 90},
        {-9.467, 2.03, 0, 0.899, 0.899, 0.899, 0},
        {9.440, -2.180, 0, 0.914, 0.914, 0.914, 90},
        {-6.34, 0, 1.339, 0.914, 0.914, 0.914, 0},
        {-6.947, 4.742, 0, 0.923, 0.923, 0.923, 90}, // Y
        {7.354, 2.605, 0, 1.330, 1.330, 1.330, 0} // X
    };

    G4RotationMatrix* rotch = new G4RotationMatrix;
    rotch->rotateY(90 * degree);

    vector<G4VPhysicalVolume*> pChromoTerr(48);

    for (unsigned int i = 0; i < countof(chromosomePositionSizeRotation); i++)
    {
        G4double* p = &chromosomePositionSizeRotation[i][0];
        G4double* size = &chromosomePositionSizeRotation[i][3];
        G4double rotation = chromosomePositionSizeRotation[i][6];
        G4ThreeVector* pos = new G4ThreeVector(p[0] * micrometer, p[1] * micrometer, p[2] * micrometer);
        G4RotationMatrix* rot = rotation == 0 ? 0 : rotch;

        ostringstream ss;
        ss << "box" << (i / 2) + 1 << (i % 2 ? 'l' : 'r');
        G4String name = ss.str();
        ss.str("");
        ss.clear();

        G4String subComponentName3 = "chromoTerr";
        G4Box* gChromoTerr = new G4Box(name,
                                    size[0] * micrometer,
                                    size[1] * micrometer,
                                    size[2] * micrometer);

        G4LogicalVolume* lChromoTerr = CreateLogicalVolume(subComponentName3,gChromoTerr);

        pChromoTerr[i] = CreatePhysicalVolume(subComponentName3,lChromoTerr, rot, pos, pNucleus);


    }

    /**************************************************************************/
    //                 Subcomponent 4: Box containing the chromatin flowers
    /**************************************************************************/

    G4String subComponentName4 = "flower";
    G4Tubs* gFlower = new G4Tubs(subComponentName4,
                                     0 * nanometer,
                                     399 * nanometer,
                                     20 * nanometer,
                                     0 * degree,
                                     360 * degree);

    G4LogicalVolume* lFlower = CreateLogicalVolume(subComponentName4,gFlower);

    //Loading flower box position for each chromosome territory
    G4String name;


    for (int k = 0; k < 22; k++)
    {
        ostringstream oss;
        oss << "chromo" << k + 1 << ".dat";
        name = oss.str();
        oss.str("");
        oss.clear();

        LoadChromosome(name.c_str(), pChromoTerr[k * 2], lFlower);
        LoadChromosome(name.c_str(), pChromoTerr[k * 2 + 1], lFlower);

    }

    LoadChromosome("chromoY.dat", pChromoTerr[44], lFlower);
    LoadChromosome("chromoX.dat", pChromoTerr[45], lFlower);


    /**************************************************************************/
    //                 Subcomponent 5: chromatin fibers
    /**************************************************************************/
if (fBuildChromatinFiber)
   {
        // chromatin fiber envelope
        G4String subComponentName5 = "chromatinFiber";
        G4Tubs* gEnv = new G4Tubs(subComponentName5,
                                  0,
                                  15.4 * nanometer,
                                  80.5 * nanometer,
                                  0 * degree,
                                  360 * degree);

       lEnv = CreateLogicalVolume(subComponentName5,
                                                    gEnv);


        // Chromatin fiber position
        for (G4int i = 0; i < 7; i++)
        {
            G4RotationMatrix* rotFiber = new G4RotationMatrix;
            rotFiber->rotateX(90 * degree);
            rotFiber->rotateY(i * 25.72 * degree);
            G4ThreeVector posFiber = G4ThreeVector(0, 152 * nanometer, 0);
            posFiber.rotateZ(i * 25.72 * degree);
            new G4PVPlacement(rotFiber,
                              posFiber,
                              lEnv,
                              "physi env",
                              lFlower,
                              false,
                              0);

            rotFiber = new G4RotationMatrix;
            rotFiber->rotateX(90 * degree);
            rotFiber->rotateY((7 + i) * 25.72 * degree);
            posFiber = G4ThreeVector(0, 152 * nanometer, 0);
            posFiber.rotateZ((7 + i) * 25.72 * degree);
            new G4PVPlacement(rotFiber,
                              posFiber,
                              lEnv,
                              "physi env",
                              lFlower,
                              false,
                              0);

            rotFiber = new G4RotationMatrix;
            rotFiber->rotateX(90 * degree);
            rotFiber->rotateY((25.72 + (i - 14) * 51.43) * degree);
            posFiber = G4ThreeVector(-36.5 * nanometer, 312 * nanometer, 0);
            posFiber.rotateZ((i - 14) * 51.43 * degree);
            new G4PVPlacement(rotFiber,
                              posFiber,
                              lEnv,
                              "physi env",
                              lFlower,
                              false,
                              0);

            rotFiber = new G4RotationMatrix;
            rotFiber->rotateX(90 * degree);
            rotFiber->rotateY(180 * degree);
            rotFiber->rotateY((i - 21) * 51.43 * degree);
            posFiber = G4ThreeVector(-103 * nanometer, 297 * nanometer, 0);
            posFiber.rotateZ((i - 21) * 51.43 * degree);
            new G4PVPlacement(rotFiber,
                              posFiber,
                              lEnv,
                              "physi env",
                              lFlower,
                              false,
                              0);

        }

    }


    /**************************************************************************/
    //                 Subcomponent 6: Bases
    /**************************************************************************/

    if (fBuildBases)
    {
        // Histones
        G4String Subcomponent6a = "Histone";
        G4Tubs* sHistone = new G4Tubs(Subcomponent6a,
                                      0,
                                      3.25 * nanometer,
                                      2.85 * nanometer,
                                      0 * degree,
                                      360 * degree);
        G4LogicalVolume* lHistone = CreateLogicalVolume(Subcomponent6a, sHistone);

        //Base pair
        G4String Subcomponent6b = "BasePair1";
        G4String Subcomponent6c = "BasePair2";
        G4Orb* sBp1 = new G4Orb(Subcomponent6b, 0.17 * nanometer);
        G4LogicalVolume* lBp1 = CreateLogicalVolume(Subcomponent6b, sBp1);

        G4Orb* sBp2 = new G4Orb(Subcomponent6c, 0.17 * nanometer);
        G4LogicalVolume* lBp2 = CreateLogicalVolume(Subcomponent6b, sBp2);

        //Phosphodiester group
        G4String Subcomponent6d = "Sugar";

        G4Orb* sSugar_48em1_nm = new G4Orb(Subcomponent6d, 0.48 * nanometer);

        G4ThreeVector posi(0.180248 * nanometer,
                           0.32422 * nanometer,
                           0.00784 * nanometer);
        G4UnionSolid* uniDNA = new G4UnionSolid("move",
                                                sSugar_48em1_nm,
                                                sSugar_48em1_nm,
                                                0,
                                                posi);

        G4ThreeVector posi2(-0.128248 * nanometer,
                            0.41227 * nanometer,
                            0.03584 * nanometer);

        G4UnionSolid* uniDNA2 = new G4UnionSolid("move2",
                                                 sSugar_48em1_nm,
                                                 sSugar_48em1_nm,
                                                 0,
                                                 posi2);

        //*******************************************************************
        // Phosphodiester group : Position
        //*******************************************************************

        for (G4int n = 2; n < 200; n++)
        {
            G4double SP1[2][3] = {
                { (-0.6 * nanometer) * cos(n * 0.26),
                    0, (0.6* nanometer) * sin(n * 0.26) },
                { (0.6 * nanometer) * cos(n * 0.26),
                    0, (-0.6 * nanometer) * sin(0.26 * n) }
            };
            G4double matriceSP1[3][3] = {
                { cos(n * 0.076), -sin(n * 0.076), 0 },
                { sin(n * 0.076), cos(n * 0.076), 0 },
                { 0, 0, 1 }
            };
            G4double matriceSP2[2][3];

            for (G4int i = 0; i < 3; i++)
            {
                G4double sumSP1 = 0;
                G4double sumSP2 = 0;
                for (G4int j = 0; j < 3; j++)
                {
                    sumSP1 += matriceSP1[i][j] * SP1[0][j];
                    sumSP2 += matriceSP1[i][j] * SP1[1][j];
                }
                matriceSP2[0][i] = sumSP1;
                matriceSP2[1][i] = sumSP2;
            }

            G4double heliceSP[3] = {
                (4.85 * nanometer) * cos(n * 0.076),
                (4.85 * nanometer) * sin(n * 0.076),
                (n * 0.026 * nanometer)
            };

            for (G4int i = 0; i < 3; i++)
            {
                matriceSP2[0][i] += heliceSP[i];
                matriceSP2[1][i] += heliceSP[i];
            }
            G4ThreeVector posSugar1(matriceSP2[0][2],
                                    matriceSP2[0][1],
                                    (matriceSP2[0][0]) - (4.25 * nanometer));
            G4ThreeVector posSugar2(matriceSP2[1][2],
                                    matriceSP2[1][1],
                                    (matriceSP2[1][0]) - (5.45 * nanometer));

            ostringstream ss;
            ss << "sugar_" << n;
            name = ss.str().c_str();
            ss.str("");
            ss.clear();

            //  snprintf(name, countof(name), "sugar %d", n);
            uniDNA = new G4UnionSolid(name,
                                      uniDNA,
                                      sSugar_48em1_nm,
                                      0,
                                      posSugar1);

            ss << "sugar_" << n;
            name = ss.str().c_str();
            ss.str("");
            ss.clear();

            //  snprintf(name, countof(name), "sugar %d", n);
            uniDNA2 = new G4UnionSolid(name,
                                       uniDNA2,
                                       sSugar_48em1_nm,
                                       0,
                                       posSugar2);
        }

        G4String SubComponent6e = "lSugar2";
        G4String SubComponent6f = "lSugar4";

        G4LogicalVolume* lSphere3 = CreateLogicalVolume(SubComponent6e, uniDNA);

        G4LogicalVolume* lSphere4 = CreateLogicalVolume(SubComponent6f, uniDNA2);

        /**************************************************************************
         Base pair Position
         **************************************************************************/
        for (G4int n = 0; n < 200; n++)
        {
            G4double bp1[2][3] = {
                { (-0.34 * nanometer) * cos(n * 0.26),
                    0, (0.34* nanometer) * sin(n * 0.26) },
                { (0.34 * nanometer) * cos(n * 0.26),
                    0, (-0.34 * nanometer) * sin(0.26 * n) }
            };
            G4double matriceBP1[3][3] = {
                { cos(n * 0.076), -sin(n * 0.076), 0 },
                {sin(n * 0.076), cos(n * 0.076), 0 },
                { 0, 0, 1 }
            };
            G4double matriceBP2[2][3];

            for (G4int i = 0; i < 3; i++)
            {
                G4double sumBP1 = 0;
                G4double sumBP2 = 0;
                for (G4int j = 0; j < 3; j++)
                {
                    sumBP1 += matriceBP1[i][j] * bp1[0][j];
                    sumBP2 += matriceBP1[i][j] * bp1[1][j];
                }
                matriceBP2[0][i] = sumBP1;
                matriceBP2[1][i] = sumBP2;
            }
            G4double heliceBP[3] = {
                (4.8 * nanometer) * cos(n * 0.076),
                (4.8 * nanometer) * sin(n * 0.076),
                n * 0.026 * nanometer
            };

            for (G4int i = 0; i < 3; i++)
            {
                matriceBP2[0][i] += heliceBP[i];
                matriceBP2[1][i] += heliceBP[i];
            }
            G4ThreeVector position1(matriceBP2[0][2],
                                    matriceBP2[0][1],
                                    matriceBP2[0][0] - (4.25 * nanometer));
            G4ThreeVector position2(matriceBP2[1][2],
                                    matriceBP2[1][1],
                                    matriceBP2[1][0] - (5.45 * nanometer));

            new G4PVPlacement(0,
                              position1,
                              lBp1,
                              "physi blue sphere",
                              lSphere3,
                              false,
                              0);
            new G4PVPlacement(0,
                              position2,
                              lBp2,
                              "physi pink sphere",
                              lSphere4,
                              false,
                              0);
        }

        /****************************************************************************/
        //                 Initial position of different elements
        /****************************************************************************/
        // DNA and histone positions
        for (int j = 0; j < 90; j++)
        {
            // DNA (bp-SP)
            G4RotationMatrix* rotStrand1 = new G4RotationMatrix;
            rotStrand1->rotateZ(j * -51.43 * degree);
            G4ThreeVector posStrand1(-2.7 * nanometer,
                                     9.35 * nanometer,
                                     (-69.9 * nanometer) + (j * 1.67 * nanometer));
            posStrand1.rotateZ(j * 51.43 * degree);
            new G4PVPlacement(rotStrand1,
                              posStrand1,
                              lSphere3,
                              "physi sugar 2",
                              lEnv,
                              false,
                              0);

            G4RotationMatrix* rotStrand2 = new G4RotationMatrix;
            rotStrand2->rotateZ(j * -51.43 * degree);
            G4ThreeVector posStrand2(-2.7 * nanometer,
                                     9.35 * nanometer,
                                     (-68.7 * nanometer) + (j * 1.67 * nanometer));
            posStrand2.rotateZ(j * 51.43 * degree);
            new G4PVPlacement(rotStrand2,
                              posStrand2,
                              lSphere4,
                              "physi sugar 4",
                              lEnv,
                              false,
                              0);

            // histones
            G4RotationMatrix* rotHistone = new G4RotationMatrix;
            rotHistone->rotateY(90 * degree);
            rotHistone->rotateX(j * (-51.43 * degree));
            G4ThreeVector posHistone(0.0,
                                     9.35 * nanometer,
                                     (-74.15 + j * 1.67) * nanometer);
            posHistone.rotateZ(j * 51.43 * degree);
            new G4PVPlacement(rotHistone,
                              posHistone,
                              lHistone,
                              "PV histone",
                              lEnv,
                              false,
                              0);
        }
    }

	InstantiateChildren(fEnvelopePhys);

	return fEnvelopePhys;
}



