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
/*
 *
 *  Created on: 15 Feb 2016
 *      Author: Nick Henthorn
 */

#ifndef TsSolenoidFiber_hh
#define TsSolenoidFiber_hh

#include "TsVGeometryComponent.hh"
#include "G4NistManager.hh"

using namespace CLHEP;
using namespace std;


class G4VPhysicalVolume;
class G4Sphere;
class G4Cons;
class G4LogicalVolume;
class G4VPhysicalVolume;
class DetectorMessenger;
class G4Orb;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


class TsSolenoidFiber : public TsVGeometryComponent
{
public:
    TsSolenoidFiber(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
            TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
    ~TsSolenoidFiber();

    G4VPhysicalVolume* Construct();

private:

//    G4LogicalVolume* fEnvelopeLog;
//    G4VPhysicalVolume* fEnvelopePhys;

    void BuildHistones(vector<pair<G4ThreeVector, G4RotationMatrix*>> &HistoneDetails,
                       G4double ChromatinRadius,
                       G4double ChromatinLength);

    void BuildDNA(vector<pair<G4ThreeVector, G4RotationMatrix*>> &HistoneDetails);


    G4LogicalVolume *lBase1, *lBase2, *lBack1, *lBack2;


    void SetDNAVolumes(G4bool BuildHalfCyl,
                       G4bool BuildQuartCyl,
                       G4bool BuildSphere);
    void GenerateDNAPath(vector<pair<G4ThreeVector, G4RotationMatrix*>> &HistoneDetails,
                         vector<G4ThreeVector> &path);
    void SegmentDNAPath(vector<G4ThreeVector> &path);
    void PlaceDNASphere(vector<G4ThreeVector> &path);
    void PlaceDNA(vector<G4ThreeVector> &path);
    void ApplyRotation(G4ThreeVector &Rotated,
                       G4ThreeVector &Position,
                       G4RotationMatrix *Rot);
    void Bezier(G4ThreeVector &start,
                G4ThreeVector &MidPoint1,
                G4ThreeVector & MidPoint2,
                G4ThreeVector &end,
                vector<G4ThreeVector> &path,
                G4int nSteps);
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
