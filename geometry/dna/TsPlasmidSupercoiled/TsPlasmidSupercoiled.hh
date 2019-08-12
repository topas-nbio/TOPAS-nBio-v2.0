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
 *  Developed by Nicholas Henthorn, PRECISE, University of Manchester
 *  nicholas.henthorn@manchester.ac.uk
 *  https://tinyurl.com/y7xuxw9d
 *  17/10/2018, LINK TO PUBLICATION
 */

#ifndef TsPlasmidSupercoiled_hh
#define TsPlasmidSupercoiled_hh

#include "TsVGeometryComponent.hh"
#include "G4NistManager.hh"
#include "DNACoordinates.hh"
#include "DNA.hh"

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


class TsPlasmidSupercoiled : public TsVGeometryComponent
{
public:
    TsPlasmidSupercoiled(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
            TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
    ~TsPlasmidSupercoiled();

    G4VPhysicalVolume* Construct();

private:

    
    void PlaceDNA(vector<DNA*> &DNApt,
                  G4bool BuildHalfCyl,
                  G4bool BuildQuartCyl,
                  G4bool BuildSphere);
    
    void Bezier(G4ThreeVector &start,
                G4ThreeVector &MidPoint1,
                G4ThreeVector &MidPoint2,
                G4ThreeVector &end,
                vector<G4ThreeVector> &path,
                G4int nSteps);
    
    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
