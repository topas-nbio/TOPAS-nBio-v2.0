/*
*
*  Ingram & Henthorn et al
*  Hi-C TOPAS Geometry
*
*/

#ifndef TsHiC_hh
#define TsHiC_hh

#include "ChromosomeParser.hh"

#include "TsVGeometryComponent.hh"

#include "G4VisAttributes.hh"
#include "TsParameterManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4Box.hh"
#include "G4Ellipsoid.hh"
#include "G4Orb.hh"

class TsHiC : public TsVGeometryComponent
{
public:
    TsHiC(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
            TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
    ~TsHiC();

    G4VPhysicalVolume* Construct();

private:

    
    G4double transX=0.0, transY=0.0, transZ=0.0;
    G4double ExtraCytoplasmSize=2.5*um;
    void ResolveParameters();
    
    
    G4VPhysicalVolume* BuildCytoplasm(G4ThreeVector Dimensions,
                                    G4double ExtraCytoplasm,
                                    G4RotationMatrix* rot);
    G4VPhysicalVolume* BuildNucleus(G4ThreeVector Dimensions,
                                    G4VPhysicalVolume* pMother);
    void BuildChromosomes(G4VPhysicalVolume* pMother,
                          std::map<G4int,std::vector<ChromObj>> Beads);
    
    //Chromosome Colours
    std::map<G4int,G4Color> VolumeColor;

};

#endif
