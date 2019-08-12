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

#ifndef GEOIMPORT_HH
#define GEOIMPORT_HH

#include <map>
#include <fstream>
#include <algorithm>

#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "G4Orb.hh"
#include "G4Ellipsoid.hh"
#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4SystemOfUnits.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"

// ************************************
// Useful structures
// ************************************

/*!
 * \brief The Molecule struct
 */
struct Molecule1
{
    Molecule1(std::string name, int copyNumber, G4ThreeVector position, double radius, double waterRadius, std::string material, int strand)
    {
        fName = name;
        fMaterial = material;
        fCopyNumber = copyNumber;
        fPosition = position;
        fRadius = radius;
        fRadiusWater = waterRadius;
        fStrand = strand;
    }

    std::string fName;
    std::string fMaterial;

    int fCopyNumber;
    int fStrand;

    G4ThreeVector fPosition;

    double fRadius;
    double fRadiusWater;

    // To sort the molecules in function of their z coordinate
    bool operator<(const Molecule1& str) const
    {
        return (fPosition.z() < str.fPosition.z() );
    }
};


/*!
 * \brief The PhysGeoImport class
 * This class deals with the import of a DnaFabric geometry into this Geant4 user application
 */
class PhysGeoImport
{
public:
    PhysGeoImport(); ///< Construcror
    PhysGeoImport(bool isVisu); ///< Constructor with argument
    ~PhysGeoImport(); ///< Destructor

    void SetFactor(double factor){fFactor=factor;} ///< Set the scale factor
    double GetFactor() const {return fFactor;} ///< Get the scale factor
    G4String GetGeoName() const {return fGeoName;} ///< Get the geometry name

    G4LogicalVolume* CreateLogicVolumeDNA(const std::string& fileName);


private:

    bool fIsVisu;
    double fFactor;
    double fSize;

    std::string fGeoName;
    std::string fNucleusName;
    std::string fNucleusType;
    std::map<std::string, double> fNucleusData;
    std::map<std::string, double> fRadiusMap;
    std::map<std::string, double> fWaterRadiusMap;

    // Vector to contain all the molecule structures listed within the imput file
    std::vector<Molecule1> fMolecules;

    // To check if this is the first voxel of the chromosome
    std::map<int, bool> fFirstMap;

    // Materials
    std::vector<G4Material*> fMaterialVect;
    G4Material* fpWater;

    G4Material* fVacuum;

    void ParseFile(const std::string& fileName);
    G4VSolid* CreateCutSolid(G4Orb *solidOrbRef,
                             Molecule1 &molRef,
                             std::vector<Molecule1> &molList,
                             G4bool in);
    void DefineMaterial();
};

#endif // GEOIMPORT_HH
