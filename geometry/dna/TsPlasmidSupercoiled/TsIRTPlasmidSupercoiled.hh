//
// ********************************************************************
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * TOPAS collaboration.                                             *
// * Use or redistribution of this code is not permitted without the  *
// * explicit approval of the TOPAS collaboration.                    *
// * Contact: Joseph Perl, perl@slac.stanford.edu                     *
// *                                                                  *
// ********************************************************************
//

#ifndef TsIRTPlasmidSupercoiled_hh
#define TsIRTPlasmidSupercoiled_hh

#include "TsVGeometryComponent.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Orb.hh"
#include "G4Point3D.hh"

#include <map>
#include <vector>

class G4VSolid;
class G4Material;
class G4LogicalVolume;

class TsIRTPlasmidSupercoiled : public TsVGeometryComponent
{
public:
	TsIRTPlasmidSupercoiled(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
						 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	~TsIRTPlasmidSupercoiled();
	
	G4VPhysicalVolume* Construct();
	
public:
	G4LogicalVolume* CreateLogicVolume(G4String fileName, G4int copyNumber,
									   G4RotationMatrix* rot, G4ThreeVector* trans);
	
private:
	
	struct TempMolecule
	{
		TempMolecule(std::string name, int copyNumber, G4ThreeVector position, double radius, double waterRadius, std::string material, int strand)
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
		
		bool operator<(const TempMolecule& str) const
		{
			return (fPosition.z() < str.fPosition.z() );
		}
	};
	
	std::vector<std::pair<G4int, std::pair<G4int, G4ThreeVector>>> fDNASpecies;
	std::string fGeoName;
	
	std::map<std::string, G4double> fRadiusMap;
	std::map<std::string, G4double> fWaterRadiusMap;
	
	std::vector<TempMolecule> fMolecules;
	G4Orb* fSolid;
	
	void ParseFile(G4String fileName);
	G4VSolid* CreateCutSolid(G4Orb *solidOrbRef,
							 TempMolecule &molRef,
							 std::vector<TempMolecule> &molList,
							 G4bool in);
	G4double fSize;
	G4double fOffsetX;
	G4double fOffsetY;
	G4double fOffsetZ;
	G4double fXMin;
	G4double fXMax;
	G4double fYMin;
	G4double fYMax;
	G4double fZMin;
	G4double fZMax;
 
};

#endif

