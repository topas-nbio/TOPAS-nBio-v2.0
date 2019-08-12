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

#ifndef GEOCALCULATION_HH
#define GEOCALCULATION_HH

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include <map>

struct DNAPlacementData;

// vector contains vector which contains DNA volume positions for one base pair
// Schema:
//
//     Vector                SubVector containing DNAPlacementData
//     |
//     |nucleosome 1 --------| Position informations about base pair 1 (4 volumes)
//     |                     | Position informations about base pair 2 (4 volumes)
//     |                     | Position informations about base pair 3 (4 volumes)
//     |                     | etc.
//     |
//     |nucleosome 2 --------| Position informations about base pair 1 (4 volumes)
//     |                     | Position informations about base pair 2 (4 volumes)
//     |                     | Position informations about base pair 3 (4 volumes)
//     |                     | etc.
//
typedef std::vector<std::vector<DNAPlacementData> > DNAPosData;

class GeoCalculation
{
public:
    GeoCalculation(G4int verbose=0, G4double factor=1.);

    ~GeoCalculation();

    // Initialize the GeoCalculation class by setting up all the geometrical parameters.
    void Initialize();

    // Getters
    std::vector<G4ThreeVector>* GetNucleosomePosition(){return fPosNucleo;}
    DNAPosData* GetAllDNAVolumePositions(){return fPosDNA;}
    std::vector<DNAPlacementData>* GetDNAVolumePositionsForNucleosome(G4int nucl){return &fPosDNA->at(nucl);}
    std::map<G4ThreeVector, G4double>* GetPosAndRadiusMap(){return fPosAndRadiusMap;}
    G4double GetSugarTHFRadiusWater(){return fSugarTHFRadiusWater;}
    G4double GetSugarTMPRadiusWater(){return fSugarTMPRadiusWater;}
    G4double GetSugarTHFRadius(){return fSugarTHFRadius;}
    G4double GetSugarTMPRadius(){return fSugarTMPRadius;}
    G4double GetBaseRadiusWater(){return fBaseRadiusWater;}
    G4double GetBaseRadius(){return fBaseRadius;}
    G4double GetFiberPitch(){return fFiberPitch;}
    G4double GetFiberDeltaAngle(){return fFiberDeltaAngle;}
    G4double GetFiberNbNuclPerTurn(){return fFiberNbNuclPerTurn;}
    G4double GetBpNum(){return fBpNum;}
    G4double GetHistoneHeight(){return fHistoneHeight;}
    G4double GetHistoneRadius(){return fHistoneRadius;}

    // Setters
    void SetSugarTHFRadiusWater(G4double radius){fSugarTHFRadiusWater=radius;}
    void SetSugarTMPRadiusWater(G4double radius){fSugarTMPRadiusWater=radius;}
    void SetSugarTHFRadius(G4double radius){fSugarTHFRadius=radius;}
    void SetSugarTMPRadius(G4double radius){fSugarTMPRadius=radius;}
    void SetBaseRadiusWater(G4double radius){fBaseRadiusWater=radius;}
    void SetBaseRadius(G4double radius){fBaseRadius=radius;}

private:
    G4int fVerbose;
    G4double fFactor;

    // Elementary parameters
    G4double fSugarTHFRadius;
    G4double fSugarTMPRadius;
    G4double fBaseRadius;
    G4double fSugarTHFRadiusWater;
    G4double fSugarTMPRadiusWater;
    G4double fBaseRadiusWater;

    G4ThreeVector fPosSugarTMP1;
    G4ThreeVector fPosSugarTHF1;
    G4ThreeVector fPosBase1;
    G4ThreeVector fPosBase2;
    G4ThreeVector fPosSugarTHF2;
    G4ThreeVector fPosSugarTMP2;

    // Important numbers

    G4int fNucleoNum;
    G4int fBpNum;

    //***************************************
    // Calculation parameters
    //***************************************

    // Fiber parameters
    G4int fHistoneNum;
    G4double fHistoneRadius;
    G4double fHistoneHeight;
    G4double fFiberPitch;
    G4double fFiberCentralRadius;
    G4double fFiberNbNuclPerTurn;
    G4double fFiberDeltaAngle;
    std::vector<std::vector<G4double> > fFiberHelixMatVect;
    std::vector<std::vector<std::vector<G4double> > > fRotFiberMatVect;

    // DNA around histone parameters
    //
    // First helix parameters
    G4int bpNumAroundHistone;//const G4int bpNumAroundHistone = 154;
    G4double angleBpAroundHistone;

    // Second helix parameters (big simple helix)
    G4double secondHelixPitch;
    G4double centralRadius;
    G4double nbBasePairPerTurn;
    G4double deltaAngle ;

    // DNA linker parameters
    //
    G4int bpNumForLinker;//const G4int bpNumForLinker = 46;
    G4double linkerCentralRadius;
    G4double linkerHeightPerBp;
    G4double linkerArcCircleMat[3];
    G4double nbBasePairPerTurnForLinker;
    G4double deltaLinkerAngle;

    // Pos containers
    std::vector<G4ThreeVector>* fPosNucleo;
    DNAPosData* fPosDNA;
    std::map<G4ThreeVector, G4double>* fPosAndRadiusMap;

    //***************************************
    // Methods
    //***************************************

    // Calculate the different DNA volume positions based on the histone number given as arguement.
    // Positions are saved in a DNAPosData structure.
    // See beginning of the .hh file to get output structure informations (DNAPosData typedef).
    DNAPosData* CalculateDNAPosition(G4int histoneNum, G4ThreeVector& posSugarTMP1, G4ThreeVector &posSugarTHF1, G4ThreeVector &posBase1,
                                     G4ThreeVector &posBase2, G4ThreeVector &posSugarTHF2, G4ThreeVector &posSugarTMP2);

    // Calculate the nucleosome positions.
    // Useful to create the geometry in Geant4: we can create only once the volumes to make a nucleosome and then move it to the different
    // positions calculated here.
    std::vector<G4ThreeVector>* CalculateNucleosomePosition(G4int nucleoNum);

    // Use the DNA positons previously calculated to generate a map (map[key]=content).
    // Key: coordinates of the volume
    // Content: radius of the volume
    std::map<G4ThreeVector, G4double>* GenerateCoordAndRadiusMap(DNAPosData *dnaPosData);
    G4double GetAngleToXAxis(G4ThreeVector t);
};

#endif // GEOCALCULATION_HH
