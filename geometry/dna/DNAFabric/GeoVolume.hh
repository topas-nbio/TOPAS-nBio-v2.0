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

#ifndef GEOVOLUME_HH
#define GEOVOLUME_HH

#include <map>
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VSolid.hh"
#include "G4Orb.hh"
#include "G4Box.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

struct DNAPlacementData;

class GeoVolume
{
public:
    GeoVolume(G4int verbose=0, G4double factor=1.);

    // Getters
    G4double GetSugarTHFRadiusWater(){return fSugarTHFRadiusWater;}
    G4double GetSugarTMPRadiusWater(){return fSugarTMPRadiusWater;}
    G4double GetSugarTHFRadius(){return fSugarTHFRadius;}
    G4double GetSugarTMPRadius(){return fSugarTMPRadius;}
    G4double GetBaseRadiusWater(){return fBaseRadiusWater;}
    G4double GetBaseRadius(){return fBaseRadius;}
    std::map<G4String, std::vector<std::vector<double> > >* GetDNAMoleculesPositions(){return fpDnaMoleculePositions;}

    // Setters
    void SetSugarTHFRadiusWater(G4double radius){fSugarTHFRadiusWater=radius;}
    void SetSugarTMPRadiusWater(G4double radius){fSugarTMPRadiusWater=radius;}
    void SetSugarTHFRadius(G4double radius){fSugarTHFRadius=radius;}
    void SetSugarTMPRadius(G4double radius){fSugarTMPRadius=radius;}
    void SetBaseRadiusWater(G4double radius){fBaseRadiusWater=radius;}
    void SetBaseRadius(G4double radius){fBaseRadius=radius;}
    void SetFiberPitch(G4double fiberPitch){fFiberPitch=fiberPitch;}
    void SetFiberNbNuclPerTurn(G4double nbNuclPerTurn){fFiberNbNuclPerTurn=nbNuclPerTurn;}
    void SetFiberDeltaAngle(G4double fiberDeltaAngle){fFiberDeltaAngle=fiberDeltaAngle;}
    void SetNucleoNum(G4double nucleoNum){fNucleoNum=nucleoNum;}
    void SetBpNum(G4double bpNum){fBpNum=bpNum;}
    void SetHistoneHeight(G4double histoneHeight){fHistoneHeight=histoneHeight;}
    void SetHistoneRadius(G4double histoneRadius){fHistoneRadius=histoneRadius;}

    G4LogicalVolume *BuildLogicFiber(std::vector<std::vector<DNAPlacementData> > *dnaVolPos,
                                     std::vector<G4ThreeVector> *posNucleo,
                                     std::map<G4ThreeVector, G4double> *posAndRadiusMap,
                                     G4bool isVisu=false);
    ~GeoVolume();
private:
    G4int fVerbose;
    G4double fFactor;
    G4Material* fWater;


    // Elementary parameters
    G4double fSugarTHFRadiusWater;
    G4double fSugarTMPRadiusWater;
    G4double fSugarTHFRadius;
    G4double fSugarTMPRadius;
    G4double fBaseRadiusWater;
    G4double fBaseRadius;

    G4double fFiberPitch;
    G4double fFiberNbNuclPerTurn;
    G4double fFiberDeltaAngle;

    G4double fHistoneHeight;
    G4double fHistoneRadius;

    G4int fNucleoNum;
    G4int fBpNum;

    // moleculeName: x, y, z, copyNumber, strand
    std::map<G4String, std::vector<std::vector<double> > >* fpDnaMoleculePositions;

    // Cut algorithm to avoid overlaps.
    // Idea: we must have a reference and a target. The reference is the solid we are considering and which could be cut if an overlap is
    // detected with the target solid. In a geometry, it implies we have to go through all the target solids for each reference solid.
    // This method will return the cutted spherical reference solid.
    G4VSolid *CreateCutSolid(G4Orb *solidOrbRef,
                             G4ThreeVector& posRef,
                             std::map<G4ThreeVector, G4double> *tarMap,
                             G4String volName = "",
                             G4bool in = false);

    // Create the solid and logical volumes required to build DNA around one histone.
    // Return a map as:
    // Key: name of the volume (backbone1, backbone2, base1, base2)
    // Content: vector of corresponding logical volumes
    std::map<G4String, std::vector<G4LogicalVolume *> >* CreateNucleosomeCuttedSolidsAndLogicals(std::vector<DNAPlacementData> *nucleosomeVolumePositions,
                                                                                                 std::map<G4ThreeVector, G4double> *posAndRadiusMap,
                                                                                                 G4bool isVisu=false);

    void CalculateMeanVol(std::map<G4String, std::vector<G4LogicalVolume *> > *logicSolidsMap);

};

#endif // GEOVOLUME_HH
