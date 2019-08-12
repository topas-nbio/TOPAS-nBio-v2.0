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

#ifndef GEOMANAGER_HH
#define GEOMANAGER_HH

#include "GeoCalculation.hh"
#include "GeoVolume.hh"


class GeoManager
{
public:
    GeoManager(G4int verbose=0, G4double factor=1.);
    ~GeoManager();

    void Initialize();

    G4LogicalVolume* BuildLogicFiber(G4bool isVisu=false);

    //void RemoveMoleculeFromChemStage(const G4String& molName, G4int copyNum, G4int strand)
    //{geoDnaMolecules->RemoveMoleculeFromChemStage(molName, copyNum, strand);}
    //void ResetRemoveMoleculeList(){geoDnaMolecules->ResetRemoveMoleculeList();}

    // Getters
    G4int GetVerbose(){return fVerbose;}
    G4double GetFactor(){return fFactor;}
    std::map<G4String, std::vector<std::vector<double> > >* GetDNAMoleculesPositions();

    // Setters
    void SetVerbose(G4int verbose){fVerbose = verbose;}
    void SetFactor(G4double factor){fFactor = factor;}

private:
    G4int fVerbose;
    G4double fFactor;

    GeoCalculation* geoCalculation;
    GeoVolume* geoVolume;



};

#endif // GEOMANAGER_HH
