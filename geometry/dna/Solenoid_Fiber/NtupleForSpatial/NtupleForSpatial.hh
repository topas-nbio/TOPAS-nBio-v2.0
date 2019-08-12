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
 *  Created on: 19 Sep 2017
 *      Author: Nick Henthorn
 */

#ifndef NtupleForSpatial_hh
#define NtupleForSpatial_hh

#include "TsVNtupleScorer.hh"

#include "G4VProcess.hh"

#include "HitPoint.hh"
#include "ClusterAlgoSpatial.hh"
#include "WriteDamageSpecSpatial.hh"


#include "G4PSDirectionFlag.hh"

#include "TsVGeometryComponent.hh"
//#include "TsGeometryManager.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UIcommand.hh"
#include "G4Tokenizer.hh"
#include "G4SystemOfUnits.hh"

#include "G4VPhysicalVolume.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"

class NtupleForSpatial : public TsVNtupleScorer
{
public:
    NtupleForSpatial(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM,
                    TsExtensionManager* eM, G4String scorerName, G4String quantity, G4String outFileName,
                    G4bool isSubScorer);

    virtual ~NtupleForSpatial();

    G4bool ProcessHits(G4Step*,G4TouchableHistory*);

    void UserHookForEndOfEvent();

    void UserHookForEndOfRun();

protected:
    // Output variables
    G4double fEneDep=-1.0;
    G4double fPosX=-1000.0;
    G4double fPosY=-1000.0;
    G4double fPosZ=-1000.0;

    G4int fEventID=-1;
    G4int fRunID=-1;
    G4int fStrand=-1;
    G4int fBasePair=-1;

    G4bool fIsIonis=false;
    G4bool fIsBack=false;
    G4bool fIsBase=false;


private:
    G4double EMin=5.0;
    G4double EMax = 37.5;
    G4double EThresh = 17.5;

    G4bool UseErange = true;
    G4bool UseEthresh = false;
    G4bool UseIonisation = false;

    G4double SensitiveFraction=0.15;

    void ReadOutput(vector<HitPoint*> &Hits, string Filename);

    //G4bool Fibre=false;
    //G4bool Cell=false;

    G4bool SingleTrack=false;

    G4String SensitiveVolumeName = "Nucleus";
};
#endif
