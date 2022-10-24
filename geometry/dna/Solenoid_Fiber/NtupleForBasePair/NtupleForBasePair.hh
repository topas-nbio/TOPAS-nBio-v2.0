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

#ifndef NtupleForBasePair_hh
#define NtupleForBasePair_hh

#include "TsVNtupleScorer.hh"

#include "G4VProcess.hh"

#include "HitPoint.hh"
#include "StartClusteringBasePair.hh"
#include "WriteDamageSpecBasePair.hh"


#include "G4PSDirectionFlag.hh"

#include "TsVGeometryComponent.hh"
// #include "TsGeometryManager.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UIcommand.hh"
#include "G4Tokenizer.hh"
#include "G4SystemOfUnits.hh"

#include "G4VPhysicalVolume.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"

class NtupleForBasePair : public TsVNtupleScorer
{
public:
    NtupleForBasePair(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM,
                    TsExtensionManager* eM, G4String scorerName, G4String quantity, G4String outFileName,
                    G4bool isSubScorer);

    virtual ~NtupleForBasePair();

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
    //G4double EMin=5.0;
    //G4double EMax = 37.5;
    //G4double EThresh = 17.5;

    //G4double SensitiveFraction=0.15;

    void ReadOutput(vector<HitPoint*> &Hits, string Filename);

    //G4bool Fiber=false;
    //G4bool Cell=false;

    G4bool SingleTrack=true;

    G4String Back_strand1_name="Backbone1";
    G4String Back_strand2_name="Backbone2";
    G4String Base_strand1_name="Base1";
    G4String Base_strand2_name="Base2";

};
#endif
