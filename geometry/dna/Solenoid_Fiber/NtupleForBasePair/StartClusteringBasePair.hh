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
//  StartClustering.h
//  Fibre
//
//  Created by Nick Henthorn on 31/03/2017.

#ifndef StartClusteringBasePair_h
#define StartClusteringBasePair_h

#include "ClusterAlgoBasePair.hh"
#include "HitPoint.hh"

#include "TsVScorer.hh"

using namespace std;
class StartClusteringBasePair{
public:
    StartClusteringBasePair(TsParameterManager *fPm);
    ~StartClusteringBasePair(){}
    void Cluster(vector<HitPoint*> &Hits);

    //the following are public so "ReadData.cc" can use them
    void SplitList(vector<HitPoint*>&Backbones, vector<HitPoint*>&Bases, vector<HitPoint*>&Hits);
    void CheckHits(vector<HitPoint*>&Hits);

private:
    G4bool EnergyAccept(HitPoint*Hit);
    
    G4bool EneRange=true;
    G4bool EneThresh=false;
    G4bool EneIonis = false;

    G4double Emin=5.0*eV;
    G4double Emax=37.5*eV;
    G4double EThresh=17.5*eV;
};


#endif /* StartClustering_h */
