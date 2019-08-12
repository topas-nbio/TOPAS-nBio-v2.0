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
//  WriteDamageSpecSpatial.h
//  Clustering
//
//  Created by Nick Henthorn on 9/5/17.

#ifndef WriteDamageSpecSpatial_h
#define WriteDamageSpecSpatial_h

#include "HitPoint.hh"
#include "Clusters.hh"
#include <map>

#include "TsVScorer.hh"
#include "G4ParticleTable.hh"

#include "TsDamagePhaseSpaceStore.hh"
#include "AddChromosomeDomain.hh"


class WriteDamageSpecSpatial{
public:
    WriteDamageSpecSpatial(TsParameterManager *fPm);
    ~WriteDamageSpecSpatial();

    void Write(vector<HitPoint*> &Hits);

private:
    void SortRawHits(vector<HitPoint*> &RawHits,
                     map<int,vector<HitPoint*>> &HitsMap);

    bool CanReadIn(vector<string> &PreviousHeader,
                   vector<string> &PreviousHists,
                   int &nEvents,
                   string filename);

    void WritePreviousEvents(vector<string> PreviousHeader,
                             vector<string> PreviousData,
                             int nEvents,
                             ofstream &outfile);

    void WriteHeader(ofstream &outfile,
                     int nEvents);

    void WriteEvents(map<int, map<int,Clusters*>> &ClusterMap,
                     ofstream &outfile);

    void ClearMap(map<int, map<int, Clusters*>>  &ClusterMap);

    string Target;
    bool Option_WriteHeader=true;
    bool Option_EvtByEvt=true;

    G4int shape;
    G4ThreeVector TargetSize;
    G4String TopasVersion;
    G4String ProgramName;
    G4int ParticlePDG;
    G4double BeamEnergySpread;
    G4String EnergyDistribution;
    G4double BeamEnergy;
    G4double BeamStDev;
    G4bool StoreDamages=false;
    G4bool AddTerritories=false;
    G4bool WriteSSB=false;



};

#endif /* WriteDamageSpecSpatial_h */










void WriteHeader(ofstream &outfile, map<int, map<int, Clusters*>>  &ClusterMap);
