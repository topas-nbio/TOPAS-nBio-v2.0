// Extra Class for use by NtupleForSpatial
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
//  WriteDamageSpecSpatial.cpp
//  Clustering
//
//  Created by Nick Henthorn on 9/5/17.

#include <stdio.h>
#include "WriteDamageSpecSpatial.hh"

using namespace std;
using namespace CLHEP;


WriteDamageSpecSpatial::WriteDamageSpecSpatial(TsParameterManager *fPm)
{
    Target = fPm->GetStringParameter("SensitiveVolumeName");
    Option_WriteHeader = fPm->GetBooleanParameter("WriteHeader");
    Option_EvtByEvt = fPm->GetBooleanParameter("EventByEvent");

    G4String Particle = fPm->GetStringParameter("BeamParticle");
    BeamEnergySpread = fPm->GetUnitlessParameter("BeamEnergySpread");
    EnergyDistribution = fPm->GetStringParameter("BeamEnergyDistribution");

    if (BeamEnergySpread>0.0){EnergyDistribution="Gaussian";}
    BeamEnergy = fPm->GetDoubleParameter("BeamEnergy", "Energy");
    BeamStDev = BeamEnergy*(BeamEnergySpread/100.);

    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();;
    ParticlePDG=particleTable->FindParticle(Particle)->GetPDGEncoding();

    //set up Parameters for Damage Output
    shape = fPm->GetIntegerParameter("Shape");
    double x = fPm->GetDoubleParameter("HLX", "Length") /um;
    double y = fPm->GetDoubleParameter("HLY", "Length") /um;
    double z = fPm->GetDoubleParameter("HLZ", "Length") /um;

    TargetSize=G4ThreeVector(x,y,z);

    TopasVersion = fPm->GetTOPASVersion();

    if (fPm->ParameterExists("StoreDamagesInMemory")){
      StoreDamages = fPm->GetBooleanParameter("StoreDamagesInMemory");
    }

    if (fPm->ParameterExists("AddChromosomeDomains")){
      AddTerritories = fPm->GetBooleanParameter("AddChromosomeDomains");
    }

    if (fPm->ParameterExists("ProgramName")){
      ProgramName = fPm->GetStringParameter("ProgramName");
    }

    if (fPm->ParameterExists("WriteSSB")){
      WriteSSB = fPm->GetBooleanParameter("WriteSSB");
    }


}

WriteDamageSpecSpatial::~WriteDamageSpecSpatial()
{}


//start code with raw hits -- convert to cluster map
void WriteDamageSpecSpatial::Write(vector<HitPoint*> &Hits)
{
    string filename = "SDDFormatDamage_Spatial.txt";
    int nEvents=0;

    //read in previous data, if exists
    vector<string> PreviousHeader;
    vector<string> PreviousHits;
    bool PreviousData = CanReadIn(PreviousHeader,
                                  PreviousHits,
                                  nEvents,
                                  filename);

    map<int,vector<HitPoint*>> HitsMap;
    SortRawHits(Hits,
                HitsMap);


    //sort current data into clusters
    map<int, map<int, Clusters*>> ClusterMap;
    for (auto HitData:HitsMap){
        map<int, Clusters*> EvtClusters;

        for (size_t i=0;i<HitData.second.size();i++){
            Clusters*clust;
            int clusterID=HitData.second[i]->GetClusterID();
            if (EvtClusters.find(clusterID)!=EvtClusters.end()){
                clust=EvtClusters[clusterID];
                EvtClusters[clusterID]->AddPoint(HitData.second[i]);
            }else{
                clust=new Clusters();
                clust->SetClusterID(clusterID);
                clust->AddPoint(HitData.second[i]);
                EvtClusters[clusterID]=clust;
            }
        }
        ClusterMap[HitData.first] = EvtClusters;
    }

    //add chromosome domains
    if (AddTerritories){
      cout<<"DAMAGE MODEL: Adding Chromosome Domains"<<endl;
      AddChromosomeDomain*domain = new AddChromosomeDomain();
      domain->Start(ClusterMap, TargetSize.x());
      delete domain;
    }

    //add new events to event counter
    for (auto data:ClusterMap){
        for (auto d:data.second){
            //SSB or CompSSB
            if (d.second->GetClusterID()<=0){
                nEvents+=d.second->GetNbPts();
            }
            //DSB
            if (d.second->GetClusterID()>0){
                nEvents++;
            }
        }
    }

    ofstream outfile;
    outfile.open(filename);

    //if exists -- write previous data
    if (PreviousData){
        G4cout<<G4endl<<"Previous Data Exists: "<<filename<<G4endl;
        G4cout<<"Adding new data"<<G4endl<<G4endl;
        WritePreviousEvents(PreviousHeader,
                            PreviousHits,
                            nEvents,
                            outfile);
    }

    //else -- write new header
    else {
        //dont write header for distributed simulation
        if (Option_WriteHeader){
            WriteHeader(outfile,
                        nEvents);
        }

    }

    //write new events
    WriteEvents(ClusterMap,
                outfile);



}

//Sort raw hits into grouped map
//grouped by primaryID or RunID
void WriteDamageSpecSpatial::SortRawHits(vector<HitPoint *> &RawHits,
                                         map<int, vector<HitPoint *> > &HitsMap)
{
    //if exposure -- put all hits into one map element
    if (!Option_EvtByEvt){
        for (size_t i=0;i<RawHits.size();i++){
            HitsMap[RawHits[i]->GetRunID()].push_back(RawHits[i]);
        }
    }
    //if event-by-event -- put hits into map of primary ID
    else {
        for (size_t i=0;i<RawHits.size();i++){
            HitsMap[RawHits[i]->GetEventID()].push_back(RawHits[i]);
        }
    }
}

//read in previousData if exists
bool WriteDamageSpecSpatial::CanReadIn(vector<string> &PreviousHeader,
                                       vector<string> &PreviousHits,
                                       int &nEvents,
                                       string filename)
{
    ifstream infile(filename);
    if (!infile){
        return false;
    }

    string line, progName;
    bool readDat=false;
    bool readHeader=false;
    bool readEvents=false;
    bool readTargetName=false;

    while (getline(infile, line)){
        if (line=="# Header data" ||
            line=="# Header data\r"){
            readHeader=true;
        }

        if (readHeader){
            PreviousHeader.push_back(line);
        }

        //if different target to previous data, dont store previous
        if (readTargetName){
            bool Match=false;

            if (line==Target){Match = true;}
            else if (line+"\t"==Target){Match = true;}
            else if (line+"\r"==Target){Match = true;}
            else if (line+"\t\r"==Target){Match = true;}

            if (!Match){
                G4cout
                <<"Previous Target Doesn't Match: '"
                <<line<<"' : '"<<Target<<"'"
                <<endl;

                PreviousHeader.clear();
                PreviousHits.clear();
                return false;
            }
            readTargetName=false;
        }

        if (readEvents){
            int CurrEvents=0;
            stringstream stream(line.data());
            stream>>CurrEvents;
            nEvents+=CurrEvents;
            readEvents=false;
        }
        if (readDat){
            PreviousHits.push_back(line);
        }

        if (line == "# Line-by-line hit data listing" ||
            line == "# Line-by-line hit data listing\r"){
            readDat=true;
        }
        if (line == "# Event count" ||
            line == "# Event count\r"){
            readHeader=false;
            readEvents=true;
        }
        if (line == "# Target model" ||
            line == "# Target model\r"){
            readTargetName=true;
        }
    }

    return true;
}

//write out previous events
void WriteDamageSpecSpatial::WritePreviousEvents(vector<string> PreviousHeader,
                                                 vector<string> PreviousData,
                                                 int nEvents,
                                                 ofstream &outfile)
{
    for (size_t i=0;i<PreviousHeader.size();i++){
        outfile<<PreviousHeader[i]<<endl;
    }

    //how are events grouped
    string grouping;
    if (Option_EvtByEvt){
        grouping = "Single Event";
    } else {
        grouping = "Exposure";
    }

    outfile
    <<nEvents<<endl
    <<"# File event grouping"<<endl
    <<grouping<<endl
    <<"#"<<endl
    <<"# Line-by-line hit data listing"<<endl;

    for  (size_t i=0;i<PreviousData.size();i++){
        outfile<<PreviousData[i]<<endl;
    }
}

//write header
void WriteDamageSpecSpatial::WriteHeader(ofstream &outfile,
                                         int nEvents)
{
    //code name
    string ProgName = ProgramName+"_spatialDBSCAN_Topas"+TopasVersion;

    //target model
    string TargetName = Target;

    double density=0.0;
    if (shape == 0){
      density = (2.0*3087.0) / ( (4.0/3.0)*pi*TargetSize.x()*TargetSize.y()*TargetSize.z());
    }
    if (shape == 1){
      TargetSize*=2.0;
      density = (2.0*3087.0) / ( TargetSize.x()*TargetSize.y()*TargetSize.z() );
    }

    //particle ID
    int particleNum=ParticlePDG;

    //beam details
    string beamDist = EnergyDistribution;
    double beamMeanEne = BeamEnergy;
    double beamStDev = BeamStDev;
    string BeamDetails = beamDist+"\t"+to_string(beamMeanEne)+"\t"+to_string(beamStDev);
    if (beamDist=="Mono"){BeamDetails="";}

    //how are events grouped
    string grouping;
    if (Option_EvtByEvt){
        grouping = "Single Event";
    } else {
        grouping = "Exposure";
    }


    outfile
    <<"# Header data" <<endl
    <<ProgName<<endl
    <<"# Target model"<<endl
    <<TargetName<<endl
    <<"# Bounding volume"<<endl
    <<shape<<"\t"<<TargetSize.x()<<"\t"<<TargetSize.y()<<"\t"<<TargetSize.z()<<endl
    <<"# DNA Density"<<endl
    <<density<<endl  //the Mbp per volume
    <<"# O2 concentration"<<endl
    <<" "<<endl
    <<"# Incident particle"<<endl
    <<particleNum<<endl
    <<"# Mean Energy (MeV)"<<endl
    <<beamMeanEne<<endl
    <<"# Distribution function"<<endl
    <<BeamDetails<<endl
    <<"# Event count"<<endl
    <<nEvents<<endl
    <<"# File event grouping"<<endl
    <<grouping<<endl
    <<"#"<<endl
    <<"# Line-by-line hit data listing"
    <<endl;
}

void WriteDamageSpecSpatial::WriteEvents(map<int, map<int,Clusters*>> &ClusterMap,
                                         ofstream &outfile)
{
    int nDSBCount=0;
    for (auto data:ClusterMap){
        int NewDat=1;
        for (auto d:data.second){
            //SSBs
            if (d.second->GetClusterID()==0 && WriteSSB){
                vector<HitPoint*> SSBs = d.second->GetHits();
                int nSSB = SSBs.size();
                int Cause=-1;
                for (int i=0;i<nSSB;i++){
                    int nBack=0, nBase=0;

                    if (SSBs[i]->GetIsBase()){nBase=1;}
                    if (SSBs[i]->GetIsBack()){nBack=1;}

                    //int strand = SSBs[i]->GetStrand();
                    //int PosDelim=0;
                    //if (strand==1 && nBack==1){PosDelim=1;}
                    //if (strand==1 && nBase==1){PosDelim=2;}
                    //if (strand==2 && nBase==1){PosDelim=3;}
                    //if (strand==2 && nBack==1){PosDelim=4;}

                    if (SSBs[i]->GetIsPhys()){Cause=0;}
                    if (SSBs[i]->GetIsChem()){Cause=1;}
                    if (SSBs[i]->GetIsChem() && SSBs[i]->GetIsPhys()){Cause=2;}

                    G4int Chromosome=SSBs[i]->GetChromID();
                    G4int ChromosomeCopy=SSBs[i]->GetChromCopy();
                    G4int Chromatid=SSBs[i]->GetChromatid();
                    G4String C_Details="     ";
                    if (AddTerritories){
                      C_Details=to_string(Chromosome)+" "+to_string(ChromosomeCopy)+" "+to_string(Chromatid);
                    }

                    G4ThreeVector pos = SSBs[i]->GetPosition()/1000.0;
                    outfile
                    <<NewDat<<"\t"
                    <<pos.x()<<"\t"<<pos.y()<<"\t"<<pos.z()<<"\t"
                    <<nBase<<" "<<nBack<<" "<<0<<"\t"
                    <<Cause<<"\t"
                    <<C_Details<<"\t" //chrom ID 3 space separated fields
                    <<" "<<"\t"     //chromPos 1 field
                    <<" "<<"\t"     //full damage spec
                    <<" "           //full DNA sequence
                    <<endl;

                    if (StoreDamages){
                      //Add SSB to damage store
                      damage_event DE = damage_event();
                      DE.SetNewEvent(NewDat);
                      DE.SetPosition(pos);
                      DE.SetNumberOfBases(nBase);
                      DE.SetNumberOfBackbones(nBack);
                      DE.SetNumberOfDSB(0);
                      DE.SetCauseOfDamage(0);

                      if (AddTerritories){
                        DE.SetChromosomeNumber(Chromosome);
                        DE.SetChromosomeCopyNumber(ChromosomeCopy);
                        DE.SetChromatidNumber(Chromatid);
                      }

                      TsDamagePhaseSpaceStore::Instance()->AddEvent(DE);
                    }



                    NewDat=0;
                }

            }

            //SSB with base on opp strand within 10 bp
            else if (d.second->GetClusterID()<0 && WriteSSB){
                int nBack = d.second->GetnBacks();
                int nBase = d.second->GetnBases();

                int nHits=d.second->GetHits().size();
                bool isChem=false;
                bool isPhys=false;
                for (int i=0;i<nHits;i++){
                    if (d.second->GetHits()[i]->GetIsChem()){isChem=true;}
                    if (d.second->GetHits()[i]->GetIsPhys()){isPhys=true;}
                }

                int Cause=-1;
                if (isPhys){Cause=0;}
                if (isChem){Cause=1;}
                if (isPhys && isChem){Cause=2;}

                G4ThreeVector pos = d.second->GetCentre()/1000.0;

                G4int Chromosome=d.second->GetChromID();
                G4int ChromosomeCopy=d.second->GetChromCopy();
                G4int Chromatid=d.second->GetChromatid();

                cout<<Chromosome<<" "<<ChromosomeCopy<<" "<<Chromatid<<endl;

                G4String C_Details="     ";
                if (AddTerritories){
                  C_Details=to_string(Chromosome)+" "+to_string(ChromosomeCopy)+" "+to_string(Chromatid);
                }

                outfile
                <<NewDat<<"\t"
                <<pos.x()<<"\t"<<pos.y()<<"\t"<<pos.z()<<"\t"
                <<nBase<<" "<<nBack<<" "<<0<<"\t"
                <<Cause<<"\t"
                <<C_Details<<"\t" //chrom ID 3 space separated fields
                <<" "<<"\t"     //chromPos 1 field
                <<" "<<"\t"     //full damage spec
                <<" "           //full DNA sequence
                <<endl;

                if (StoreDamages){
                  //Add SSB to damage store
                  damage_event DE = damage_event();
                  DE.SetNewEvent(NewDat);
                  DE.SetPosition(pos);
                  DE.SetNumberOfBases(nBase);
                  DE.SetNumberOfBackbones(nBack);
                  DE.SetNumberOfDSB(0);
                  DE.SetCauseOfDamage(0);

                  if (AddTerritories){
                    DE.SetChromosomeNumber(Chromosome);
                    DE.SetChromosomeCopyNumber(ChromosomeCopy);
                    DE.SetChromatidNumber(Chromatid);
                  }

                  TsDamagePhaseSpaceStore::Instance()->AddEvent(DE);
                }


                NewDat=0;
            }

            //DSBs
            else if (d.second->GetClusterID()>0) {
                nDSBCount++;
                int nBack = d.second->GetnBacks();
                int nBase = d.second->GetnBases();

                int nHits=d.second->GetHits().size();
                bool isChem=false;
                bool isPhys=false;
                for (int i=0;i<nHits;i++){
                    if (d.second->GetHits()[i]->GetIsChem()){isChem=true;}
                    if (d.second->GetHits()[i]->GetIsPhys()){isPhys=true;}
                }

                int Cause=-1;
                if (isPhys){Cause=0;}
                if (isChem){Cause=1;}
                if (isPhys && isChem){Cause=2;}

                G4int Chromosome=d.second->GetChromID();
                G4int ChromosomeCopy=d.second->GetChromCopy();
                G4int Chromatid=d.second->GetChromatid();
                G4String C_Details="     ";
                if (AddTerritories){
                  C_Details=to_string(Chromosome)+" "+to_string(ChromosomeCopy)+" "+to_string(Chromatid);
                }


                G4ThreeVector pos = d.second->GetCentre()/1000.0;
                outfile
                <<NewDat<<"\t"
                <<pos.x()<<"\t"<<pos.y()<<"\t"<<pos.z()<<"\t"
                <<nBase<<" "<<nBack-2<<" "<<1<<"\t"
                <<Cause<<"\t"
                <<C_Details<<"\t" //chrom ID 3 space separated fields
                <<" "<<"\t"     //chromPos 1 field
                <<" "<<"\t"     //full damage spec
                <<" "           //full DNA sequence
                <<endl;


                if (StoreDamages){
                  //Add DSB to damage store
                  damage_event DE = damage_event();
                  DE.SetNewEvent(NewDat);
                  DE.SetPosition(pos);
                  DE.SetNumberOfBases(nBase);
                  DE.SetNumberOfBackbones(nBack-2);
                  DE.SetNumberOfDSB(1);
                  DE.SetCauseOfDamage(0);

                  if (AddTerritories){
                    DE.SetChromosomeNumber(Chromosome);
                    DE.SetChromosomeCopyNumber(ChromosomeCopy);
                    DE.SetChromatidNumber(Chromatid);
                  }

                  TsDamagePhaseSpaceStore::Instance()->AddEvent(DE);
                }



                NewDat=0;
            }
        }
    }

    outfile.close();


    G4cout<<G4endl;
    G4cout<<"DAMAGE MODEL: "<<nDSBCount<<" DSBs scored"<<G4endl;
    G4cout<<"DAMAGE MODEL: "<<TsDamagePhaseSpaceStore::Instance()->GetDamagePhaseSpace().size()<<" records in the store"<<G4endl;
    G4cout<<G4endl;

    ofstream outfile2;
    outfile2.open("DSB.dat", ios_base::app);
    outfile2<<nDSBCount<<endl;
    outfile2.close();

    ClearMap(ClusterMap);
}


//clear the data
void WriteDamageSpecSpatial::ClearMap(map<int, map<int, Clusters *> > &ClusterMap)
{
    for (auto data:ClusterMap){
        for (auto d:data.second){
            delete d.second;
            d.second=NULL;
        }
    }

    ClusterMap.clear();
}
