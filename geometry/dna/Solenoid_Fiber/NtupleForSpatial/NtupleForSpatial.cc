// Scorer for NtupleForSpatial
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

#include "NtupleForSpatial.hh"
#include "G4Threading.hh"



NtupleForSpatial::NtupleForSpatial(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM,
                                 TsScoringManager* scM, TsExtensionManager* eM,
                                 G4String scorerName, G4String quantity, G4String outFileName,
                                 G4bool isSubScorer)
: TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{
    //Define ntuple columns

G4cout<<"BUILDING SCORER"<<G4endl;
    fNtuple->RegisterColumnD(&fPosX, "Position X", "nm");
    fNtuple->RegisterColumnD(&fPosY, "Position Y", "nm");
    fNtuple->RegisterColumnD(&fPosZ, "Position Z", "nm");
    fNtuple->RegisterColumnD(&fEneDep, "Energy", "eV");
    fNtuple->RegisterColumnI(&fEventID, "Event ID");
    fNtuple->RegisterColumnI(&fRunID, "Run ID");
    fNtuple->RegisterColumnI(&fStrand, "Strand Num");
    fNtuple->RegisterColumnB(&fIsIonis, "IsIonisation");
    fNtuple->RegisterColumnI(&fBasePair, "BP");
    fNtuple->RegisterColumnB(&fIsBack, "IsBack");
    fNtuple->RegisterColumnB(&fIsBase, "IsBase");



//Set up some run parameters

    //The sensitive volume name
    if (fPm->ParameterExists("SensitiveVolumeName")){
      SensitiveVolumeName = fPm->GetStringParameter("SensitiveVolumeName");
    }

    //How much of the volume is Sensitive
    if (fPm->ParameterExists("SensitiveFraction")){
      SensitiveFraction = fPm->GetUnitlessParameter("SensitiveFraction");
    }

    //The track grouping
    if (fPm->ParameterExists("EventByEvent")){
      if(fPm->GetBooleanParameter("EventByEvent")){
        SingleTrack=true;
      } else {
        SingleTrack=false;
      }
    }

    //The damage mode
    if (fPm->ParameterExists("DamageMode")){
      G4String DMode = fPm->GetStringParameter("DamageMode");
      if (DMode == "Threshold"){
        UseEthresh=true;
      } else if (DMode == "Ionisation"){
        UseIonisation=true;
      } else {
        UseErange=true;
      }
    } else {
      UseErange=true;
    }

    //Set the energy values
    if (fPm->ParameterExists("MinEneThresh")){
      EThresh = fPm -> GetDoubleParameter("MinEneThresh", "Energy") /eV;
    }
    if (fPm->ParameterExists("MinEneRange")){
      EMin = fPm -> GetDoubleParameter("MinEneRange", "Energy") /eV;
    }
    if (fPm->ParameterExists("MaxEneRange")){
      EMax = fPm -> GetDoubleParameter("MaxEneRange", "Energy") /eV;
    }

}


NtupleForSpatial::~NtupleForSpatial() {}


G4bool NtupleForSpatial::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
    if (!fIsActive) {
        fSkippedWhileInactive++;
        return false;
    }

    ResolveSolid(aStep);

    fEneDep = aStep->GetTotalEnergyDeposit()/eV;

    G4String Volume = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();

    // if (Volume !="World"){
    //   G4cout<<Volume<<" "<<fEneDep<<G4endl;
    // }

    if (fEneDep>0.0){
        if (Volume!=SensitiveVolumeName){
          return false;
        }

        if (G4UniformRand()>SensitiveFraction){
          return false;
        }

        //if using erange
        if (UseErange){
          if (fEneDep<EMin){return false;}
          G4double proba = (fEneDep - EMin)/ (EMax-EMin);
          G4double rando = G4UniformRand();
          if (rando > proba){return false;}
        }

        //if using ethresh
        if (UseEthresh){
          if (fEneDep<EThresh){return false;}
        }

        //if using Ionisation
        G4String ProcessName = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
        G4String find = "G4DNAIonisation";
        if (ProcessName.find(find) != std::string::npos){
            fIsIonis=true;
        } else {
            fIsIonis=false;
        }
        if (ProcessName=="hIoni" || ProcessName=="eIoni"){
          fIsIonis=true;
        }
        if (UseIonisation && !fIsIonis){return false;}


        if (G4UniformRand()>0.5){
            fStrand=1;
        } else {
            fStrand=2;
        }

        fIsBase=false;
        fIsBack=true;
        fBasePair=-1;


        //get position
        G4ThreeVector pos = aStep->GetPreStepPoint()->GetPosition()/nm;
        fPosX=pos.x();
        fPosY=pos.y();
        fPosZ=pos.z();

        //set Run id's
        fEventID=GetEventID();
        fRunID=GetRunID();
        // ---------! /Generic Scoring !---------//

        //Put the units back to default
        fEneDep*=eV;
        fPosX*=nm;
        fPosY*=nm;
        fPosZ*=nm;

        //record hit
        fNtuple -> Fill();

        return true;

    }
    return false;
}


void NtupleForSpatial::UserHookForEndOfEvent()
{
}

void NtupleForSpatial::UserHookForEndOfRun()
{
  if (G4Threading::IsMasterThread()){

     string Filename = "OutputFile";
    // string name = Filename+".phsp";
    remove("OutputFile.phsp");
    //Filename+=".phsp";
    fNtuple->SetFileName(Filename);
    fNtuple->Write();

     //G4cout<<"ntuple size: "<<fNtuple->GetNumberOfEntries()<<G4endl;

     //Read in the Ntuple & set up hit store
     vector<HitPoint*> Hits;
     ReadOutput(Hits, Filename);

    //G4cout<<"New Run "<<Hits.size()<<G4endl;

     if (Hits.size()>0){
       //split the hits by event ID (if run on exposure put all hits into same event id)
       map<int, vector<HitPoint*>> EventHits;
       if (SingleTrack){
         int nHits = Hits.size();
         for (int i=0;i<nHits;i++){
           EventHits[Hits[i]->GetEventID()].push_back(Hits[i]);
         }
       } else {
         int nHits = Hits.size();
         for (int i=0;i<nHits;i++){
           EventHits[Hits[i]->GetRunID()].push_back(Hits[i]);
         }
       }

       for (auto data:EventHits){
         //run cluster algorithm
         ClusterAlgoSpatial *clustering = new ClusterAlgoSpatial();
         clustering->RunDBSCAN(data.second);
         delete clustering;
       }

       EventHits.clear();

        //write damage details -- splits data by evtID or runID depending on "SingleTrack"
        WriteDamageSpecSpatial * DamSpec = new WriteDamageSpecSpatial(fPm);
         DamSpec->Write(Hits);
         delete DamSpec;
       }

     //clear lists
     for (size_t i=0;i<Hits.size();i++){
         delete Hits[i];
         Hits[i]=NULL;
     }
     Hits.clear();


  } else {
    G4cout<<"Worker ntuple size: "<<fNtuple->GetNumberOfEntries()<<G4endl;
  }



}

void NtupleForSpatial::ReadOutput(vector<HitPoint*> &Hits, string filename)
{
  //string filename = fPm->GetStringParameter(GetFullParmName("OutputFile"))+".phsp";
  G4cout<<G4endl;
  G4cout<<"DAMAGE MODEL: Reading "<<filename+".phsp"<<G4endl;
  G4cout<<G4endl;
  ifstream infile;
  infile.open(filename+".phsp");
  string line;

  while (getline(infile, line)){
    string field;
    stringstream stream(line.data());

    G4double x = 0.;  //x (nm)
    G4double y = 0.;  //y (nm)
    G4double z = 0.;  //z (nm)
    G4double Energy = 0.;  //edep (eV)
    G4int PrimaryID = 0;  //eventID
    G4int RunID = 0;  //runID
    G4int strand = 0;  //strand
    G4bool IsIonis;  //is ionis
    G4int BP = -1;  //the base pair
    G4bool IsBack;
    G4bool IsBase;

    stream>>x>>y>>z>>Energy>>PrimaryID>>RunID>>strand>>IsIonis>>BP>>IsBack>>IsBase;

    HitPoint * hit = new HitPoint();
    hit->SetEnergy(Energy);
    hit->SetPosition(G4ThreeVector(x,y,z));
    hit->SetStrand(strand);
    hit->SetIsIonis(IsIonis);
    hit->SetIsPhys(true);
    hit->SetEventID(PrimaryID);
    hit->SetRunID(RunID);
    hit->SetBP(BP);
    hit->SetIsBack(IsBack);
    hit->SetIsBase(IsBase);

    Hits.push_back(hit);
  }

}
