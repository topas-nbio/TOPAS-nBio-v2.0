// Scorer for NtupleForBasePair
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

#include "NtupleForBasePair.hh"
#include "G4Threading.hh"



NtupleForBasePair::NtupleForBasePair(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM,
                                     TsScoringManager* scM, TsExtensionManager* eM,
                                     G4String scorerName, G4String quantity, G4String outFileName,
                                     G4bool isSubScorer)
: TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{
    //Define ntuple columns

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
    if (fPm->ParameterExists("EventByEvent")){
        if(fPm->GetBooleanParameter("EventByEvent")){
            SingleTrack=true;
        } else {
            SingleTrack=false;
        }
    }

    //name the sensitive volumes
    if (fPm->ParameterExists("Base_strand1_name")){
    	Base_strand1_name=fPm->GetStringParameter("Base_strand1_name");
    }
    if (fPm->ParameterExists("Base_strand2_name")){
    	Base_strand2_name=fPm->GetStringParameter("Base_strand2_name");
    }
    if (fPm->ParameterExists("Back_strand1_name")){
    	Back_strand1_name=fPm->GetStringParameter("Back_strand1_name");
    }
    if (fPm->ParameterExists("Back_strand2_name")){
    	Back_strand2_name=fPm->GetStringParameter("Back_strand2_name");
    }

}


NtupleForBasePair::~NtupleForBasePair() {
}


G4bool NtupleForBasePair::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
    if (!fIsActive) {
        fSkippedWhileInactive++;
        return false;
    }


    ResolveSolid(aStep);

    G4StepPoint *theStepPoint = aStep->GetTrack()->GetStep()->GetPreStepPoint();
    G4TouchableHandle theTouchable = theStepPoint->GetTouchableHandle();

    fEneDep = aStep->GetTotalEnergyDeposit()/eV;

    G4String Volume = theTouchable->GetVolume()->GetName();

    if (fEneDep<=0.0){
      return false;
    }

    if (Volume.find(Base_strand1_name) == std::string::npos &&
        Volume.find(Base_strand2_name) == std::string::npos &&
        Volume.find(Back_strand1_name) == std::string::npos &&
        Volume.find(Back_strand2_name) == std::string::npos){
        return false;
    }

    //get strand and volume
    if (Volume.find(Base_strand1_name) != std::string::npos || Volume.find(Back_strand1_name) != std::string::npos){
        fStrand=1;
        if (Volume.find(Base_strand1_name) != std::string::npos){
            fIsBase=true;
            fIsBack=false;
        } else {
            fIsBase=false;
            fIsBack=true;
        }
    } else {
        fStrand=2;
        if (Volume.find(Base_strand2_name) != std::string::npos){
            fIsBase=true;
            fIsBack=false;
        } else {
            fIsBase=false;
            fIsBack=true;
        }
    }

    //get copy num (basepair ID)
    fBasePair=theTouchable->GetCopyNumber();


    G4String ProcessName = aStep->GetTrack()->GetStep()->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    G4String find = "G4DNAIonisation";
    if (ProcessName.find(find) != std::string::npos){
        fIsIonis=true;
    } else {
        fIsIonis=false;
    }

    if (ProcessName=="hIoni" || ProcessName=="eIoni"){
      fIsIonis=true;
    }

    //get position
    G4ThreeVector pos = theStepPoint->GetPosition()/nm;
    fPosX=pos.x();
    fPosY=pos.y();
    fPosZ=pos.z();

    //set Run id's
    fEventID=GetEventID();
    fRunID=GetRunID();

    //Put the units back to default
    fEneDep*=eV;
    fPosX*=nm;
    fPosY*=nm;
    fPosZ*=nm;

    //record hit
    fNtuple -> Fill();

    return true;
}


void NtupleForBasePair::UserHookForEndOfEvent()
{
}

void NtupleForBasePair::UserHookForEndOfRun()
{
    if (G4Threading::IsMasterThread()){

        string Filename = "OutputFile";
        //Filename+=".phsp";
        fNtuple->SetFileName(Filename);
        fNtuple->Write();

        //G4cout<<"ntuple size: "<<fNtuple->GetNumberOfEntries()<<G4endl;

        //Read in the Ntuple & set up hit store
        vector<HitPoint*> Hits;
        ReadOutput(Hits, Filename);

        //BasePair Clustering
        if (Hits.size()>0){
            //split the hits by event ID or run ID
            map<int, vector<HitPoint*>> EventHits;
            if (SingleTrack){
                for (size_t i=0;i<Hits.size();i++){
                    EventHits[Hits[i]->GetEventID()].push_back(Hits[i]);
                }
            } else {
                for (size_t i=0;i<Hits.size();i++){
                    EventHits[Hits[i]->GetRunID()].push_back(Hits[i]);
                }
            }

            for (auto data:EventHits){
                //start clustering (sum energy in volumes, and cluster)
                StartClusteringBasePair*clustering = new StartClusteringBasePair(fPm);
                clustering->Cluster(data.second);
                //delete clustering;
            }

            //write damage details
            WriteDamageSpecBasePair * DamSpec = new WriteDamageSpecBasePair(fPm);
            DamSpec->Write(Hits);
            //delete DamSpec;

            EventHits.clear();
        }


        //clear lists
        //for (size_t i=0;i<Hits.size();i++){
            //delete Hits[i];
            //Hits[i]=NULL;
        //}
        Hits.clear();

    } else {
        G4cout<<"Worker ntuple size: "<<fNtuple->GetNumberOfEntries()<<G4endl;
    }



}

void NtupleForBasePair::ReadOutput(vector<HitPoint*> &Hits, string filename)
{
    //string filename = fPm->GetStringParameter(GetFullParmName("OutputFile"))+".phsp";
    G4cout<<"Reading In "<<filename+".phsp"<<G4endl;
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
        G4int EventID = 0;  //eventID
        G4int RunID = 0;  //runID
        G4int strand = 0;  //strand
        G4bool IsIonis;  //is ionis
        G4int BP = -1;  //the base pair
        G4bool IsBack;
        G4bool IsBase;

        stream
        >>x
        >>y
        >>z
        >>Energy
        >>EventID
        >>RunID
        >>strand
        >>IsIonis
        >>BP
        >>IsBack
        >>IsBase;

        HitPoint * hit = new HitPoint();
        hit->SetEnergy(Energy);
        hit->SetPosition(G4ThreeVector(x,y,z));
        hit->SetStrand(strand);
        hit->SetIsIonis(IsIonis);
        hit->SetIsPhys(true);
        hit->SetEventID(EventID);
        hit->SetRunID(RunID);
        hit->SetBP(BP);
        hit->SetIsBack(IsBack);
        hit->SetIsBase(IsBase);

        Hits.push_back(hit);
    }

}
