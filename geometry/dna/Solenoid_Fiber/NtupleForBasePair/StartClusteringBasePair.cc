// Extra Class for use by NtupleForBasePair
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
//  StartClustering.cpp
//  Fibre
//
//  Created by Nick Henthorn on 31/03/2017.

#include <stdio.h>
#include "StartClusteringBasePair.hh"

StartClusteringBasePair::StartClusteringBasePair(TsParameterManager *fPm)
{
    //set up damage induction model
    if (fPm->ParameterExists(fPm->GetStringParameter("DamageMode"))){
        string DamageMode = fPm->GetStringParameter("DamageMode");
        if (DamageMode=="Ionisation"){
            EneRange=false;
            EneThresh=false;
            EneIonis = true;
        }
        if (DamageMode=="Range"){
            EneRange=true;
            EneThresh=false;
            EneIonis = false;
        }
        if (DamageMode=="Threshold"){
            EneRange=false;
            EneThresh=true;
            EneIonis = false;
        }
        
    }

    
    //set up energy values
    if (fPm->ParameterExists(fPm->GetDoubleParameter("MinEneRange", "Energy"))){
        Emin = fPm->GetDoubleParameter("MinEneRange", "Energy");
    }
    if (fPm->ParameterExists(fPm->GetDoubleParameter("MaxEneRange", "Energy"))){
        Emax = fPm->GetDoubleParameter("MaxEneRange", "Energy");
    }
    if (fPm->ParameterExists(fPm->GetDoubleParameter("MinEneThresh", "Energy"))){
        EThresh = fPm->GetDoubleParameter("MinEneThresh", "Energy");
    }
    
}

void StartClusteringBasePair::Cluster(vector<HitPoint*> &Hits){
    int HitsSize = Hits.size();

    if (HitsSize>0){
        //Check for multiple hits
        CheckHits(Hits);

        if (Hits.size()>0){
            //Split the hit list to backbone and bases
            vector<HitPoint*> BackboneHits;
            vector<HitPoint*> BaseHits;
            SplitList(BackboneHits, BaseHits, Hits);

            if (BackboneHits.size()>0){
                ClusterAlgoBasePair * Clust = new ClusterAlgoBasePair();
                Clust->RunDBSCAN(BackboneHits);
                Clust->ClusterBases(BaseHits);

                Clust->PotentialClusters(BackboneHits, BaseHits);


                delete Clust;
                Hits.clear(); //dont delete, ptrs sent to local analysis and deleted later

            }
        }
    }


}

//if multiple hits in the same volume, sum
void StartClusteringBasePair::CheckHits(vector<HitPoint*>&Hits)
{
    G4int nHits = Hits.size();
    vector<HitPoint*> HitsToDelete;
    G4int todelete=0;

    for (G4int i=0;i<nHits;i++){
        HitPoint*hit1=Hits[i];
        bool toMerge=false;
        vector<HitPoint*> MergeList;
        if (!hit1->GetChecked()){
            for (G4int j=i+1;j<nHits;j++){
                if (i!=j){
                    HitPoint*hit2=Hits[j];
                    if (!hit2->GetChecked()){
                        bool bothBases=false, bothBacks=false;
                        if (hit1->GetIsBack() && hit2->GetIsBack()){bothBacks=true;}
                        if (hit1->GetIsBase() && hit2->GetIsBase()){bothBases=true;}

                        if (bothBacks||bothBases){
                            bool sameStrand=false, sameBP=false;
                            if (hit1->GetStrand()==hit2->GetStrand()){sameStrand=true;}
                            if (hit1->GetBP()==hit2->GetBP()){sameBP=true;}
                            if (sameBP && sameStrand){
                                //G4cout<<"Merge found. Prim "<<hit1->GetPrimaryID()<<" "<<hit2->GetPrimaryID()<<" BP: "<<hit1->GetBP()<<" "<<hit2->GetBP()<<" EDEP "<<hit1->GetEnergy()<<" "<<hit2->GetEnergy()<<G4endl;
                                //merge found
                                toMerge=true;
                                hit2->SetChecked(true);
                                hit2->SetMarkedForDelete(true);
                                todelete++;
                                MergeList.push_back(hit2);
                            }
                        }
                    }
                }

            }
        }
        if (toMerge){
            //do the merge
            int nMerges=MergeList.size();
            G4double En=hit1->GetEnergy();
            bool isionis=hit1->GetIsIonis();
            G4ThreeVector avgPos=hit1->GetPosition();
            for (int k=0;k<nMerges;k++){
                En+=MergeList[k]->GetEnergy();
                if (MergeList[k]->GetIsIonis()){isionis=true;}
                avgPos+=MergeList[k]->GetPosition();
            }
            avgPos/=((double)nMerges+1.0); //+1 for the point that receives the merges
            hit1->SetEnergy(En);
            hit1->SetPosition(avgPos);
            hit1->SetIsIonis(isionis);
        }
        toMerge=false;
        hit1->SetChecked(true);
    }

    if (todelete>0){
        //delete the hit that was merged from the list
        int deleted=0;
        nHits=Hits.size();

        int i=Hits.size()-1;
        while (deleted!=todelete){
            //cout<<i<<endl;
            if (Hits[i]->GetMarkedForDelete()){

                delete Hits[i];
                Hits[i]=NULL;
                Hits.erase(Hits.begin()+i);

                deleted++;
                //G4cout<<i<<G4endl;
            }
            i--;
            if (i<0){ cout<<"Might be stuck "<<i<<endl; i=Hits.size()-1;} //hopefully never uses, else stuck in loop
        }
        //cout<<endl;
        //G4cout<<"Todelete "<<todelete<<" deleted "<<deleted<<G4endl<<G4endl;
    }
}

//check the hit is accepted (edep in volume)
bool StartClusteringBasePair::EnergyAccept(HitPoint*Hit)
{
    G4double edep=Hit->GetEnergy();

    //energy within a range
    if (EneRange){
        G4double EMinDamage = Emin/eV;
        G4double EMaxDamage = Emax/eV;

        if (edep<EMinDamage){
            Hit->SetMarkedForDelete(true);
            return false;
        }
        if (edep>EMaxDamage){return true;}

        G4double proba = (edep - EMinDamage)/ (EMaxDamage-EMinDamage);
        if (G4UniformRand() > proba){
            Hit->SetMarkedForDelete(true);
            return false;
        }
        else {return true;}
    }

    //energy above a threshold
    else if (EneThresh){
        if (edep>EThresh/eV){
            return true;
        } else {
            Hit->SetMarkedForDelete(true);
            return false;
        }
    }

    //ionisation inside DNA
    else if (EneIonis){
        if (Hit->GetIsIonis()){return true;}
        else {
            Hit->SetMarkedForDelete(true);
            return false;
        }
    }

    //PANIC -> no option selected
    else{
        G4cout<<"NO OPTION SELECTED TO DETERMINE DNA DAMAGE"<<G4endl;
        Hit->SetMarkedForDelete(true);
        return false;
    }

}

//take the hit list and split into backbones and bases
void StartClusteringBasePair::SplitList(vector<HitPoint*>&Backbones, vector<HitPoint*>&Bases, vector<HitPoint*>&Hits)
{
    int nHits=Hits.size();
    for (int i=0;i<nHits;i++){
        EnergyAccept(Hits[i]);
    }

    //remove the hits that dont meet energy criteria
    int n=Hits.size()-1;
    while (n>-1){
        if (Hits[n]->GetMarkedForDelete()){
            delete Hits[n];
            Hits[n]=NULL;
            Hits.erase(Hits.begin()+n);
        } else {
            if (Hits[n]->GetIsBack()){Backbones.push_back(Hits[n]);}
            else if (Hits[n]->GetIsBase()){Bases.push_back(Hits[n]);}
        }
        n--;
    }
}
