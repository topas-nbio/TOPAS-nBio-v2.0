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
//  AddAddChromosomeDomain.cpp
//  Clustering
//
//  Created by Nick Henthorn on 27/10/2017.
/*
 Using the methodology of:
 A Model for Interphase Chromosomes and Evaluation of Radiation-Induced Aberrations
 W . R. Holley, I. S. Mian, S. J. Park, B. Rydberg and A. Chatterjee
 */


#include "AddChromosomeDomain.hh"


void AddChromosomeDomain::Start(map<int, map<int, Clusters*>> &ClusterMap,
                                double NucRad)
{
    Verbosity=0;    //>0 will output stuff to screen

    NuclearRadius=NucRad;   //radius of bounding sphere (um)

    MaxGap=0.0;          //each domain must touch at least one other
    MaxOverlapFraction=0.9; //The max overlap as a fraction of the radius of the smallest domain in the comparison

    //set up 46 spherical domains
    ConstructTerritories();

    for (auto data:ClusterMap){
        for (auto d:data.second){
            //isolated damage
            if (d.second->GetClusterID()==0){
                for (size_t i=0;i<d.second->GetHits().size();i++){
                    double x_hit=d.second->GetHits()[i]->GetPosition().x()/1000.0;
                    double y_hit=d.second->GetHits()[i]->GetPosition().y()/1000.0;
                    double z_hit=d.second->GetHits()[i]->GetPosition().z()/1000.0;

                    int chromosome=-1, copy=-1, chromatid=-1;
                    AddChromID(x_hit,
                               y_hit,
                               z_hit,
                               chromosome,
                               copy,
                               chromatid);
                    d.second->GetHits()[i]->SetChromID(chromosome);
                    d.second->GetHits()[i]->SetChromCopy(copy);
                    d.second->GetHits()[i]->SetChromatid(chromatid);
                }

            }
            //DSB and comp SSB
            else {
                double x_hit=d.second->GetCentre().x()/1000.0;
                double y_hit=d.second->GetCentre().y()/1000.0;
                double z_hit=d.second->GetCentre().z()/1000.0;

                int chromosome=-1, copy=-1, chromatid=-1;
                AddChromID(x_hit,
                           y_hit,
                           z_hit,
                           chromosome,
                           copy,
                           chromatid);

                d.second->SetChromID(chromosome);
                d.second->SetChromCopy(copy);
                d.second->SetChromatid(chromatid);

            }
        }
    }

    cout<<endl;


}

void AddChromosomeDomain::AddChromID(double x_hit,
                                     double y_hit,
                                     double z_hit,
                                     int &chromosome,
                                     int &copy,
                                     int &chromatid)
{
        vector<int> PossibleDomains(0);

        //check if the point is within a domain
        for (auto d:DomainCentre){
            double LHS = pow((x_hit-d.second[0]),2.0) +
                         pow((y_hit-d.second[1]),2.0) +
                         pow((z_hit-d.second[1]),2.0);
            double RHS = DomainRadii[d.first]*DomainRadii[d.first];

            if (LHS<=RHS){
                PossibleDomains.push_back(d.first);
            }
        }

        //if one domain, set it to that
        if (PossibleDomains.size()==1){
            ConvertIDtoChromID(PossibleDomains[0],
                               chromosome,
                               copy,
                               chromatid);
        }

        //if not in a domain, pick closest one
        else if (PossibleDomains.size()==0){
            double minSep=numeric_limits<double>::max();
            int ClosestID=-1;
            for (auto d:DomainCentre){
                double separation = sqrt(pow((x_hit-d.second[0]),2.0) +
                                         pow((y_hit-d.second[1]),2.0) +
                                         pow((z_hit-d.second[2]),2.0));
                separation-=DomainRadii[d.first];

                if (separation<minSep){
                    minSep=separation;
                    ClosestID=d.first;
                }
            }
            ConvertIDtoChromID(ClosestID,
                               chromosome,
                               copy,
                               chromatid);
        }

        //if in overlap of 2 or more domains, randomly pick
        else if (PossibleDomains.size()>0){
            map<int, double> ProbInDomain;
            double sum=0.0;

            //calculate distance from hit to centre of each domain
            for (size_t j=0;j<PossibleDomains.size();j++){
                double separation = sqrt(pow((x_hit-DomainCentre[PossibleDomains[j]][0]),2.0) +
                                         pow((y_hit-DomainCentre[PossibleDomains[j]][1]),2.0) +
                                         pow((z_hit-DomainCentre[PossibleDomains[j]][2]),2.0));
                sum+=separation;

                ProbInDomain[PossibleDomains[j]]=separation;
            }


            //normalise separations to sum to 1
            for (auto d:ProbInDomain){
                double renorm = d.second/sum;
                ProbInDomain[d.first] = renorm;
            }


            //change to cumulative
            double newSum=0.0;
            for (auto d:ProbInDomain){
                newSum+=d.second;
                ProbInDomain[d.first]=newSum;
            }

            //randomly choose, with probs weighted by separation to centre of domain
            double Rando = G4UniformRand();
            int ID=-1;
            for (auto d:ProbInDomain){
                if (Rando<=d.second){
                    ID=d.first;
                    break;
                }
            }
            ConvertIDtoChromID(PossibleDomains[0],
                               chromosome,
                               copy,
                               ID);
        }

}

void AddChromosomeDomain::ConstructTerritories()
{
    //Chromosome Size (Mbp)
    //http://www.ncbi.nlm.nih.gov/mapview/maps.cgi?ORG=hum&MAPS=ideogr,est,loc&LINKS=ON&VERBOSE=ON&CHR=1
    int MBP[46];

    //chrom 1->23 copy 1 = 0->22 (23=X)
    MBP[0]=250; MBP[1]=242; MBP[2]=198; MBP[3]=190; MBP[4]=182; MBP[5]=171; MBP[6]=159;
    MBP[7]=145; MBP[8]=138; MBP[9]=134; MBP[10]=135; MBP[11]=133; MBP[12]=114;
    MBP[13]=107; MBP[14]=102; MBP[15]=90; MBP[16]=83; MBP[17]=80; MBP[18]=59;
    MBP[19]=64; MBP[20]=47; MBP[21]=51; MBP[22]=156;

    //chrom 1->23 copy 2 = 23->45 (45=Y)
    MBP[23]=250; MBP[24]=242; MBP[25]=198; MBP[26]=190; MBP[27]=182; MBP[28]=171; MBP[29]=159;
    MBP[30]=145; MBP[31]=138; MBP[32]=134; MBP[33]=135; MBP[34]=133; MBP[35]=114;
    MBP[36]=107; MBP[37]=102; MBP[38]=90; MBP[39]=83; MBP[40]=80; MBP[41]=59;
    MBP[42]=64; MBP[43]=47; MBP[44]=51; MBP[45]=57;
    double c_Genom=0.0;
    for (int i=0;i<46;i++){
        c_Genom+=(double)MBP[i];
    }

    //Domain Radius
    //Radial Arrangement of Chromosome Territories in Human Cell Nuclei: A Computer Model Approach Based on Gene Density Indicates a Probabilistic Global Positioning Code
    // Kreth et al
    //Radius = R_Nuc * pow((v*(c_Chrom/c_Genom)),(1.0/3.0));
    //v = 1;

    for (int i=0;i<46;i++){
        DomainRadii[i]=NuclearRadius*pow(1.0*(double)MBP[i]/c_Genom, (1.0/3.0));

        if (Verbosity>1){
            cout<<endl<<"Domain "<<i<<", Radius = "<<DomainRadii[i]<<" um"<<endl;
        }

    }

    if (Verbosity>1){
        cout<<endl<<"max Overlap "<<MaxOverlap<<" um"<<endl<<endl;
    }

    //randomly place the territory in spherical nucleus
    bool PlacedAll=false;
    int AttemptCount=0;
    while (!PlacedAll){
        PlacedAll=RandomlyPlace();

        AttemptCount++;
        if (Verbosity>0 && !PlacedAll){
            cout<<"Failed To Place All Domains"<<endl;
            cout<<"Trying Again "<<AttemptCount<<endl;
            cout<<endl;
        }
    }

}

bool AddChromosomeDomain::RandomlyPlace()
{
    int failCount=0;
    int failLimit=10000;

    for (auto DRadiiData:DomainRadii){
        bool acceptPosition=false;
        double x=0.0, y=0.0, z=0.0;
        while (!acceptPosition){
            //Generate a random radii
            double RandRadius = (NuclearRadius-DRadiiData.second) * G4UniformRand();

            //Generate a random Theta
            double Theta =G4UniformRand() * (180.0);

            //Generate a random Phi
            double Phi =G4UniformRand() * (360.0);

            x = RandRadius*sin(Theta)*cos(Phi);
            y = RandRadius*sin(Theta)*sin(Phi);
            z = RandRadius*cos(Theta);

            acceptPosition=CheckDomain(DRadiiData.first,
                                       x,
                                       y,
                                       z,
                                       DRadiiData.second);

            if (failCount>failLimit){
                if (Verbosity>0){
                    cout<<"Failed too much, Clearing"<<endl;
                }
                DomainCentre.clear();
                return false;
            }

            failCount++;
        }

        if (Verbosity>1){
            cout<<"Placed Domain: "<<DRadiiData.first<<" x:"<<x<<" y:"<<y<<" z:"<<z<<endl<<endl;
        }

        DomainCentre[DRadiiData.first][0]=x;
        DomainCentre[DRadiiData.first][1]=y;
        DomainCentre[DRadiiData.first][2]=z;

    }

    return true;

}


bool AddChromosomeDomain::CheckDomain(int DomainID,
                                   double x,
                                   double y,
                                   double z,
                                   double Radius)
{
    //bool OverlapCondition=false;
    //bool TouchingCondition=false;

    //check overlap with other domains
    for (auto position:DomainCentre){

        //dont check yourself (you may wreck yourself)
        if (position.first!=DomainID){
            double x_prev = position.second[0];
            double y_prev = position.second[1];
            double z_prev = position.second[2];

            double separation = sqrt(pow((x-x_prev),2.0) +
                                     pow((y-y_prev),2.0) +
                                     pow((z-z_prev),2.0));
            double Radius_prev = DomainRadii[position.first];

            double overlap = separation - Radius - Radius_prev;
            // <0 -> overlap
            // >0 -> no overlap
            // =0 -> touching

            //set max overlap to a fraction of the smaller domain
            double minRadius=Radius;
            if (Radius_prev<Radius){minRadius=Radius_prev;}
            MaxOverlap=minRadius*MaxOverlapFraction;

            //overlapping by more than max
            if (overlap<=-MaxOverlap){

                //OverlapCondition=false;
                //TouchingCondition=false;
                if (Verbosity>1){
                    cout<<"Checking Domain: "<<DomainID<<", Against: "<<position.first<<endl;
                    cout<<"Proposed Centre "<<x<<" "<<y<<" "<<z<<" Radius "<<Radius<<endl;
                    cout<<"Comparis Centre "<<position.second[0]<<" "<<position.second[1]<<" "<<position.second[2]<<" Radius "<<Radius_prev<<endl;
                    cout<<"Too much overlap, "<<-overlap<<", max: "<<MaxOverlap<<endl<<endl;
                }
                return false; //overlapping too much with a volume. Start again
            }

            //not touching another domain
            if (overlap>0.0){
                if (Verbosity>1){
                    cout<<"Checking Domain: "<<DomainID<<" Against "<<position.first<<endl;
                    cout<<"Proposed Centre "<<x<<" "<<y<<" "<<z<<" Radius "<<Radius<<endl;
                    cout<<"Comparis Centre "<<position.second[0]<<" "<<position.second[1]<<" "<<position.second[2]<<" Radius "<<Radius_prev<<endl;
                    cout<<"Not touching a domain"<<endl<<endl;
                }
            //} else {
            //    TouchingCondition=true;
            }
        }
    }

    return true;
    //    if (TouchingCondition){
    //        return true;
    //    } else {
    //        return false;
    //    }

}

void AddChromosomeDomain::ConvertIDtoChromID(int id, int &Chromosome, int &Copy, int &chromatid)
{
    map<int, int> ChromID;
    for (int i=0;i<46;i++){
        int TempID=-1;
        if (i>22){
            TempID=i+1-23;
        } else {
            TempID=i+1;
        }
        ChromID[i]=TempID;
    }
    Chromosome=ChromID[id];

    map<int,int> CopyID;
    for (int i=0;i<46;i++){
        int CID=1;
        if (i>22){
            CID=2;
        }
        CopyID[i]=CID;
    }
    Copy=CopyID[id];

    if ( G4UniformRand()>0.5){
        chromatid=2;
    }else{
        chromatid=1;
    }

}
