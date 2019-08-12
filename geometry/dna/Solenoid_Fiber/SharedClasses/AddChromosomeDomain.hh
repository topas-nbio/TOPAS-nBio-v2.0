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
//  AddChromosomeDomain.h
//  Clustering
//
//  Created by Nick Henthorn on 27/10/2017.

#ifndef AddChromosomeDomain_h
#define AddChromosomeDomain_h


#include <stdio.h>
#include "HitPoint.hh"
#include "Clusters.hh"
#include <vector>
#include <map>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

class AddChromosomeDomain{
public:
    AddChromosomeDomain(){}
    ~AddChromosomeDomain(){}

    void Start(map<int, map<int, Clusters*>> &ClusterMap,
               double NucRad);

    map<int,double> GetDomainRadii(){return DomainRadii;}
    map<int,double[3]> GetDomainCentre(){return DomainCentre;}

private:
    void AddChromID(double x_hit,
                    double y_hit,
                    double z_hit,
                    int &chromosome,
                    int &copy,
                    int &chromatid);

    void ConstructTerritories();
    bool RandomlyPlace();
    bool CheckDomain(int DomainID,
                     double x,
                     double y,
                     double z,
                     double Radius);

    map<int,double> DomainRadii;
    map<int,double[3]> DomainCentre;

    void ConvertIDtoChromID(int id, int &Chromosome, int &Copy, int &chromatid);

    int Verbosity=0;

    double NuclearRadius=2.5;   //size of bounding volume (um)
    double MaxOverlap=0.100;    //how much overlap between domains (um) [0.675 from literature]
    double MaxGap=0.0;          //each domain must touch at least one other
    double MaxOverlapFraction = 0.1;    //fractional overlap for the smallest domain

};



#endif /* AddChromosomeDomain_h */
