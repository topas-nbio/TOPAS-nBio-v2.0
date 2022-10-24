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
// Created by John-William Warmenhoven on 27/10/18.
// DaMaRiS is developed at the University of Manchester.
// See README for references.
//

#pragma once

#include <G4MoleculeDefinition.hh>

class DrDSBEnd_Generic : public G4MoleculeDefinition
{
private:
    static DrDSBEnd_Generic * theInstance;
    //@@@@ Static pointer to a map containing all the instances
    static std::map<G4String, DrDSBEnd_Generic*> theInstanceMap;
    DrDSBEnd_Generic (){}
    ~DrDSBEnd_Generic (){}

public:
    static DrDSBEnd_Generic* Definition(G4String);
};
