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
// Created by John-William Warmenhoven.
// DaMaRiS is developed at the University of Manchester.
// See README for references.
//
#pragma once

#include <G4MoleculeDefinition.hh>

class DrClock : public G4MoleculeDefinition
{
private:
	static DrClock* theInstance;
	DrClock(){}
	~DrClock(){}

public:
	static DrClock* Definition();
};
