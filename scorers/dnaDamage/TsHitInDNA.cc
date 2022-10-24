// Extra Class for TsScoreDNADamageSBS
//
// ********************************************************************
// *																  *
// * This file is part of the TOPAS-nBio extensions to the			  *
// *   TOPAS Simulation Toolkit.									  *
// * The TOPAS-nBio extensions are freely available under the license *
// *   agreement set forth at: https://topas-nbio.readthedocs.io/	  *
// *																  *
// ********************************************************************
//
// Authors: Alejandro Bertolet, Jan Schuemann

#include "TsHitInDNA.hh"

TsHitInDNA::TsHitInDNA()
{
	fEventID = -1;

	fEdep = 0.;
	fParticleName = "";

	fStrandID = -1;
	fComponentID = -1;
	fPos = G4ThreeVector(0., 0., 0.);

	fChromosomeID = -1;
	fBasePairID = -1;

	fDamageType = -1;
}

TsHitInDNA::~TsHitInDNA() {}
