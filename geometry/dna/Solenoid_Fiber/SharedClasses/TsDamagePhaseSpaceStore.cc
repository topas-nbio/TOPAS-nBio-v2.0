// Extra Class for use by BlahBlahBlah
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
 * TsDamagePhaseSpaceStore.cc
 *
 *  Created on: Oct 12, 2017
 *      Author: JWarmenhoven & NHenthorn
 */

#include "TsDamagePhaseSpaceStore.hh"

TsDamagePhaseSpaceStore* TsDamagePhaseSpaceStore::fpgTsDamagePhaseSpaceStore(0);

TsDamagePhaseSpaceStore* TsDamagePhaseSpaceStore::Instance(){
		if(!fpgTsDamagePhaseSpaceStore) fpgTsDamagePhaseSpaceStore = new TsDamagePhaseSpaceStore();
		return fpgTsDamagePhaseSpaceStore;
}

void TsDamagePhaseSpaceStore::Destroy(){
	if(!fpgTsDamagePhaseSpaceStore) delete fpgTsDamagePhaseSpaceStore;
	fpgTsDamagePhaseSpaceStore = NULL;
}

TsDamagePhaseSpaceStore::TsDamagePhaseSpaceStore(){}
TsDamagePhaseSpaceStore::~TsDamagePhaseSpaceStore(){}
