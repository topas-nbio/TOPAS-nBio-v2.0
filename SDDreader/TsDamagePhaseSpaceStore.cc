// Extra Class for use by TOPAS
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
// Created by John-William Warmenhoven and Nicholas Henthorn on 10/2017.
// Please cite: https://doi.org/10.1667/RR15209.1
//

#include "TsDamagePhaseSpaceStore.hh"

TsDamagePhaseSpaceStore* TsDamagePhaseSpaceStore::fpgTsDamagePhaseSpaceStore(0);

TsDamagePhaseSpaceStore* TsDamagePhaseSpaceStore::Instance(){
		if(!fpgTsDamagePhaseSpaceStore) fpgTsDamagePhaseSpaceStore = new TsDamagePhaseSpaceStore();
		return fpgTsDamagePhaseSpaceStore;
}

void TsDamagePhaseSpaceStore::Destroy(){
	if(!fpgTsDamagePhaseSpaceStore) delete fpgTsDamagePhaseSpaceStore;
	fpgTsDamagePhaseSpaceStore = nullptr;
}

TsDamagePhaseSpaceStore::TsDamagePhaseSpaceStore(){}

TsDamagePhaseSpaceStore::~TsDamagePhaseSpaceStore(){
    delete fDamageHeader;
    for(auto level_1: fDamageExposures){
        for(auto level_2: level_1){
            for(auto level_3: level_2){
                delete level_3;
            }
        }
    }
}

DrDamageHeader* TsDamagePhaseSpaceStore::GetDamageHeaderPtr() {
    if (fDamageHeader != nullptr) return fDamageHeader;
    else {
        fDamageHeader = new DrDamageHeader();
        return fDamageHeader;
    }
}

std::vector<std::vector<std::vector<DrDamageEvent*> > >& TsDamagePhaseSpaceStore::GetDamageEventStore() {
    return fDamageExposures;
}
