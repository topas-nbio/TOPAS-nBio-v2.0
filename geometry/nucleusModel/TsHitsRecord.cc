// Extra Class for TsNucleusScoring
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
// Authors: Hongyu Zhu

#include "TsHitsRecord.hh"

TsHitsRecord::TsHitsRecord()
{
        fx = 0. ;
        fy = 0. ;
        fz = 0. ;
        fPos = G4ThreeVector(fx, fy, fx);
        fEdep = 0. ;
        fParticleEnergy = 0. ;
        fParticleName  = "";
        fEventID = -1;
        fRunID = -1;
        fVolumeName = "";
        fChromosomeID = 0;
        fBasePairID   = 0;
        fIsDirectDamage = true;
        fProcess = "";

        flagDSB = false;
        flagSSB = false;
        DSBType       = -1;

        flagComplexity = false;
        flagSSB_P = false;
        flagDSB_P = false;
        flagDSB_PP = false;
}
    
TsHitsRecord::~TsHitsRecord()
{}
