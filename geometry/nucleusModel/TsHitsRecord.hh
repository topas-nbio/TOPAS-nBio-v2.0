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

#ifndef TsHitsRecord_hh
#define TsHitsRecord_hh


#include "G4UnitsTable.hh"
#include <stdint.h>

class TsHitsRecord 
{
public:
    TsHitsRecord();
    ~TsHitsRecord();

    // Get method
    inline G4ThreeVector GetPosition() {return fPos;}
    inline G4double GetPosX() {return fx;}
    inline G4double GetPosY() {return fy;}
    inline G4double GetPosZ() {return fz;}   
    inline G4double GetEdep() {return fEdep;}
    inline G4double GetParticleEnergy() {return fParticleEnergy;}
    inline G4int GetChromosomeID() {return fChromosomeID;}
    inline G4int GetBasePairID() {return fBasePairID;}
    inline G4int GetEventID() {return fEventID;}
    inline G4int GetRunID() {return fRunID;}
	inline G4int GetVoxelID() {return fVoxelID;}
    inline G4int GetDSBType() {return DSBType;}
    inline G4bool GetIsDirectDamage () {return fIsDirectDamage;}
	inline G4bool GetFlagDSB() {return flagDSB;}
	inline G4bool GetFlagSSB() {return flagSSB;}
    inline G4bool GetFlagSSB_P() {return flagSSB_P;}
    inline G4bool GetFlagDSB_P() {return flagDSB_P;}
    inline G4bool GetFlagDSB_PP() {return flagDSB_PP;}
    inline G4String GetParticleName() {return fParticleName;}
	inline G4String GetVolumeName() {return fVolumeName;}
	inline G4String GetProcess() {return fProcess;}


    // Set method
    inline void SetPosition (G4ThreeVector aPos) {fPos = aPos;}
	inline void SetPosX(G4double aPosX) {fx = aPosX;}
	inline void SetPosY(G4double aPosY) {fy = aPosY;}
	inline void SetPosZ(G4double aPosZ) {fz = aPosZ;}
    inline void SetEdep (G4double aEdep) {fEdep = aEdep;}
    inline void SetParticleEnergy (G4double aParticleEnergy) {fParticleEnergy = aParticleEnergy;}
    inline void SetParticleName (G4String aParticleName) {fParticleName = aParticleName;}
    inline void SetChromosomeID(G4int aID) {fChromosomeID = aID;}
    inline void SetBasePairID(G4int aID) {fBasePairID = aID;}
    inline void SetEventID (G4int aEventID) {fEventID = aEventID;}
    inline void SetRunID (G4int aRunID) {fRunID = aRunID;}
	inline void SetVoxelID(G4int aVoxelID) {fVoxelID = aVoxelID;}
    inline void SetVolumeName (G4String aVolumeName) {fVolumeName = aVolumeName;}
    inline void SetIsDirectDamage (G4bool aIsDirectDamage) {fIsDirectDamage = aIsDirectDamage;}
	inline void SetFlagDSB (G4bool aflag) {flagDSB = aflag;}
	inline void SetFlagSSB (G4bool aflag) {flagSSB = aflag;}
    inline void SetFlagSSB_P (G4bool aflag) {flagSSB_P  = aflag;}
    inline void SetFlagDSB_P (G4bool aflag) {flagDSB_P  = aflag;}
    inline void SetFlagDSB_PP(G4bool aflag) {flagDSB_PP = aflag;}
    inline void SetDSBType (G4int type) {DSBType= type;}
    inline void SetProcess (G4String p) {fProcess = p;}

    
protected:

    
private:
	G4double fEdep ;  //edep (eV)
	G4double fParticleEnergy;
	G4String fParticleName;
	G4String fVolumeName ;
	G4int fChromosomeID;
	G4int fBasePairID;
	G4int fEventID ;
	G4int fRunID ;
	G4int fVoxelID ;
	G4double fx ;  //x (nm)
	G4double fy ;  //y (nm)
	G4double fz ;  //z (nm)
	G4ThreeVector fPos;
	G4bool fIsDirectDamage; // true: direct damage; false: indirect damage
	G4String fProcess;

	G4bool flagDSB;
	G4bool flagSSB;
	G4int  DSBType;        // -1: not a DSB; 0: dir DSB; 1: indir DSB; 2: hybrid DSB

	G4bool flagComplexity;
	G4bool flagSSB_P;   // SSB+
	G4bool flagDSB_P;   // DSB+
	G4bool flagDSB_PP;  // DSB++
};

#endif
