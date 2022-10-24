//
// ********************************************************************
// *																  *
// * This file is part of the TOPAS-nBio extensions to the			  *
// *   TOPAS Simulation Toolkit.									  *
// * The TOPAS-nBio extensions are freely available under the license
// *   agreement set forth at: https://topas-nbio.readthedocs.io/	  *
// *																  *
// ********************************************************************
//
// Authors: Alejandro Bertolet, Jan Schuemann

#ifndef TsHitInDNA_hh
#define TsHitInDNA_hh

#include "G4UnitsTable.hh"

class TsHitInDNA
{
public:
	TsHitInDNA();
	~TsHitInDNA();

	// Get methods
	inline G4int GetEventID()						{ return fEventID; }
	inline G4double GetEdep()						{ return fEdep; }
	inline G4String GetParticleName()				{ return fParticleName; }
	inline std::vector<G4int> GetHierarchicalIDs()	{ return fHierarchicalIDs; }
	inline G4int GetStrandNumber()					{ return fStrandID; }
	inline G4int GetDNAComponentID()				{ return fComponentID; }
	inline G4ThreeVector GetPosition()				{ return fPos; }
	inline G4String GetProcess() 					{ return fProcess;}
	inline G4int GetChromosomeID()					{ return fChromosomeID; }
	inline G4int GetBasePairID()					{ return fBasePairID; }
	inline G4int GetDamageType()					{ return fDamageType; }

	// Set methods
	inline void SetEventID(G4int id)						{ fEventID = id; }
	inline void SetEdep(G4double edep)						{ fEdep = edep; }
	inline void SetParticleName(G4String name)				{ fParticleName = name; }
	inline void SetHierarchicalIDs(std::vector<G4int> ids)	{ fHierarchicalIDs = ids; }
	inline void SetStrandNumber(G4int stnum)				{ fStrandID = stnum; }
	inline void SetDNAComponentID(G4int cid)				{ fComponentID = cid; }
	inline void SetPosition(G4ThreeVector pos)				{ fPos = pos; }
	inline void SetProcess (G4String p)				 		{ fProcess = p;}
	inline void SetChromosomeID(G4int chid)					{ fChromosomeID = chid; }
	inline void SetBasePairID(G4int bpid)					{ fBasePairID = bpid; }
	inline void SetDamageType(G4int d)						{ fDamageType = d; }

private:
	G4int fEventID;

	G4double fEdep;
	G4String fParticleName;

	std::vector<G4int> fHierarchicalIDs;
	G4int fStrandID;
	G4int fComponentID;
	G4ThreeVector fPos;
	G4String fProcess;

	G4int fChromosomeID;
	G4int fBasePairID;

	G4int fDamageType;
};

#endif
