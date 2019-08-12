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
 * TsDamagePhaseSpaceStore.hh
 *
 *  Created on: Oct 12, 2017
 *      Author: JWarmenhoven & NHenthorn
 */

#ifndef INCLUDE_TSDAMAGEPHASESPACESTORE_HH_
#define INCLUDE_TSDAMAGEPHASESPACESTORE_HH_

#include <vector>
#include "G4ios.hh"
#include "G4UnitsTable.hh"

struct damage_event{

	//----------------------------Set Functions---------------------------------
	inline void SetNewEvent(G4int _newEvent) {newEvent = _newEvent;}
	inline void SetPosition(G4ThreeVector _position) {position = _position;}
	inline void SetNumberOfBases(G4int _bases) {damageClass[0] = _bases;}
	inline void SetNumberOfBackbones(G4int _backbones)
			{damageClass[1] = _backbones;}
	inline void SetNumberOfDSB(G4int _DSBs) {damageClass[2] = _DSBs;}
	inline void SetCauseOfDamage(G4int _causeOfDamage)
			{causeOfDamage = _causeOfDamage;}
	inline void SetChromosomeNumber(G4int _chromosomeNumber)
			{chromosomeID[0] = _chromosomeNumber;}
	inline void SetChromosomeCopyNumber(G4int _chromosomeCopyNumber)
			{chromosomeID[1] = _chromosomeCopyNumber;}
	inline void SetChromatidNumber(G4int _chromatidNumber)
			{chromosomeID[2] = _chromatidNumber;}
	inline void SetChromosomePosition(G4double _chromosomePosition)
	 		{chromPosition = _chromosomePosition;}
	inline void SetFullBreakSpec(G4String _fullBreakSpec)
			{Full_Break_Spec = _fullBreakSpec;}
	inline void SetDNASequence(G4String _DNASequence)
			{DNA_Sequence = _DNASequence;}
	//----------------------------Get Functions---------------------------------
	inline G4int GetNewEvent() const {return newEvent;}
	inline G4ThreeVector GetPosition() const {return position;}
	inline G4int GetNumberOfBases() const {return damageClass[0];}
	inline G4int GetNumberOfBackbones() const {return damageClass[1];}
	inline G4int GetNumberOfDSB() const {return damageClass[2];}
	inline G4int GetCauseOfDamage() const {return causeOfDamage;}
	inline G4int GetChromosomeNumber() const {return chromosomeID[0];}
	inline G4int GetChromosomeCopyNumber() const {return chromosomeID[1];}
	inline G4int GetChromatidNumber() const {return chromosomeID[2];}
	inline G4double GetChromosomePosition() const {return chromPosition;}
	inline G4String GetFullBreakSpec() const {return Full_Break_Spec;}
	inline G4String GetDNASequence() const {return DNA_Sequence;}
	//------------Data fields because John likes POD Structures-----------------
	G4int newEvent{-1};															// [0] Same event [1] Start of damages associated with a new particle [2] Start of damages associated with a new dose
	G4ThreeVector position{0.,0.,0.};											// Geometric position in nucleus
	std::vector<G4int> damageClass{-1,-1,-1};									// [Number of base damages, number of single strand breaks, number of DSBs]
	G4int causeOfDamage{-1};													// [0] direct damage [1] indirect damage [2] combination
	std::vector<G4int> chromosomeID{0,0,0};										// [Chromosome number, copy number, chromatid number]
	G4double chromPosition{-1.};												// Fraction along chromosome
	G4String Full_Break_Spec{""};												// strand number [space] base pair number [space] type of break /
	G4String DNA_Sequence{""};													// 5' to 3' only
};

class TsDamagePhaseSpaceStore
{
public:
	static TsDamagePhaseSpaceStore* Instance();
	inline std::vector<damage_event> GetDamagePhaseSpace() const
			{return damagePhaseSpace;}
	inline void AddEvent(damage_event _event)
			{damagePhaseSpace.push_back(_event);}
	void Destroy();
	std::vector<damage_event> damagePhaseSpace;									// Here for now because Nick accesses it directly in current version of code, move to protected

protected:
    static TsDamagePhaseSpaceStore* fpgTsDamagePhaseSpaceStore;
    TsDamagePhaseSpaceStore();
	~TsDamagePhaseSpaceStore();
};

#endif /* INCLUDE_TSDAMAGEPHASESPACESTORE_HH_ */
