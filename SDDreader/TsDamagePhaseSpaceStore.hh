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

#pragma once

#include "DrReaderSDD.hh"

//@@@@ This struct is only here as legacy code to maintain backwards
//@@@@ compatibility with older versions of code used to develop the
//@@@@ SDD. Feel free to remove if your code does not need to use.
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

	// [0] Same event [1] Start of damages associated with a new particle [2] Start of damages associated with a new dose
    G4int newEvent{-1};
	// Geometric position in nucleus
    G4ThreeVector position{0.,0.,0.};
	// [Number of base damages, number of single strand breaks, number of DSBs]
    std::vector<G4int> damageClass{-1,-1,-1};
	// [0] direct damage [1] indirect damage [2] combination
    G4int causeOfDamage{-1};
	// [Chromosome number, copy number, chromatid number]
    std::vector<G4int> chromosomeID{0,0,0};
	// Fraction along chromosome
    G4double chromPosition{-1.};
	// strand number [space] base pair number [space] type of break /
    G4String Full_Break_Spec{""};
	// 5' to 3' only
    G4String DNA_Sequence{""};
};

class TsDamagePhaseSpaceStore
{
public:
	static TsDamagePhaseSpaceStore* Instance();
	void Destroy();

    typedef std::vector<std::vector<std::vector<DrDamageEvent*> > > dmgstore;
	dmgstore& GetDamageEventStore();
    DrDamageHeader* GetDamageHeaderPtr();

    //-----------------------------------------------------------
    //@@@@ LEGACY CODE
    inline std::vector<damage_event> GetDamagePhaseSpace() const
    {return damagePhaseSpace;}
    inline void AddEvent(damage_event _event)
    {damagePhaseSpace.push_back(_event);}
    std::vector<damage_event> damagePhaseSpace;
    //@@@@ LEGACY CODE
    //-----------------------------------------------------------

protected:
    static TsDamagePhaseSpaceStore* fpgTsDamagePhaseSpaceStore;
    TsDamagePhaseSpaceStore();
	~TsDamagePhaseSpaceStore();
	DrDamageHeader* fDamageHeader{nullptr};
	dmgstore fDamageExposures;
};
