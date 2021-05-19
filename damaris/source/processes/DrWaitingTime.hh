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

#include <G4VITDiscreteProcess.hh>

class DrWaitingTime : public G4VITDiscreteProcess {
public:
	DrWaitingTime();
	virtual ~DrWaitingTime();

    void StartTracking(G4Track*);
	G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);
	G4bool IsApplicable(const G4ParticleDefinition&);
	G4double PostStepGetPhysicalInteractionLength(const G4Track&,G4double,
													G4ForceCondition*);
private:
	G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*);
	void Process(const G4Track&);
protected:
	void GenerateWaitingTime(const G4Track&);

	G4VParticleChange aParticleChange;
    struct WaitingTimeState : public G4ProcessState {
		WaitingTimeState();
        virtual ~WaitingTimeState(){}
        G4double fPreviousTimeAtPreStepPoint{-1};
        G4bool fHasJustProc{false};
    };
};
