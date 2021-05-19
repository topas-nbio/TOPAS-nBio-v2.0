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

class DrFRAP: public G4VITDiscreteProcess{
public:
	DrFRAP();
	~DrFRAP();

	G4int fVerbose;

    void StartTracking(G4Track*);
	G4VParticleChange* PostStepDoIt(const G4Track& track, const G4Step& step);
	G4bool IsApplicable(const G4ParticleDefinition&);
	G4double PostStepGetPhysicalInteractionLength(const G4Track& track, G4double previousStepSize, G4ForceCondition* condition);
private:
	DrFRAP(const DrFRAP &right);
	G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*);

protected:
	G4VParticleChange* FRAP(const G4Track& track);
	G4VParticleChange aParticleChange;

	struct FRAPState : public G4ProcessState{
		FRAPState();
        virtual ~FRAPState(){}
        G4double fPreviousTimeAtPreStepPoint;
    };
};
