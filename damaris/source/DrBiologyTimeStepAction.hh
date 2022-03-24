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

#include <G4UserTimeStepAction.hh>

class TsParameterManager;

class DrBiologyTimeStepAction : public G4UserTimeStepAction
{
public:
	DrBiologyTimeStepAction(TsParameterManager*, G4String);
	virtual ~DrBiologyTimeStepAction();

	DrBiologyTimeStepAction(const DrBiologyTimeStepAction& other);
	DrBiologyTimeStepAction& operator=(const DrBiologyTimeStepAction& other);

	virtual void StartProcessing(){;}
	virtual void UserPreTimeStepAction();
	virtual void UserPostTimeStepAction();
	virtual void UserReactionAction(const G4Track&,
									const G4Track&,
									const std::vector<G4Track*>*);
	virtual void EndProcessing(){;}

private:
	TsParameterManager* fPm;
    void DebugPositions();
    void PlotPositions();
    void ExtractDisplacements();
    void RunFinalCheckBreaks();
};
