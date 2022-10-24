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

class DrReportMSD : public G4VITDiscreteProcess {
public:
  DrReportMSD();
  ~DrReportMSD();

  G4int fVerbose;

  void StartTracking(G4Track *);
  G4VParticleChange *PostStepDoIt(const G4Track &track, const G4Step &step);
  G4bool IsApplicable(const G4ParticleDefinition &);
  G4double PostStepGetPhysicalInteractionLength(const G4Track &track, G4double previousStepSize, G4ForceCondition *condition);

private:
  DrReportMSD(const DrReportMSD &right);
  G4double GetMeanFreePath(const G4Track &, G4double, G4ForceCondition *);

protected:
  G4VParticleChange *RecordMSD(const G4Track &track);
  G4VParticleChange aParticleChange;
  struct ReportMSDState : public G4ProcessState {
    ReportMSDState();
    virtual ~ReportMSDState() {}
    G4double fPreviousTimeAtPreStepPoint;
  };
};
