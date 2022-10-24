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
// Created by John-William Warmenhoven on 29/10/18.
// DaMaRiS is developed at the University of Manchester.
// See README for references.
//
//

#pragma once

#include <G4VITDiscreteProcess.hh>
#include "DrBreakMolecule.hh"

class G4MoleculeDefinition;

class DrProteinKinetics_Generic: public G4VITDiscreteProcess
{
public:
    DrProteinKinetics_Generic(G4String, G4MoleculeDefinition*, std::vector<G4MoleculeDefinition*>, G4double, G4bool, G4int);
    ~DrProteinKinetics_Generic();
    G4int fVerbose;
    void StartTracking(G4Track*);
    G4VParticleChange* PostStepDoIt(const G4Track& track, const G4Step& step);
    G4bool IsApplicable(const G4ParticleDefinition&);
    G4double PostStepGetPhysicalInteractionLength(const G4Track& track, G4double previousStepSize, G4ForceCondition* condition);

private:
    DrProteinKinetics_Generic(const DrProteinKinetics_Generic &right);
    G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition*);
    G4MoleculeDefinition* fFromMolecule;
    std::vector<G4MoleculeDefinition*> fToMolecule;
    G4double fProcessTime;
    unsigned fNumberSecondaries;
    G4String fProcessName;
    G4bool fCheckAddLesions;
    //0 = no, 1 = backbone, 2 = base
    G4int fCleanAddLesions;

protected:
    G4VParticleChange* PKGenericAction(const G4Track& track);
    G4VParticleChange aParticleChange;
    struct PKGenericState : public G4ProcessState
    {
        PKGenericState();
        virtual ~PKGenericState(){}
        G4double fPreviousTimeAtPreStepPoint;
    };
};
