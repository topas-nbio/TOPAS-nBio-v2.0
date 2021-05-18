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

#include <G4MoleculeDefinition.hh>
#include "DrProteinKinetics_Generic.hh"

class TsParameterManager;
class DrDefinitions {

    //-----------------------------------------------
    //                  Functions
    //-----------------------------------------------

protected:
    DrDefinitions();
    static DrDefinitions *definitionsInstance;

private:
    G4String GetFullParmName(const G4String&);
    void SetParameters();
    std::vector<G4String> RegisterMolecules(G4String);

public:
    ~DrDefinitions();
    static DrDefinitions *Instance();

    inline void ResetCurrentExplicitBinNuber() { fCurrentExplicitBinNuber = -1; }
    inline void IncrimentCurrentExplicitBinNuber() { fCurrentExplicitBinNuber++; }

    //-----------------Init Functions-----------------
    void Initialise(TsParameterManager*);
    void InitialiseParameters();
    void InitialisePathway();
    std::vector<std::pair<G4double, std::vector<G4String>>> InitialiseReactions();
    void InitialiseBinning();

    //-----------------Get Functions-----------------
    inline G4double GetJumpDiff() const { return fDiffusionCoefficientForJump; }
    inline G4double GetTrapDiff() const { return fDiffusionCoefficientForTrapped; }
    inline G4double GetMinumumWaitingTime() const { return fMinWaitingTime; }
    inline G4bool GetIsSubDiffusion() const { return fIsSubDiffusion; }
    inline G4double GetBoundingCellOrNucleusRadius() const { return fBoundingRadius; }
    inline std::vector<std::pair<G4MoleculeDefinition*,DrProteinKinetics_Generic*> >
            GetProcList() const {return fProcList;}
    inline G4String GetDamageFile() const { return fDamageFileName; }
    inline G4int GetSelectFromExposure() const { return fSelectFromExposure; }
    inline G4int GetBiologyRepeatNumber() const { return fBiologyRepeatNumber; }
    inline G4double GetBiologyEndTime() const { return fBiologyEndTime; }
    inline G4bool GetIsReadFromFile() const { return fIsReadFromFile; }
    inline G4int GetAlternativeRunID() const { return fAlternativeRunID; }
    inline G4double GetObserveDurationForMSD() const { return fObserveDurationForMSD; }
    inline G4double GetObserveStepSizeForMSD() const { return fObserveStepSizeForMSD; }
    inline std::vector<G4String> GetHasKu() const { return fHasKu; }
    inline std::vector<G4String> GetHasPKcs() const { return fHasPKcs; }
    inline std::vector<G4String> GetHasXRCC4() const { return fHasXRCC4; }
    inline G4double GetDSBSeparation() const { return fDSBSeparation; }
    inline G4double GetDSBTimeDelay() const { return fDSBTimeDelay; }
    inline std::vector<G4int> GetDSBPattern() const { return fDSBPattern; }
    inline G4double GetDSBOffset() const { return fDSBOffset; }
    inline G4int GetDSBOriginNumber() const { return fDSBOriginNumber; }
    inline G4int GetDSBColumnNumber() const { return fDSBColumnNumber; }
    inline G4double GetDSBColumnRadius() const { return fDSBColumnRadius; }
    inline G4String GetExplicitBinning() const { return fExplicitBinning; }
    inline std::vector<G4double> GetExplicitBins() const { return fExplicitBins; }
    inline std::vector<G4double> GetExplicitBinDifferences() const { return fExplicitBinDiff; }
    inline G4int GetCurrentExplicitBinNuber() const { return fCurrentExplicitBinNuber; }
    inline G4bool GetTurnOffTime() const {return fTurnOffTime; }
    inline std::map<G4String, G4MoleculeDefinition*> GetNameMap() const {return fNameMap; }

    //-----------------Set Functions-----------------
    inline void SetVerbosity(G4int _value) { fVerbose = _value; }
    //-----------------------------------------------

    //-----------------------------------------------
    //                  Parameters
    //-----------------------------------------------
private:
    //-----------------Diffusion Parameters-----------------
    G4bool fIsSubDiffusion;
    G4double fDiffusionCoefficientForJump;
    G4double fDiffusionCoefficientForTrapped;
    G4double fMinWaitingTime;
    G4double fObserveDurationForMSD;
    G4double fObserveStepSizeForMSD;
    //-----------------Simulation Parameters-----------------
    TsParameterManager* fPm;
    G4String fName;
    G4int fAlternativeRunID;
    G4int fBiologyRepeatNumber;
    G4int fVerbose;
    G4double fBiologyEndTime;
    G4bool fTurnOffTime;
    G4String fExplicitBinning;
    std::vector<G4double> fExplicitBins;
    std::vector<G4double> fExplicitBinDiff;
    G4int fCurrentExplicitBinNuber;
    //-----------------Geometry Parameters-----------------
    G4double fBoundingRadius;
    //-----------------Damage Parameters-----------------
    G4int fSelectFromExposure;
    G4bool fIsReadFromFile;
    G4String fDamageFileName;
    G4double fDSBSeparation;
    G4double fDSBTimeDelay;
    std::vector<G4int> fDSBPattern;
    G4double fDSBOffset;
    G4int fDSBOriginNumber;
    G4int fDSBColumnNumber;
    G4double fDSBColumnRadius;
    //-----------------Biology Parameters-----------------
    std::vector<G4String> fHasPKcs;
    std::vector<G4String> fHasKu;
    std::vector<G4String> fHasXRCC4;
    std::map<G4String, G4MoleculeDefinition*> fNameMap;
    std::vector<std::pair<G4MoleculeDefinition*, DrProteinKinetics_Generic*> > fProcList;
};
