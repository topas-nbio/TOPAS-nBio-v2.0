// Scorer for TsIRTInterPulse
//
// ********************************************************************
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * TOPAS collaboration.                                             *
// * Use or redistribution of this code is not permitted without the  *
// * explicit approval of the TOPAS collaboration.                    *
// * Contact: Joseph Perl, perl@slac.stanford.edu                     *
// *                                                                  *
// ********************************************************************
//

#include "TsScoreWithIRTCummulative.hh"
#include "TsIRT.hh"
#include "TsIRTConfiguration.hh"
#include "TsIRTUtils.hh"

#include "G4ITTrackHolder.hh"
#include "G4EventManager.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAChemistryManager.hh"
#include "G4Scheduler.hh"
#include "G4MolecularConfiguration.hh"
#include "G4Molecule.hh"
#include "G4H2O.hh"
#include "G4Electron_aq.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"
//#include <algorithm>
//#include <stdlib.h>
// #include <algorithm>

TsScoreWithIRTCummulative::TsScoreWithIRTCummulative(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                                                     G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
: TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer),
fPm(pM), fEnergyDepositPerEvent(0), fEnergyDepositPerEventEverywhere(0), fName(scorerName), fOldEvent(-1)
{
    SetUnit("");
    
    fIRT = new TsIRT(fPm, fName);
    fUtils = fIRT->GetUtils();
    fStepTimes = fIRT->GetStepTimes();
    for ( size_t u = 0; u < fStepTimes.size(); u++ ) {
        fVEnergyDeposit.push_back(0.0);
        fVEnergyDepositEverywhere.push_back(0.0);
    }
    
    fNtuple->RegisterColumnD(&fGValue, "TotalNumberOfSpeciesAtTime", "");
    fNtuple->RegisterColumnD(&fTime,    "Time", "picosecond");
    fNtuple->RegisterColumnS(&fMoleculeName, "MoleculeName");
    fNtuple->RegisterColumnD(&fEnergy, "EnergyDeposit", "eV");
    
    fPrescribedDose = fPm->GetDoubleParameter(GetFullParmName("PrescribedDose"),"Dose");

    fTimeDistribution = fPm->GetStringParameter(GetFullParmName("PulseDistribution"));
    fTimeDistribution.toLower();

    fNumberOfPulses = 1;
    if ( fPm->ParameterExists(GetFullParmName("NumberOfPulses")) ) {
        fNumberOfPulses = fPm->GetIntegerParameter(GetFullParmName("NumberOfPulses"));
        G4double pulsesFrequency = fPm->GetDoubleParameter(GetFullParmName("PulsesFrequency"),"perTime");
        if ( fNumberOfPulses > 1 && pulsesFrequency == 0 ) {
            G4cerr << "TOPAS is exiting due to an error in parameter" << G4endl;
            G4cerr << GetFullParmName("PulsesFrequency") << G4endl;
            G4cerr << "The number of pulses is bigger than 1, therefore, the frequency must be larger than zero" << G4endl;
        }
        fPulsesTimeDelay = 1./pulsesFrequency;
    }
    
    fLowLimitTo1ps = false;
    if ( fPm->ParameterExists(GetFullParmName("ForceLowTimeCutTo1ps")) )
        fLowLimitTo1ps = fPm->GetBooleanParameter(GetFullParmName("ForceLowTimeCutTo1ps"));
    
    if ( fTimeDistribution == "gaussian" ) {
        fTimeDistributionType = 1;
        fTimeMean = fPm->GetDoubleParameter(GetFullParmName("PulseTimeMean"),"Time");
        fTimeFWHM = fPm->GetDoubleParameter(GetFullParmName("PulseTimeFWHM"),"Time");
        fTimeStdv = fTimeFWHM/2.354820045;
    } else if ( fTimeDistribution == "uniform" ) {
        fTimeDistributionType = 2;
        fTimeMean = fPm->GetDoubleParameter(GetFullParmName("PulseTimeMean"),"Time");
        fTimeFWHM = fPm->GetDoubleParameter(GetFullParmName("PulseTimeFWHM"),"Time");
        fTimeStdv = fTimeFWHM/2.354820045;
    } else if ( fTimeDistribution == "none" ) {
        fTimeDistributionType = 0;
    } else if ( fTimeDistribution == "discrete" ) {
        fTimeDistributionType = 3;
        fTimeValues = fPm->GetDoubleVector(GetFullParmName("TimeValues"),"Time");
        fTimeWeights = fPm->GetUnitlessVector(GetFullParmName("TimeWeights"));
        fNbOfTimes = fPm->GetVectorLength(GetFullParmName("TimeValues"));
        G4double sum = 0.0;
        for ( int i = 0; i < fNbOfTimes; i++ )
            sum += fTimeWeights[i];
        
        if ( sum != 1.0 )
            for ( int i = 0; i < fNbOfTimes; i++ )
                fTimeWeights[i] /= sum;
        
        fTimeTops.push_back(fTimeWeights[0]);
        for ( int i = 1; i < fNbOfTimes; i++ )
            fTimeTops.push_back(fTimeTops[i-1] + fTimeWeights[i]);
        
    } else if ( fTimeDistribution == "exponential" ) {
        fTimeDistributionType = 4;
        fTimeMean = fPm->GetDoubleParameter(GetFullParmName("PulseTimeMean"),"Time");
    } else {
        G4cerr << "TOPAS is exiting due to an error in parameter" << G4endl;
        G4cerr << GetFullParmName("PulseDistribution") << G4endl;
        G4cerr << "Distribution: " << fTimeDistribution << " does not found" << G4endl;
        G4cerr << "Available distributions: Gaussian, Poisson and Uniform" << G4endl;
    }
    
    // To filter only physical interactiones within a virtual TsBox.
    fTestIsInside = false;
    if ( fPm->ParameterExists(GetFullParmName("SensitiveVolumeName"))) {
        fTestIsInside = true;
        fSensitiveVolume = fPm->GetStringParameter(GetFullParmName("SensitiveVolumeName"));
    }
    
    fTCut = 1.0 * ps;
    
    fNbOfScoredEvents = 0;
    fNbOfScoredEventsEverywhere = 0;

    fTotalDose = 0.0;
    fDosePerPulse = 0.0;
    fPulseTimeShift = 0.0;
    
    fShiftTime = 0.0 * ps;
    fMinShiftTime = 1.0 * s;
    
    G4String fileName = fPm->GetStringParameter("Sc/RootFileName") + ".bin";
    remove(fileName);
    fTimeOutFile.open(fileName, std::ios::binary | std::ios::app);
}


TsScoreWithIRTCummulative::~TsScoreWithIRTCummulative()
{
}


G4bool TsScoreWithIRTCummulative::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
    if (!fIsActive) {
        fSkippedWhileInactive++;
        return false;
    }
    
    if ( -1 < aStep->GetTrack()->GetTrackID() ) {
        G4double edep = aStep->GetTotalEnergyDeposit() ;
        if ( edep > 0 ) {
            ResolveSolid(aStep);
            edep *= aStep->GetPreStepPoint()->GetWeight();
            
            G4double mass = 0.0;
            G4double density = aStep->GetPreStepPoint()->GetMaterial()->GetDensity();
            fEnergyDepositPerEventEverywhere += edep ;

            if ( fTestIsInside ) {
                G4TouchableHistory* touchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
                const G4String& volumeName = touchable->GetVolume()->GetName();
                if( !volumeName.contains(fSensitiveVolume)) {
                    return false;
                }
            }
            
            mass = density * fSolid->GetCubicVolume();
            G4double dose = edep / mass;
            fTotalDose += dose;
            fEnergyDepositPerEvent += edep ;
            fDosePerPulse += dose;
            
            return true;
        }
    } else {
        
        fEventID = GetEventID();
        if ( fEventID != fOldEvent ) { // New history, sample a time within the pulse
            G4bool resample;
            if ( fTimeDistributionType == 1 ) {
                resample = true;
                while ( resample ) {
                    fShiftTime = G4RandGauss::shoot(fTimeMean+fPulseTimeShift, fTimeStdv);
                    if ( fShiftTime > 0 )
                        resample = false;
                }
            } else if ( fTimeDistributionType == 2 ) {
                resample = true;
                while (resample) {
                    fShiftTime = G4RandFlat::shoot(fTimeMean-fTimeFWHM*0.5, fTimeMean+fTimeFWHM*0.5);
                    if ( fShiftTime > 0 )
                        resample = false;
                }
                fShiftTime += fPulseTimeShift;
            } else if ( fTimeDistributionType == 3 ) {
                G4double aRandom = G4UniformRand();
                G4int j = fNbOfTimes - 1;
                while ((fTimeTops[j] >= aRandom) && (j >= 0)) {
                    j--;
                }
                fShiftTime = fTimeValues[j];
                fShiftTime += fPulseTimeShift;
            } else if ( fTimeDistributionType == 4 ) {
                resample = true;
                while (resample) {
                    fShiftTime = G4RandExponential::shoot(fTimeMean + fPulseTimeShift);
                    if ( fShiftTime >= 1.0*ps )
                        resample = false;
                }
            } else {
                fShiftTime = 1.0 * ps;
            }
            
            fOldEvent = fEventID;
            if ( fLowLimitTo1ps && fTimeDistributionType != 3 && (fNbOfScoredEvents == 0 && fNbOfScoredEventsEverywhere == 0))
                fShiftTime = 1.0 * ps;
            
            
        }
        
        return true;
    }
    return false;
}


void TsScoreWithIRTCummulative::UserHookForPreTimeStepAction() {
    if (G4Scheduler::Instance()->GetNbSteps() == 2) {
        G4TrackManyList* trackList = G4ITTrackHolder::Instance()->GetMainList();
        G4ManyFastLists<G4Track>::iterator it_begin = trackList->begin();
        G4ManyFastLists<G4Track>::iterator it_end = trackList->end();
        
        for(;it_begin!=it_end;++it_begin){
            if ( fTestIsInside ) {
                const G4String& volumeName = (*it_begin)->GetVolume()->GetName();
                if ( volumeName.contains(fSensitiveVolume) ) {
                    G4double time = fShiftTime;
                    fIRT->AddMolecule(*it_begin, time, 0, G4ThreeVector());
                }
            }
            else {
                fIRT->AddMolecule(*it_begin, fShiftTime, 0, G4ThreeVector());
            }
        }
        
        G4Scheduler::Instance()->Stop();
    }
}


void TsScoreWithIRTCummulative::UserHookForEndOfEvent() {
    if ( fEnergyDepositPerEvent > 0 ) {
        G4cout << "- current dose " << fTotalDose/gray << " Gy out "
        << fPrescribedDose/gray << " deposited within pulse at time "
        << fShiftTime/ps << " ps " << G4endl;
        
        fVEnergyDepositPerEvent.push_back(std::make_pair(fShiftTime,fEnergyDepositPerEvent));
        fNbOfScoredEvents++;
    }
    
    if ( fEnergyDepositPerEventEverywhere > 0 ) {
        fVEnergyDepositPerEventEverywhere.push_back(std::make_pair(fShiftTime,fEnergyDepositPerEventEverywhere));
        fNbOfScoredEventsEverywhere++;
    }
    
    fEnergyDepositPerEvent = 0.0;
    fEnergyDepositPerEventEverywhere = 0.0;
    
    if(fNumberOfPulses > 1 && fDosePerPulse >= fPrescribedDose/fNumberOfPulses) {
        fDosePerPulse = 0.0;
        fPulseTimeShift += fPulsesTimeDelay;
        G4cout << "-- New Pulse at " << fPulseTimeShift/ps << " ps " << G4endl;
    }
    
    if(fTotalDose > fPrescribedDose )  {
        G4int tBin;
        G4cout << " --- IRT start for event " << GetEventID() << G4endl;
        
        for ( size_t t = 0; t < fVEnergyDepositPerEvent.size(); t++ ) {
            if ( fTimeDistributionType != 3 )
                tBin = fUtils->FindBin(fVEnergyDepositPerEvent[t].first, fStepTimes);
            else
                tBin = fUtils->FindBin(fVEnergyDepositPerEvent[t].first, fStepTimes);
            
            for ( int i = tBin; i < (int)fStepTimes.size(); i++)
                fVEnergyDeposit[i] += fVEnergyDepositPerEvent[t].second;
            
            G4double saveTime = fVEnergyDepositPerEvent[t].first/ps;
            G4double saveEdep = fVEnergyDepositPerEvent[t].second/joule;
            fTimeOutFile.write(reinterpret_cast<char*>(&saveTime), sizeof saveTime);
            fTimeOutFile.write(reinterpret_cast<char*>(&saveEdep), sizeof saveEdep);
            
            G4cout << " --- energy deposit at time " << fVEnergyDepositPerEvent[t].first/ps
                   << " ps " << fVEnergyDepositPerEvent[t].second/eV << G4endl;
        }
        
        fIRT->runIRT();
        std::map<G4String, std::map<G4double, G4int>> irt = fIRT->GetGValues();
        
        G4cout << " --- IRT ends for event " << GetEventID() << G4endl;
        for ( auto& nameTimeAndGvalue : irt ) {
            G4String name = nameTimeAndGvalue.first;
            fMoleculeName = name;
            for ( auto& timeAndGvalue : (nameTimeAndGvalue.second) ) {
                G4double time = timeAndGvalue.first;
                G4double gvalue = timeAndGvalue.second;
                fGValue = gvalue;
                fTime = time;
                tBin = fUtils->FindBin(time, fStepTimes);
                fEnergy = fVEnergyDeposit[tBin];
                fNtuple->Fill();
            }
        }
        
        irt.clear();
        
        fIRT->Clean();
        
        G4cout << " ------ Needed " << fNbOfScoredEvents
        << " to achieve an accumulated dose of " << fTotalDose/gray
        << " Gy" << G4endl;
        
        fNtuple->Write();
        fTimeOutFile.close();
        Output();
        Clear();
        G4RunManager::GetRunManager()->AbortRun(true);
    }
}


