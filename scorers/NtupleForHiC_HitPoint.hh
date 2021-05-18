/*
*
*  Ingram & Henthorn et al
*  Hi-C TOPAS Geometry
*
*/

#ifndef NtupleForHiC_HitPoint_hh
#define NtupleForHiC_HitPoint_hh

#include <globals.hh>
#include <G4SystemOfUnits.hh>
#include "Randomize.hh"
#include <G4ThreeVector.hh>



using namespace std;

class NtupleForHiC_HitPoint{
public:
    NtupleForHiC_HitPoint(){};
        ~NtupleForHiC_HitPoint(){};
        
        void SetIsBack(G4bool opt){isBack=opt;}
        G4bool GetIsBack(){return isBack;}
        
        void SetIsBase(G4bool opt){isBase=opt;}
        G4bool GetIsBase(){return isBase;}
        
        void SetBP(G4int bp){BP=bp;}
        G4int GetBP(){return BP;}
        
        void SetClusterID(G4int val){ClusterID=val;}
        G4int GetClusterID(){return ClusterID;}
        
        void SetIsVisited(G4bool visit){Visited=visit;}
        G4bool GetIsVisited(){return Visited;}
        
        
        void SetEventID(G4int PID){EventID=PID;}
        G4int GetEventID(){return EventID;}
        
        void SetEnergy(G4double edep){EnergyeV=edep;}
        G4double GetEnergy(){return EnergyeV;}
        
        void SetStrand(G4int val){Strand=val;}
        G4int GetStrand(){return Strand;}
        
        void SetHitID(G4int val){HitID=val;}
        G4int GetHitID(){return HitID;}
        
        void SetChecked(G4bool opt){IsCheckedForDuplicates = opt;}
        G4bool GetChecked(){return IsCheckedForDuplicates;}
        
        void SetIsIonis(G4bool opt){IsIonisation = opt;}
        G4bool GetIsIonis(){return IsIonisation;}
        
        void SetPosition(G4ThreeVector pos){Position=pos;}
        G4ThreeVector GetPosition(){return Position;}
        
        void SetMarkedForDelete(G4bool opt){MarkedForDelete=opt;}
        G4bool GetMarkedForDelete(){return MarkedForDelete;}
        
        void SetIsChem(G4bool opt){isChem=opt;}
        G4bool GetIsChem(){return isChem;}
        
        void SetIsPhys(G4bool opt){isPhys=opt;}
        G4bool GetIsPhys(){return isPhys;}
    
        void SetIsChargeMig(G4bool opt){isChargeMig=opt;}
        G4bool GetIsChargeMig(){return isChargeMig;}
        
        void SetProcess(G4String proc){Process=proc;}
        G4String GetProcess(){return Process;}
        
        void SetChromID(G4int id){chromID=id;}
        G4int GetChromID(){return chromID;}
        
        void SetChromCopy(G4int val){chromCopy=val;}
        G4int GetChromCopy(){return chromCopy;}
        
        void SetChromatid(G4int val){chromatid=val;}
        G4int GetChromatid(){return chromatid;}
        
        void SetRunID(G4int RID){RunID=RID;}
        G4int GetRunID(){return RunID;}
        
        //New Elements - SDD V1.0
        void SetExposureTime(G4double TimeInMins){ExposureTime=TimeInMins;}
        G4double GetExposureTime(){return ExposureTime;}
        
        void SetLocalTime(G4double TimeInMins){LocalTime=TimeInMins;}
        G4double GetLocalTime(){return LocalTime;}
        
        void SetHitTimeInMinutes(G4double TimeInMinutes){HitTimeInMinutes=TimeInMinutes;}
        G4double GetHitTimeInMinutes(){return HitTimeInMinutes;}
        
        void SetParticleThatCausedDamage(G4int pa){ParticleThatCausedDamage=pa;}
        G4int GetParticleThatCausedDamage(){return ParticleThatCausedDamage;}
        void SetParticleThatCausedDamage_string(G4String name){ParticleThatCausedDamage_string=name;}
        G4String GetParticleThatCausedDamage_string(){return ParticleThatCausedDamage_string;}
        
        void SetEnergyOfParticleThatCausedDamage(G4double ene){EnergyOfParticleThatCausedDamage=ene;}
        G4double GetEnergyOfParticleThatCausedDamage(){return EnergyOfParticleThatCausedDamage;}
        
        void SetTranslationOfParticleThatCausedDamage(G4ThreeVector pa){TranslationOfParticleThatCausedDamage=pa;}
        G4ThreeVector GetTranslationOfParticleThatCausedDamage(){return TranslationOfParticleThatCausedDamage;}
        
        void SetDirectionOfParticleThatCausedDamage(G4ThreeVector pa){DirectionOfParticleThatCausedDamage=pa;}
        G4ThreeVector GetDirectionOfParticleThatCausedDamage(){return DirectionOfParticleThatCausedDamage;}
        
        void SetChromosomePosition(G4double po){chromosomeposition=po;}
        G4double GetChromosomePosition(){return chromosomeposition;}
        
        void SetDNAVolume(G4int vo /* 0) Missing, 1) A, 2) C, 3) T, 4) G*/){DNAVolume=vo;}
        G4int GetDNAVolume(){return DNAVolume;}
        
        void SetInHeterochromatin(G4bool opt){Heterochromatin=opt;}
        G4bool GetInHeterochromatin(){return Heterochromatin;}
        
        void SetInEuchromatin(G4bool opt){Euchromatin=opt;}
        G4bool GetInEuchromatin(){return Euchromatin;}
        
        void SetAcceptHit(G4bool opt){acceptHit=opt;}
        G4bool GetAcceptHit(){return acceptHit;}
        
        void SetVolume(G4String vol){volumeName=vol;}
        G4String GetVolume(){return volumeName;}
        
        void SetVolumeCopy(G4int i){VolumeCopy=i;}
        G4int GetVolumeCopy(){return VolumeCopy;}
        
        
    private:
        G4double EnergyeV=0.0;
        G4ThreeVector Position=G4ThreeVector(0,0,0);
        G4int Strand=0;
        
        G4int ClusterID=0;
        
        G4int EventID=0;
        G4int RunID=0;
        
        G4bool Visited=false;
        
        G4int HitID=0;
        
        G4bool isBack=false;
        G4bool isBase=false;
        G4int BP=0;
        
        G4bool IsIonisation=false;
        
        G4bool IsCheckedForDuplicates=false;
        
        G4bool MarkedForDelete=false;
        
        G4bool isPhys=false;
        G4bool isChem=false;
        G4bool isChargeMig=false;
        
        G4String Process;
        
        G4int chromID=-1;
        G4int chromCopy=-1;
        G4int chromatid=-1;
        
        G4double chromosomeposition=0.0;
        
        G4double HitTimeInMinutes=0.0;
        G4double ExposureTime=0.0;
        G4double LocalTime=0.0;
        
        G4int DNAVolume=0;
        
        G4int ParticleThatCausedDamage=0;
        G4String ParticleThatCausedDamage_string;
        G4double EnergyOfParticleThatCausedDamage=0.0;
        G4ThreeVector TranslationOfParticleThatCausedDamage;
        G4ThreeVector DirectionOfParticleThatCausedDamage;
        
        G4bool Heterochromatin=false;
        G4bool Euchromatin=false;
        
        G4bool acceptHit=false;
        
        G4String volumeName;
        
        G4int VolumeCopy=-1;
};


#endif 
