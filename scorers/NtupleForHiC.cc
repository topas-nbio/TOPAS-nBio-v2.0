// Scorer for NtupleForHiC
/*
*
*  Ingram & Henthorn et al
*  Hi-C TOPAS Geometry
*
*/


#include "NtupleForHiC.hh"



NtupleForHiC::NtupleForHiC(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM,
                                 TsScoringManager* scM, TsExtensionManager* eM,
                                 G4String scorerName, G4String quantity, G4String outFileName,
                                 G4bool isSubScorer)
: TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{
    //Define ntuple columns
    fNtuple->RegisterColumnD(&fPosX, "Pos_X", "nm");
    fNtuple->RegisterColumnD(&fPosY, "Pos_Y", "nm");
    fNtuple->RegisterColumnD(&fPosZ, "Pos_Z", "nm");
    fNtuple->RegisterColumnI(&fCause, "Cause");
    fNtuple->RegisterColumnI(&fStrand, "DNA_Strand");
    fNtuple->RegisterColumnI(&fChromosome, "ChromosomeID");
    fNtuple->RegisterColumnI(&fBeadID, "BeadID");
    fNtuple->RegisterColumnI(&fEventID, "Event_ID");
    fNtuple->RegisterColumnI(&fRunID, "Run_ID");
    
    
    //Set Some Scoring Parameters
    //How much of the volume is Sensitive
    if (fPm->ParameterExists("Sc/HiCScorer/SensitiveFraction")){
      SensitiveFraction = fPm->GetUnitlessParameter("Sc/HiCScorer/SensitiveFraction");
    }
    
    //Set the damage energy values
    if (fPm->ParameterExists("Sc/HiCScorer/MinEneRange")){
      EMin = fPm->GetDoubleParameter("Sc/HiCScorer/MinEneRange", "Energy") /eV;
    }
    if (fPm->ParameterExists("Sc/HiCScorer/MaxEneRange")){
      EMax = fPm->GetDoubleParameter("Sc/HiCScorer/MaxEneRange", "Energy") /eV;
    }
    
}


NtupleForHiC::~NtupleForHiC() {}


G4bool NtupleForHiC::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
    if (!fIsActive) {
        fSkippedWhileInactive++;
        return false;
    }

    ResolveSolid(aStep);

    G4double EneDep = aStep->GetTotalEnergyDeposit()/eV;
    G4String Vol = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();
    
    //check Edep is in a chromosome volume
    G4String SensVol="Chromosome";
    if (!(Vol.find(SensVol) != G4String::npos)){return false;}
    
    //check Edep is sufficient to make a break
    if (EneDep<EMin){return false;}
    G4double proba = (EneDep - EMin)/ (EMax-EMin);
    G4double rando1 = G4UniformRand();
    if (rando1>proba){return false;}
    
    //check Edep is "in" DNA
    G4double rando2 = G4UniformRand();
    if (rando2>=SensitiveFraction){
      return false;
    }
    
    //set damage cause (here, direct)
    //0=direct, 1=indirect, 2=combo, 3=charge migrations
    fCause=0;
    
    //get position
    G4ThreeVector pos = aStep->GetPreStepPoint()->GetPosition()/nm;
    fPosX=pos.x();
    fPosY=pos.y();
    fPosZ=pos.z();

    //set Run id's
    fEventID=GetEventID();
    fRunID=GetRunID();
    
    //set random DNA strand
    G4double rando3 = G4UniformRand();
    fStrand=1;
    if (rando3>0.5){fStrand=2;}
    
    //set chromosome ID
    fChromosome=stoi(Vol.substr(Vol.find("_") + 1));
    
    //set the bead copy number (used for fraction along chromosome)
//    G4int fBeadID_1=aStep->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo();
//    G4int fBeadID_2=aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber();
//    G4int fBeadID_3=aStep->GetTrack()->GetVolume()->GetCopyNo();
//    G4int fBeadID_4=aStep->GetTrack()->GetTouchable()->GetCopyNumber();
//    G4String fBeadID_5=Vol.substr(Vol.find("-") + 1);
    
    //TOPAS appends copy number to end of name, use "-" to identify
    fBeadID=stoi(Vol.substr(Vol.find("-") + 1));
    
//    std::cout<<Vol<<" "<<fChromosome<<" "<<fBeadID_1<<" "<<fBeadID_2<<" "<<fBeadID_3<<" "<<fBeadID_4<<" "<<fBeadID_5<<std::endl;
    
    //Put the units back to default
    fPosX*=nm;
    fPosY*=nm;
    fPosZ*=nm;

    //record hit
    fNtuple -> Fill();

    return true;

        
}


void NtupleForHiC::UserHookForEndOfEvent()
{
}

void NtupleForHiC::UserHookForEndOfRun()
{
    fNtuple->Write();

}
