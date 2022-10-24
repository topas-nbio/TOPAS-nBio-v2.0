// Scorer for TsIRTStrandBreaks
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

#include "TsScoreDNADamageWithIRT.hh"
#include "TsIRT.hh" 
#include "TsIRTConfiguration.hh"

#include "G4VSolid.hh"
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

TsScoreDNADamageWithIRT::TsScoreDNADamageWithIRT(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
												 G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
: TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer),
fPm(pM), fEnergyDepositPerEvent(0), fName(scorerName), fOutputFileName(outFileName)
{
	SetUnit("");
	
	fIRT = new TsIRT(fPm, fName);
	
	fNtuple->RegisterColumnD(&fTime, "Time", "ps");
	fNtuple->RegisterColumnI(&fOHDNAReactions, "Yield");
	fNtuple->RegisterColumnS(&fMoleculeName, "Molecule");
	fPrescribedDose = fPm->GetDoubleParameter(GetFullParmName("PrescribedDose"),"Dose");
	fNbOfScoredEvents = 0;
	fEnergyDepositPerEvent = 0.0;
}


TsScoreDNADamageWithIRT::~TsScoreDNADamageWithIRT()
{
	delete fIRT;
}


G4bool TsScoreDNADamageWithIRT::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
	if (!fIsActive) {
		fSkippedWhileInactive++;
		return false;
	}
	
	if ( -1 < aStep->GetTrack()->GetTrackID() ) {
		G4double edep = aStep->GetTotalEnergyDeposit();
		if ( edep > 0 ) {
			fEnergyDepositPerEvent += edep * aStep->GetPreStepPoint()->GetWeight();
			G4StepPoint* preStep = aStep->GetPreStepPoint();
			G4TouchableHistory* touchable = (G4TouchableHistory*)(preStep->GetTouchable());
			const G4String& vName = touchable->GetVolume()->GetName();
			
			
			if (vName.contains("deoxyribose" ) ) {
				G4int plasmidID = 0;
				G4int baseID = touchable->GetVolume(0)->GetCopyNo();
				if ( vName.contains("water") ) // get parent copy number
					plasmidID = touchable->GetVolume(1)->GetCopyNo();
				else // get grand parent copy number
					plasmidID = touchable->GetVolume(2)->GetCopyNo();
				
				G4bool found = false;
				if ( fDirectSSB.find(plasmidID) == fDirectSSB.end() )
					found = true;
				
				if (vName.contains("deoxyribose1") ) {
					if ( found ) {
						if (fDirectSSB[plasmidID].find(baseID) == fDirectSSB[plasmidID].end()) {
							fDirectSSB[plasmidID][baseID] = edep;
						} else {
							fDirectSSB[plasmidID][baseID] += edep;
						}
					} else {
						fDirectSSB[plasmidID][baseID] = edep;
					}
				} else if (vName.contains("deoxyribose2"))  {
					if ( found ) {
						if ( fDirectSSB[plasmidID].find(baseID) == fDirectSSB[plasmidID].end() ) {
							fDirectSSB[plasmidID][-baseID] = edep;
						} else {
							fDirectSSB[plasmidID][-baseID] += edep;
						}
					} else {
						fDirectSSB[plasmidID][-baseID] = edep;
					}
				}
			}
			return true;
		}
	} /*else {
		G4double time = aStep->GetTrack()->GetGlobalTime();
		if ( time >= 1.0*ps ) {
			time = 1.0 * ps;
			const G4String& matName = aStep->GetTrack()->GetMaterial()->GetName();
 			if ( matName == "G4_WATER" )
			    fIRT->AddMolecule(aStep, time, 0, G4ThreeVector());
			
			aStep->GetTrack()->SetTrackStatus(fStopAndKill);
			fNbOfScoredEvents++;
		}
		return true;
	}*/
	return false;
}

#include "G4LogicalVolumeStore.hh"

void TsScoreDNADamageWithIRT::UserHookForEndOfEvent() {
	
	G4double mass = 1.0 * g/cm3 * G4LogicalVolumeStore::GetInstance()->GetVolume(fComponent->GetName())->GetSolid()->GetCubicVolume();
	if( fEnergyDepositPerEvent/mass > fPrescribedDose)  {

		G4double totalDose = fEnergyDepositPerEvent/mass;
		
		G4cout << "-- Direct SSB [[[ ";
		G4int totalDirectSSB = 0;
		for ( auto& plasmidDirect : fDirectSSB ) {
			for ( auto& direct : plasmidDirect.second ) {
				if ( direct.second >= 17.5 * eV ) {
					G4cout << "  SB:" << direct.first << " " << direct.second << G4endl;
					totalDirectSSB++;
				}
			}
		}
		G4cout << "]]]" << G4endl;
		fDirectSSB.clear();
		
		G4cout << " --- IRT begins for total histories " << GetEventID() << G4endl;
		G4cout << " ---   Read dna information " << G4endl;
		G4String fileName = fPm->GetStringParameter(GetFullParmName("InputFile"));
		std::ifstream dna(fileName);
		std::vector<G4int> acceptedDNA;
		if ( fPm->ParameterExists(GetFullParmName("OnlyUseDNAMoleculesWithID")) ){
			G4int* accepted = fPm->GetIntegerVector(GetFullParmName("OnlyUseDNAMoleculesWithID"));
 			G4int naccepted = fPm->GetVectorLength(GetFullParmName("OnlyUseDNAMoleculesWithID"));
			for ( int k = 0; k < naccepted; k++ )
				acceptedDNA.push_back(accepted[k]);
		}
			
		G4int basePairID, moleculeID;
		G4double x, y, z;
		G4int totalSpecies = 0;
		while(1) {
			dna >> basePairID >> moleculeID >> x >> y >> z;
			if (!dna.good()) break;
			G4bool use = true; 
			if ( acceptedDNA.size() > 0 ) {
				if ( std::find(acceptedDNA.begin(), acceptedDNA.end(), moleculeID ) == acceptedDNA.end())
					use =false;
			}	
			if (!use) continue;
			x *= nm;
			y *= nm;
			z *= nm;
			TsIRTConfiguration::TsMolecule aMol;
			aMol.id = moleculeID;
			aMol.position = G4ThreeVector(x,y,z);
			aMol.time = 1.0*ps;
			aMol.reacted = false;
			aMol.trackID = basePairID;
			aMol.isDNA = true;
			aMol.spin = 0;
			fIRT->AddMolecule(aMol);
			totalSpecies++;
		}
		
		dna.close();
		G4cout << " ---   Running IRT for " << totalSpecies << " dna and "
		<< fNbOfScoredEvents << " chemical species" << G4endl;
	
		fIRT->runIRT();
		std::map<G4String, std::map<G4double, G4int>> irt = fIRT->GetGValues();
		std::map<G4int, std::pair<G4int,G4int>> breaks = fIRT->GetReactedDNA();
		G4int totalIndirectSSB = G4int(breaks.size());
		G4int totalSSB = totalDirectSSB + totalIndirectSSB;
		G4cout << " --- IRT ends " << GetEventID() << G4endl;

		G4int numberOfBasePairs = G4int(totalSpecies/4);
		G4double perDA = numberOfBasePairs * 660/0.34;
		G4double per100eV_to_umolPerJoule = 0.1036;

		G4cout << " -- total SSB: " << totalSSB << G4endl;
		G4cout << "=========== Total direct SSB        : " << totalDirectSSB << G4endl;
		G4cout << "=========== Total indirect SSB      : " << totalIndirectSSB << G4endl;
		G4cout << "=========== Total dose in envelope  : " << totalDose/gray<< " Gy " <<G4endl;
		G4cout << "=========== Total edep in envelope  : " << fEnergyDepositPerEvent/MeV << " MeV " << G4endl;
		G4cout << "=========== SSB in plasmid       /Gy: " <<  totalSSB/totalDose /gray       << " /Gy" << G4endl;
		G4cout << "=========== SSB in plasmid    /Gy/Da: " <<  totalSSB/totalDose/gray/perDA << " /Gy/Da" << G4endl;
		G4cout << "=========== SSB in plasmid    /100eV: " << 100*totalSSB/(fEnergyDepositPerEvent/eV) << " /100 eV" << G4endl;
		G4cout << "=========== SSB in plasmid    umol/J: " << per100eV_to_umolPerJoule*100*totalSSB/(fEnergyDepositPerEvent/eV) << " umol/J" << G4endl;
		G4cout << "=========== Total events to reach a dose of " << fPrescribedDose/gray << " (Gy) = " << GetEventID() << G4endl;


		std::ofstream out(fOutputFileName += ".dnadamaged" );
		out << "# | BasePair-ID | Strand-ID| " << G4endl;
		
		for ( auto& sb : breaks ) {
			//out << sb.first << "   " << sb.second.first << "   " << sb.second.second << G4endl;
			G4int strand = 1;
			if (sb.first < 0) {strand = 2;}
			out << std::setw(14) << abs(sb.first) << std::setw(13) << strand << G4endl;
		}
		out.close();
		
		for ( auto& nameTimeAndGvalue : irt ) {
			fMoleculeName = nameTimeAndGvalue.first;
			for ( auto& timeAndGvalue : (nameTimeAndGvalue.second) ) {
				fTime = timeAndGvalue.first;
				fOHDNAReactions = timeAndGvalue.second;
				fNtuple->Fill();
			}
		}
		
		fNbOfScoredEvents++;
		irt.clear();
		
		fEnergyDepositPerEvent = 0.0;
		fNtuple->Write();
		
		Output();
		Clear();
		G4RunManager::GetRunManager()->AbortRun(true);
	}
}


void TsScoreDNADamageWithIRT::UserHookForPreTimeStepAction() {
	if (G4Scheduler::Instance()->GetNbSteps() == 2) {
		G4TrackManyList* trackList = G4ITTrackHolder::Instance()->GetMainList();
		G4ManyFastLists<G4Track>::iterator it_begin = trackList->begin();
		G4ManyFastLists<G4Track>::iterator it_end = trackList->end();
		
		for(;it_begin!=it_end;++it_begin){
			const G4String& matName = it_begin->GetMaterial()->GetName();
			if ( matName == "G4_WATER" ) {
				G4double time = it_begin->GetGlobalTime();
				fIRT->AddMolecule(*it_begin, time, 0, G4ThreeVector());
			}
		}
		
		G4Scheduler::Instance()->Stop();
	}
}
