// Extra Class for TsEmDNAChemistry
#include "TsIRT.hh"
#include "TsIRTUtils.hh"
#include "TsIRTConfiguration.hh"
#include "TsParameterManager.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4Molecule.hh"
#include "G4VSolid.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include "G4Timer.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4GeometryTolerance.hh"

#include <stdlib.h>
#include <algorithm>
#include <string>
#include <iostream>
#include <set>

TsIRT::TsIRT(TsParameterManager* pM, G4String parmName)
: fPm(pM), fName(parmName), fHighTimeScavenger(false), fVerbosity(0), fScorersInitialized(false)
{
	fReactionConf = new TsIRTConfiguration(parmName, fPm);
	fUtils = new TsIRTUtils();
	
	G4String chemistryList = fPm->GetStringParameter("Ch/ChemistryName");
	
	if ( fPm->ParameterExists("Ts/Verbosity"))
		fVerbosity = fPm->GetIntegerParameter("Ts/Verbosity");
	
	if ( fPm->ParameterExists(GetFullParmName("TimeLower")) ) {
		G4double tmin = fPm->GetDoubleParameter(GetFullParmName("TimeLower"), "Time");
		G4double tmax = fPm->GetDoubleParameter(GetFullParmName("TimeUpper"), "Time");
		G4int    tbin = fPm->GetIntegerParameter(GetFullParmName("TimeBins"));
		fStepTimes = fUtils->CreateTimeSteps(tmin, tmax, tbin, true);
		fReactionConf->SetTimeLimits(tmin, tmax);
		if ( fPm->ParameterExists(GetFullParmName("HighTimeScavenger"))) {
			fHighTimeScavenger = fPm->GetBooleanParameter(GetFullParmName("HighTimeScavenger"));
		}
		
	} else {
		fStepTimes = fUtils->CreateTimeSteps(0.1*ps, 1.0e6*ps, 100, true);
		fReactionConf->SetTimeLimits(0.1*ps, 1.0e6*ps);
	}
	
	fUseSpinScaled = false;
	if (fPm->ParameterExists(GetFullParmName("UseScaledProbabilityForSpinBehavior")))
		fUseSpinScaled = fPm->GetBooleanParameter(GetFullParmName("UseScaledProbabilityForSpinBehavior"));
	
	G4cout << "\n -------------- TOPAS IRT Warning -------------" << G4endl;
	
	if ( fUseSpinScaled )
		G4cout << " -- SpinBehavior is ScaledProbability" << G4endl;
	else
		G4cout << " -- SpinBehavior is ExplicitAssignament" << G4endl;
	
	G4cout << " ----------------------------------------------" << G4endl;
	G4cout << G4endl;
	
	
	fMoleculesIDs = fReactionConf->GetMoleculeIDs();
	std::map<G4int, G4String> moleculesNames = fReactionConf->GetMoleculeNames();
	
	if ( !fPm->ParameterExists(GetFullParmName("ReportMoleculesNamed")) ) {
		for ( int i = 1; i <= (int)moleculesNames.size(); i++ )
			fMolecules[i] = moleculesNames[i];
		
	} else {
		G4String* reportMolecules = fPm->GetStringVector(GetFullParmName("ReportMoleculesNamed"));
		G4int nbOfReportedMolecules = fPm->GetVectorLength(GetFullParmName("ReportMoleculesNamed"));
		
		for ( int i = 0; i < nbOfReportedMolecules; i++ )
			fMolecules[fMoleculesIDs[reportMolecules[i]]] = reportMolecules[i];
	}
	
	fReportDelta = false;
	if ( fPm->ParameterExists(GetFullParmName("ReportDeltaGValues"))) {
		fReportDelta = fPm->GetBooleanParameter(GetFullParmName("ReportDeltaGValues"));
	}
	
	fSortByTime = false;
	if ( fPm->ParameterExists(GetFullParmName("SortByTimeAscendant")))
		fSortByTime = fPm->GetBooleanParameter(GetFullParmName("SortByTimeAscendant"));
	
	fGlobalIndex = 0;
	
	fRCutOff = fReactionConf->GetRCutOff(1.0*us);
	fBinWidth = fRCutOff;
	if ( fVerbosity > 0 )
		G4cout << " -- Cut off -- " << fRCutOff/nm << " nm " << std::endl;
	
	if ( fPm->ParameterExists("Ch/"+chemistryList+"/SpaceBinningWidth") )
		fBinWidth = fPm->GetDoubleParameter("Ch/"+chemistryList+"/SpaceBinningWidth", "Length");
	
	fTestForContactReactions = false;
	if ( fPm->ParameterExists("Ch/"+chemistryList+"/TestForContactReactions"))
		fTestForContactReactions = fPm->GetBooleanParameter("Ch/"+chemistryList+"/TestForContactReactions");
	
	if ( fVerbosity > 0 && fTestForContactReactions )
		std::cout << " Going to test tracks by reaction at contact. " << G4endl;
	
	fTestIsInside = false;
	if (fPm->ParameterExists(GetFullParmName("OnlyIncludeChemicalSpeciesFromVirtualRegion/HLX"))) {
		fDx = fPm->GetDoubleParameter(GetFullParmName("OnlyIncludeChemicalSpeciesFromVirtualRegion/HLX"),"Length");
		fDy = fPm->GetDoubleParameter(GetFullParmName("OnlyIncludeChemicalSpeciesFromVirtualRegion/HLY"),"Length");
		fDz = fPm->GetDoubleParameter(GetFullParmName("OnlyIncludeChemicalSpeciesFromVirtualRegion/HLZ"),"Length");
		fTestIsInside = true;
	}
	
	fXMin = 1e9*nm;
	fYMin = 1e9*nm;
	fZMin = 1e9*nm;
	
	fXMax = 0e0*nm;
	fYMax = 0e0*nm;
	fZMax = 0e0*nm;
	
}


TsIRT::~TsIRT() {
	delete fReactionConf;
	delete fUtils;
}


G4bool TsIRT::Inside(G4ThreeVector p ) {
	G4double delta = 1 * nm;// * 1e-3 * nm;
	G4double dist = std::max(std::max(
									  std::abs(p.x())-fDx,
									  std::abs(p.y())-fDy),
							 std::abs(p.z())-fDz);
	if (dist > delta) return false;
	return !(dist > -delta);
}


void TsIRT::AddMolecule(TsIRTConfiguration::TsMolecule aMol) {
	if ( !fScorersInitialized )
		initializeScorers();
	
	G4ThreeVector position = aMol.position;
	
	if ( ! fSortByTime ) {
		fChemicalSpecies.push_back(aMol);
	} else {
		fVChemicalSpecies.push_back(aMol);
		fVTimes.push_back(std::make_pair(aMol.time, fGlobalIndex));
		fGlobalIndex++;
	}
	
	// count
	G4int tBin = fUtils->FindBin(aMol.time, fStepTimes);
	
	if ( -1 < tBin ) {
		for ( int tbin = tBin; tbin < (int)fStepTimes.size(); tbin++ ) {
			if ( fTheGvalue.find(aMol.id) == fTheGvalue.end() )
				fTheGvalue[aMol.id][tbin] = 1;
			else
				fTheGvalue[aMol.id][tbin]++;
		}
		
		if ( fTestIsInside ) {
			if ( Inside(aMol.position)) {
				for ( int tbin = tBin; tbin < (int)fStepTimes.size(); tbin++ ) {
					if ( fTheGvalueInVolume.find(aMol.id) == fTheGvalueInVolume.end() )
						fTheGvalueInVolume[aMol.id][tbin] = 1;
					else
						fTheGvalueInVolume[aMol.id][tbin]++;
				}
			}
		}
	}
	
	if ( fXMin > position.x() ) fXMin = position.x();
	if ( fYMin > position.y() ) fYMin = position.y();
	if ( fZMin > position.z() ) fZMin = position.z();
	
	if ( fXMax < position.x() ) fXMax = position.x();
	if ( fYMax < position.y() ) fYMax = position.y();
	if ( fZMax < position.z() ) fZMax = position.z();
}


void TsIRT::AddMolecule(G4Step* aStep, G4double time, G4int moleculeID, G4ThreeVector offset) {
	G4int pdg = -1;
	G4ThreeVector position = aStep->GetPreStepPoint()->GetPosition();
	const G4String& name = GetMolecule(aStep->GetTrack())->GetName();
	
	pdg = fMoleculesIDs[name];
	
	if ( pdg > 0 ) {
		TsIRTConfiguration::TsMolecule aMol;
		aMol.id = pdg;
		aMol.position = position + offset;
		aMol.time = time;
		aMol.reacted = false;
		// trackID != 0 reserved for DNA, < 0 = left strand, > 0 = right strand
		aMol.trackID = moleculeID;
		if ( 0 == moleculeID )
			aMol.isDNA = false;
		else
			aMol.isDNA = true;
		
		if ( pdg == 1 | pdg == 5 )
			aMol.spin = G4UniformRand() > 0.5 ? 1 : 0;
		else
			aMol.spin = -1;
		
		AddMolecule(aMol);
	}
}


void TsIRT::AddMolecule(G4Track* aTrack, G4double time, G4int moleculeID, G4ThreeVector offset) {
	G4int pdg = -1;
	G4ThreeVector position = aTrack->GetPosition();
	const G4String& name = GetMolecule(aTrack)->GetName();
	
	pdg = fMoleculesIDs[name];
	
	if ( pdg > 0 ) {
		TsIRTConfiguration::TsMolecule aMol;
		aMol.id = pdg;
		aMol.position = position + offset;
		aMol.time = time;
		aMol.reacted = false;
		// trackID != 0 reserved for DNA, < 0 = left strand, > 0 = right strand
		aMol.trackID = moleculeID;
		if ( 0 == moleculeID )
			aMol.isDNA = false;
		else
			aMol.isDNA = true;
		
		if ( pdg == 1 | pdg == 5 )
			aMol.spin = G4UniformRand() > 0.5 ? 1 : 0;
		else
			aMol.spin = -1;
		
		AddMolecule(aMol);
	}
}


void TsIRT::AddMolecule(const G4Track* aTrack, G4double time, G4int moleculeID, G4ThreeVector offset) {
	G4int pdg = -1;
	G4ThreeVector position = aTrack->GetPosition();
	const G4String& name = GetMolecule(aTrack)->GetName();
	
	pdg = fMoleculesIDs[name];
	
	if ( pdg > 0 ) {
		TsIRTConfiguration::TsMolecule aMol;
		aMol.id = pdg;
		aMol.position = position + offset;
		aMol.time = time;
		aMol.reacted = false;
		// trackID != 0 reserved for DNA, < 0 = left strand, > 0 = right strand
		aMol.trackID = moleculeID;
		if ( 0 == moleculeID )
			aMol.isDNA = false;
		else
			aMol.isDNA = true;
		
		if ( pdg == 1 | pdg == 5 )
			aMol.spin = G4UniformRand() > 0.5 ? 1 : 0;
		else
			aMol.spin = -1;
		
		AddMolecule(aMol);
	}
}


void TsIRT::Clean() {
	fChemicalSpecies.clear();
	
	fTheGvalue.clear();
	fTheGvalueInVolume.clear();
	fSpaceBinned.clear();
	fUsed.clear();
	
	fGValues.clear();
	fReactedDNA.clear();
	
	if ( fSortByTime ) {
		fVChemicalSpecies.clear();
		fVTimes.clear();
		fGlobalIndex = 0;
	}
	
	if (fReportDelta)
		fDeltaGValues.clear();
	
	if (fTestIsInside)
		fGValuesInVolume.clear();
	
	fScorersInitialized = false;
	
}


void TsIRT::FindBinIndexes(G4ThreeVector thisPos, G4double rcutOff) {
	fxiniIndex = fUtils->FindBin(fNx, fXMin, fXMax, thisPos.x()-rcutOff);
	fxendIndex = fUtils->FindBin(fNx, fXMin, fXMax, thisPos.x()+rcutOff);
	fyiniIndex = fUtils->FindBin(fNy, fYMin, fYMax, thisPos.y()-rcutOff);
	fyendIndex = fUtils->FindBin(fNy, fYMin, fYMax, thisPos.y()+rcutOff);
	fziniIndex = fUtils->FindBin(fNz, fZMin, fZMax, thisPos.z()-rcutOff);
	fzendIndex = fUtils->FindBin(fNz, fZMin, fZMax, thisPos.z()+rcutOff);
}


void TsIRT::initializeScorers() {
	for ( auto& nameAndID : fMolecules ) {
		G4int id = nameAndID.first;
		for ( int i = 0; i < (int)fStepTimes.size(); i++ ) {
			fTheGvalue[id][i] = 0.0;
		}
	}
	
	if ( fReportDelta ) {
		for  (int id = 0; id < fReactionConf->GetNumberOfReactions(); id++ ) {
			for (int i = 0; i < (int)fStepTimes.size(); i++ ) {
				G4double time = fStepTimes[i];
				fDeltaGValues[id][time] = 0.0;
			}
		}
	}
	
	if ( fTestIsInside ) {
		for ( auto& nameAndID : fMolecules ) {
			G4int id = nameAndID.first;
			for ( int i = 0; i < (int)fStepTimes.size(); i++ ) {
				fTheGvalueInVolume[id][i] = 0.0;
			}
		}
		std::cout << " Going to score G-value from species created/lost in box region of dimension: ("
		<< fDx*2/nm << ", " << fDy*2/nm << ", " << fDz*2/nm << ") nm "
		<< std::endl;
	}
	
	fScorersInitialized = true;
}


void TsIRT::contactReactions(G4int i) {
	fReactedByContact = false;
	
	// background reaction
	G4int indexOf1stOrder = fReactionConf->ContactFirstOrderAndBackgroundReactions(fChemicalSpecies[i]);
	if ( -1 < indexOf1stOrder ) {
		fReactedByContact = fReactionConf->MakeReaction(fChemicalSpecies, fSpaceBinned,
														fNx, fNy, fNz, fXMin, fXMax,
														fYMin, fYMax, fZMin, fZMax,
														fTheGvalue, fStepTimes, i,
														indexOf1stOrder, fChemicalSpecies[i].time, fUsed);
		if ( fReactedByContact ) {
			if ( fReportDelta ) {
				G4int tBin = fUtils->FindBin(fChemicalSpecies[i].time, fStepTimes);
				if ( -1 < tBin ) {
					for ( int ti = tBin; ti < (int)fStepTimes.size(); ti++ ) {
						fDeltaGValues[indexOf1stOrder][fStepTimes[ti]]++;
					}
				}
			}
			return;
		}
	} // background reaction
	
	G4ThreeVector thisPos = fChemicalSpecies[i].position;
	FindBinIndexes(thisPos, 0.0);
	
	for ( int ii = fxiniIndex; ii <= fxendIndex; ii++ ) {
		for ( int jj = fyiniIndex; jj <= fyendIndex; jj++ ) {
			for ( int kk = fziniIndex; kk <= fzendIndex; kk++ ) {
				for ( int n = 0; n < (int)fSpaceBinned[ii][jj][kk].size(); n++ ) {
					if ( fReactedByContact )
						continue;
					
					G4int j = fSpaceBinned[ii][jj][kk][n];
					if ( j == i )
						continue;
					
					if ( !fChemicalSpecies[j].reacted ) {
						G4int indexOfReaction = fReactionConf->GetReactionIndex(fChemicalSpecies[j].id,
																				fChemicalSpecies[i].id);
						if ( -1 < indexOfReaction ) {
							G4double r = fReactionConf->GetReaction(indexOfReaction).effectiveReactionRadius;
							if ( fReactionConf->GetReaction(indexOfReaction).reactionType == 4 )
								r = fReactionConf->GetReaction(indexOfReaction).effectiveTildeReactionRadius;
							
							G4double p = fReactionConf->GetReaction(indexOfReaction).probabilityOfReaction;
							G4double timeJ = fChemicalSpecies[j].time;
							
							if ( fChemicalSpecies[i].time != timeJ ) continue; /// Only zero timers<---Best match Clifford et al.
							
							if ( (fChemicalSpecies[i].position - fChemicalSpecies[j].position).mag2() < r*r ) {
								if ( !fTestIsInside )
									fReactedByContact = fReactionConf->MakeReaction(fChemicalSpecies, fSpaceBinned,
																					fNx, fNy, fNz, fXMin, fXMax,
																					fYMin, fYMax, fZMin, fZMax,
																					fTheGvalue, fStepTimes, i,
																					j, indexOfReaction, fChemicalSpecies[i].time, p, fUsed);
								else
									fReactedByContact = fReactionConf->MakeReaction(fChemicalSpecies, fSpaceBinned,
																					fNx, fNy, fNz, fXMin, fXMax,
																					fYMin, fYMax, fZMin, fZMax,
																					fTheGvalue, fTheGvalueInVolume, fStepTimes, i,
																					j, indexOfReaction, fChemicalSpecies[i].time, p, fUsed);
								
								if ( fReactedByContact ) {
									if ( fChemicalSpecies[j].isDNA )
										fReactedDNA[fChemicalSpecies[j].trackID] = std::make_pair(fChemicalSpecies[j].id,fChemicalSpecies[i].id);
									
									if ( fReportDelta ) {
										G4int tdBin = fUtils->FindBin(fChemicalSpecies[i].time, fStepTimes);
										if ( -1 < tdBin ) {
											for ( int ti = tdBin; ti < (int)fStepTimes.size(); ti++ ) {
												fDeltaGValues[indexOfReaction][fStepTimes[ti]]++;
											}
										}
									}
									return;
								}
							} // if p1 - p2 < r2
						} // -1 < indexOfReaction
					} // if reacted
				} // for n
			} // for ii
		} // for jj
	} // for ii
}


void TsIRT::runIRT() {
	fNx = G4int((fXMax-fXMin)/fBinWidth) == 0 ? 1 : G4int((fXMax-fXMin)/fBinWidth);
	fNy = G4int((fYMax-fYMin)/fBinWidth) == 0 ? 1 : G4int((fYMax-fYMin)/fBinWidth);
	fNz = G4int((fZMax-fZMin)/fBinWidth) == 0 ? 1 : G4int((fZMax-fZMin)/fBinWidth);
	
	if ( 1 < fVerbosity )
		std::cout << " Going to bin the space containing the track into "
		<< fNx << " x, " << fNy << " y and " << fNz << " z bins of width "
		<< fBinWidth/nm << " nm with limits (" << fXMin/nm << "," << fXMax/nm << "), ("
		<< fYMin/nm << ", " << fYMax/nm << "), (" << fZMin/nm << ", " << fZMax/nm << ") nm" << std::endl;
	
	if ( fSortByTime ) {
		std::sort(fVTimes.begin(), fVTimes.end());
		for ( size_t u = 0; u < fVChemicalSpecies.size(); u++ )
			fChemicalSpecies.push_back(fVChemicalSpecies[fVTimes[u].second]);
		
		fVChemicalSpecies.erase(fVChemicalSpecies.begin(), fVChemicalSpecies.end());
		fVTimes.erase(fVTimes.begin(), fVTimes.end());
		
		if ( fVerbosity > 0 )
			std::cout << " Going to pre-sort tracks by ascendant time " << G4endl;
	}
	
	for ( int t = 0; t < (int)fChemicalSpecies.size(); t++ ) {
		fUsed.push_back(false);
		G4int I = fUtils->FindBin(fNx, fXMin, fXMax, fChemicalSpecies[t].position.x());
		G4int J = fUtils->FindBin(fNy, fYMin, fYMax, fChemicalSpecies[t].position.y());
		G4int K = fUtils->FindBin(fNz, fZMin, fZMax, fChemicalSpecies[t].position.z());
		fSpaceBinned[I][J][K].push_back(t);
	}
	
	if ( 1 < fVerbosity )
		std::cout << " Processing " << fChemicalSpecies.size() << " species " << std::endl;
	
	if ( fTestForContactReactions ) { // Needed for QA
		
		for ( G4int i = 0; i < (G4int) fChemicalSpecies.size(); i++ ) {
			if ( fChemicalSpecies[i].reacted )
				continue;
			
			if ( fChemicalSpecies[i].isDNA )
				continue;
			
			contactReactions(i);
		}
		
		if ( 1 < fVerbosity )
			std::cout << "     -- After processing remained "
			<< fChemicalSpecies.size() << " species " << std::endl;
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 2.0 Sampling of the IRT
	G4double tCut = fStepTimes[fStepTimes.size()-1];
	std::vector<G4int> irtIndex;
	std::vector<std::pair<G4double, G4int > > irtValues;
	std::vector<G4int> irtMolA;
	std::vector<G4int> irtMolB;
	std::vector<G4bool> irtIsBackground;
	std::vector<G4double> irtOrigTime;
	std::vector<std::pair<G4ThreeVector, G4ThreeVector> > irtPositions;
	
	G4int index = 0;
	
	G4int chunk = 0;
	if ( 2 < fVerbosity ) {
		chunk = (G4int)fChemicalSpecies.size()/10;
		G4cout << "     -- Begin sampling of IRT for " << fChemicalSpecies.size() << " species " << G4endl;
	}
	
	for ( int i = 0; i < (G4int)fChemicalSpecies.size(); i++ ) {
		if ( fChemicalSpecies[i].reacted )
			continue;
		
		if ( fChemicalSpecies[i].isDNA )
			continue;
		
		G4ThreeVector thisPos = fChemicalSpecies[i].position;
		FindBinIndexes(thisPos, fRCutOff);
		
		fUsed[i] = true;
		for ( int ii = fxiniIndex; ii <= fxendIndex; ii++ ) {
			for ( int jj = fyiniIndex; jj <= fyendIndex; jj++ ) {
				for ( int kk = fziniIndex; kk <= fzendIndex; kk++ ) {
					for ( int n = 0; n < (int)fSpaceBinned[ii][jj][kk].size(); n++ ) {
						G4int j = fSpaceBinned[ii][jj][kk][n];
						
						if ( fChemicalSpecies[j].reacted || fUsed[j] )
							continue;
						
						G4int indexOfReaction = fReactionConf->GetReactionIndex(fChemicalSpecies[i].id,
																				fChemicalSpecies[j].id);
						if( -1 < indexOfReaction ) {
							G4double timeI = fChemicalSpecies[i].time;
							G4double timeJ = fChemicalSpecies[j].time;
							G4double dt = fChemicalSpecies[i].time - timeJ;
							//                          if ( dt < 0 ) continue; // Not afect if all initial times are the same, but affects for flash
							
							G4ThreeVector origPositionI = fChemicalSpecies[i].position;
							G4ThreeVector origPositionJ = fChemicalSpecies[j].position;
							G4bool resampleA = false, resampleB = false;
							if ( 0 < dt ) {
								fReactionConf->Diffuse(fChemicalSpecies[j],dt);
								resampleB = true;
							}
							
							if ( 0 > dt ) {
								fReactionConf->Diffuse(fChemicalSpecies[i],-dt);
								resampleA = true;
							}
							
							G4double irt = fReactionConf->GetIndependentReactionTime(fChemicalSpecies[i],
																					 fChemicalSpecies[j],
																					 indexOfReaction);
							
							if ( 0 < irt && irt <= tCut) {
								irt += fChemicalSpecies[i].time;
								irtValues.push_back(std::make_pair(irt,index));
								irtIndex.push_back(indexOfReaction);
								irtMolA.push_back(i);
								irtMolB.push_back(j);
								irtOrigTime.push_back(fChemicalSpecies[i].time);
								irtPositions.push_back(std::make_pair(fChemicalSpecies[i].position,
																	  fChemicalSpecies[j].position));
								irtIsBackground.push_back(false);
								index++;
								
							}
							
							if ( resampleB ) {
								fChemicalSpecies[j].position = origPositionJ;
								fChemicalSpecies[j].time = timeJ;
							}
							
							if ( resampleA ) {
								fChemicalSpecies[i].position = origPositionI;
								fChemicalSpecies[i].time = timeI;
							}
						}
					}
				}
			}
		}
		
		std::pair<G4int, G4double> tnIscav = fReactionConf->SampleIRTFirstOrderAndBackgroundReactions(fChemicalSpecies[i]);
		G4double tscav = tnIscav.second;
		
		if ((fHighTimeScavenger && 0 < tscav && (tscav <= tCut) ) ||
			(!fHighTimeScavenger && tscav >= fChemicalSpecies[i].time && (tscav <= tCut))) {
			tscav += fChemicalSpecies[i].time;
			irtValues.push_back(std::make_pair(tscav,index));
			irtIndex.push_back(tnIscav.first);
			irtMolA.push_back(i);
			irtMolB.push_back(i);
			irtOrigTime.push_back(fChemicalSpecies[i].time);
			irtPositions.push_back(std::make_pair(fChemicalSpecies[i].position, fChemicalSpecies[i].position));
			irtIsBackground.push_back(true);
			index++;
		}
		
		if ( 2 < fVerbosity && ( i >= chunk && i % chunk == 0 )) {
			G4cout << "     ---- current sampled times " << i << " out " << (G4int)fChemicalSpecies.size() <<  G4endl;
		}
	}
	
	// 4.0 re-arrange in ascendant order
	std::sort(irtValues.begin(), irtValues.end());
	
	if ( 2 < fVerbosity ) {
		chunk = (int)irtValues.size()/10;
		G4cout << "     -- Begin performing of reactions " << G4endl;
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////
	// 5.0 make reactions
	for ( int i = 0; i < (int)irtValues.size(); i++ ) {
		G4double irt = irtValues[i].first;
		
		if ( irt < tCut ) {
			G4int idx = irtValues[i].second;
			G4int iM = irtMolA[idx];
			G4int jM = irtMolB[idx];
			if ( fChemicalSpecies[iM].reacted || fChemicalSpecies[jM].reacted )
				continue;
			
			G4int indexOfReaction = irtIndex[idx];
			
			// 5.1 set the position of the reactants to those used to sample the irt
			if ( fReactionConf->GetReaction(indexOfReaction).reactionType == 5 ) {
				if (fUseSpinScaled) {
					if (fChemicalSpecies[iM].spin == fChemicalSpecies[jM].spin )
						continue;
					else if ( 0.5 > G4UniformRand() )
						continue;
				}
				else {
					if ( !(0 == fChemicalSpecies[iM].spin  && fChemicalSpecies[iM].spin == fChemicalSpecies[jM].spin  ))
						continue;
				}
			}
			
			G4int tBin = fUtils->FindBin(irt, fStepTimes);
			
			fChemicalSpecies[iM].position = irtPositions[idx].first;
			fChemicalSpecies[iM].time = irtOrigTime[idx];
			
			std::vector<G4ThreeVector> positions;
			std::vector<G4int> products;
			G4bool reactionIsInside = false;
			
			// 5.2 sample the position of the reactants using position approach Clifford et al 1986
			if ( !irtIsBackground[idx] ) {
				if ( !fChemicalSpecies[iM].isDNA && !fChemicalSpecies[jM].isDNA ) { // None is DNA
					fChemicalSpecies[jM].position = irtPositions[idx].second;
					fChemicalSpecies[jM].time = irtOrigTime[idx];
					
					fReactionConf->ResampleReactantsPosition(fChemicalSpecies[iM], fChemicalSpecies[jM],
															 indexOfReaction, irt);
					
					positions = fReactionConf->GetPositionOfProducts(fChemicalSpecies[iM],
																	 fChemicalSpecies[jM], indexOfReaction);
				} else { // at least one is DNA. Set product positions as the DNA molecule position.
					// Score the base pair ID, and molecule that caused the damage.
					if ( fChemicalSpecies[iM].isDNA ) {
						for ( int ip = 0; ip < 3; ip++ )
							positions.push_back(fChemicalSpecies[iM].position);
						
						fReactedDNA[fChemicalSpecies[iM].trackID] = std::make_pair(fChemicalSpecies[iM].id,fChemicalSpecies[jM].id);
					} else {
						for ( int ip = 0; ip < 3; ip++ )
							positions.push_back(fChemicalSpecies[jM].position);
						
						fReactedDNA[fChemicalSpecies[jM].trackID] = std::make_pair(fChemicalSpecies[jM].id,fChemicalSpecies[iM].id);
					}
				}
				
				products = (fReactionConf->GetReaction(indexOfReaction)).products;
				
				if ( 0 <= tBin ) {
					for ( int ti = tBin; ti < (int)fStepTimes.size(); ti++ ) {
						fTheGvalue[fChemicalSpecies[iM].id][ti]--;
						fTheGvalue[fChemicalSpecies[jM].id][ti]--;
						if ( fReportDelta )
							fDeltaGValues[indexOfReaction][fStepTimes[ti]]++;
						
					}
					
					if ( fTestIsInside ) {
						if ( Inside(fChemicalSpecies[iM].position)) {
							for ( int ti = tBin; ti < (int)fStepTimes.size(); ti++ ) {
								fTheGvalueInVolume[fChemicalSpecies[iM].id][ti]--;
								fTheGvalueInVolume[fChemicalSpecies[jM].id][ti]--;
							}
							reactionIsInside = true;
						}
					}
				}
				
				fChemicalSpecies[iM].reacted = true;
				fChemicalSpecies[jM].reacted = true;
				
			} else {
				fReactionConf->Diffuse(fChemicalSpecies[iM],irt-fChemicalSpecies[iM].time);
				for ( int ip = 0; ip < 3; ip++ )
					positions.push_back(fChemicalSpecies[iM].position);
				
				products = (fReactionConf->GetReaction(indexOfReaction)).products;
				if ( 0 <= tBin ) {
					for ( int ti = tBin; ti < (int)fStepTimes.size(); ti++ ) {
						fTheGvalue[fChemicalSpecies[iM].id][ti]--;
						if ( fReportDelta )
							fDeltaGValues[indexOfReaction][fStepTimes[ti]]++;
					}
					if ( fTestIsInside ) {
						if ( Inside(fChemicalSpecies[iM].position)) {
							for ( int ti = tBin; ti < (int)fStepTimes.size(); ti++ ) {
								fTheGvalueInVolume[fChemicalSpecies[iM].id][ti]--;
							}
							reactionIsInside = true;
						}
					}
				}
				
				fChemicalSpecies[iM].reacted = true;
			}
			
			// 5.3 create new products at positions resampled in previous step
			for ( size_t u = 0; u < products.size(); u++ ) {
				TsIRTConfiguration::TsMolecule aProd;
				aProd.id = products[u];
				aProd.position = positions[u];
				aProd.time = irt;
				aProd.reacted = false;
				aProd.trackID = 0;
				aProd.isDNA = false;
				
				if ( 1 == products[u] || 5 == products[u] )
					aProd.spin = G4UniformRand() > 0.5 ? 1 : 0;
				else
					aProd.spin = -1;
				
				fUsed.push_back(false);
				
				G4int newID = (G4int)fChemicalSpecies.size();
				
				G4int I = fUtils->FindBin(fNx, fXMin, fXMax, positions[u].x());
				G4int J = fUtils->FindBin(fNy, fYMin, fYMax, positions[u].y());
				G4int K = fUtils->FindBin(fNz, fZMin, fZMax, positions[u].z());
				
				fSpaceBinned[I][J][K].push_back(newID);
				
				// Merge this specie with the track
				fChemicalSpecies.push_back(aProd);
				
				if ( 0 <= tBin ) {
					for ( int ti = tBin; ti < (int)fStepTimes.size(); ti++ )
						fTheGvalue[aProd.id][ti]++;
					
					if ( fTestIsInside ) {
						if ( reactionIsInside) {
							for ( int ti = tBin; ti < (int)fStepTimes.size(); ti++ )
								fTheGvalueInVolume[aProd.id][ti]++;
						}
					}
				}
				
				contactReactions(newID);
				
				// then Sample reaction times of this product with existing initial species
				if ( !fReactedByContact ) {
					FindBinIndexes(positions[u], fRCutOff);
					////////////////////////////////
					for ( int ii = fxiniIndex; ii <= fxendIndex; ii++ ) {
						for ( int jj = fyiniIndex; jj <= fyendIndex; jj++ ) {
							for ( int kk = fziniIndex; kk <= fzendIndex; kk++ ) {
								for ( int n = 0; n < (int)fSpaceBinned[ii][jj][kk].size(); n++ ) {
									G4int j = fSpaceBinned[ii][jj][kk][n];
									if ( fChemicalSpecies[j].reacted )
										continue;
									
									if ( j == newID ) continue;
									G4int jIndexOfReaction = fReactionConf->GetReactionIndex(aProd.id,
																							 fChemicalSpecies[j].id);
									if( -1 < jIndexOfReaction ) {
										// temporary update position of j-Molecule for irt sampling
										G4ThreeVector origPositionI = aProd.position;
										G4ThreeVector origPositionJ = fChemicalSpecies[j].position;
										G4double timeI = aProd.time;
										G4double timeJ = fChemicalSpecies[j].time;
										G4bool resampleA = false, resampleB = false;
										G4double dt = timeI - timeJ;
										if ( 0 < dt ) {
											fReactionConf->Diffuse(fChemicalSpecies[j],dt);
											resampleB = true;
										}
										
										if ( 0 > dt ) { // same is done for initial species, see 442
											fReactionConf->Diffuse(aProd,-dt);
											resampleA = true;
										}
										
										G4double jIrt = fReactionConf->GetIndependentReactionTime(aProd,
																								  fChemicalSpecies[j],
																								  jIndexOfReaction);
										if ( 0 < jIrt && jIrt <= tCut) {
											jIrt += aProd.time;
											irtValues.push_back(std::make_pair(jIrt,index));
											irtIndex.push_back(jIndexOfReaction);
											irtMolA.push_back(newID);
											irtMolB.push_back(j);
											irtOrigTime.push_back(irt);
											irtPositions.push_back(std::make_pair(aProd.position, fChemicalSpecies[j].position));
											irtIsBackground.push_back(false);
											index++;
										}
										
										if ( resampleB ) {
											fChemicalSpecies[j].position = origPositionJ;
											fChemicalSpecies[j].time = timeJ;
										}
										
										if ( resampleA ) {
											aProd.position = origPositionI;
											aProd.time = timeI;
										}
									}
								}
							}
						}
					}
					// Check if there'll be background reaction
					std::pair<G4int, G4double> ptnIscav = fReactionConf->SampleIRTFirstOrderAndBackgroundReactions(aProd);
					G4double ptscav = ptnIscav.second;
					if ((fHighTimeScavenger && 0 < ptscav && (ptscav + aProd.time <= tCut)) ||
						(!fHighTimeScavenger && ptscav >= irt && (ptscav + aProd.time <= tCut))) {
						ptscav += aProd.time;
						irtValues.push_back(std::make_pair(ptscav,index));
						irtIndex.push_back(ptnIscav.first);
						irtMolA.push_back(newID);
						irtMolB.push_back(newID);
						irtOrigTime.push_back(irt);
						irtPositions.push_back(std::make_pair(aProd.position, aProd.position));
						irtIsBackground.push_back(true);
						index++;
					}
					//
				}
			}
			std::sort(irtValues.begin(), irtValues.end());
		} else {
			break;
		}
		
		if (fVerbosity > 2 && ( i >= chunk && i % chunk == 0 )) {
			G4cout << "     ---- current reactions performed " << i << " out " << chunk*10 << G4endl;
		}
	}
	
	for ( auto& nameAndID : fMolecules ) {
		G4String name = nameAndID.second;
		G4int id = nameAndID.first;
		for ( auto& timeAndGvalue : fTheGvalue[id] ) {
			G4int tBin = timeAndGvalue.first;
			G4double time = fStepTimes[tBin];
			fGValues[name][time] = timeAndGvalue.second;
		}
		if ( fTestIsInside ) {
			for ( auto& timeAndGvalue : fTheGvalueInVolume[id] ) {
				G4int tBin = timeAndGvalue.first;
				G4double time = fStepTimes[tBin];
				fGValuesInVolume[name][time] = timeAndGvalue.second;
			}
		}
	}
	
	irtIndex.clear();
	irtValues.clear();
	irtMolA.clear();
	irtMolB.clear();
	irtIsBackground.clear();
	irtOrigTime.clear();
	irtPositions.clear();
	fTheGvalue.clear();
	fTheGvalueInVolume.clear();
	fSpaceBinned.clear();
	fUsed.clear();
	
	fXMin = 1e9*nm;
	fYMin = 1e9*nm;
	fZMin = 1e9*nm;
	
	fXMax = 0e0*nm;
	fYMax = 0e0*nm;
	fZMax = 0e0*nm;
}


G4String TsIRT::GetFullParmName(G4String suffix ) {
	return "Sc/" + fName + "/" + suffix;
}


std::map<G4String, std::map<G4double, G4int>> TsIRT::GetGValues() {
	return fGValues;
}


std::map<G4String, std::map<G4double, G4int>> TsIRT::GetGValuesInVolume() {
	return fGValuesInVolume;
}


std::map<G4int, std::map<G4double, G4int>> TsIRT::GetDeltaGValues() {
	return fDeltaGValues;
}


std::map<G4int, std::pair<G4int,G4int>> TsIRT::GetReactedDNA() {
	return fReactedDNA;
}


std::pair<G4String, G4String> TsIRT::GetReactants(G4int ReactionIndex) {
	return fReactionConf->GetReactants(ReactionIndex);
}


std::vector<G4String> TsIRT::GetProducts(G4int ReactionIndex) {
	return fReactionConf->GetProducts(ReactionIndex);
}

