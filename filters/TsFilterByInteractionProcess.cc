// Filter for OnlyIncludeIfInteractionProcess
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

#include "TsFilterByInteractionProcess.hh"

#include "G4VProcess.hh"
#include "G4ProcessVector.hh"
#include "G4ProcessTable.hh"

TsFilterByInteractionProcess::TsFilterByInteractionProcess(G4String name, TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM,
												   TsFilterManager* fM, TsVGenerator* generator, TsVScorer* scorer, TsVFilter* parentFilter)
:TsVFilter(name, pM, mM, gM, fM, generator, scorer, parentFilter) {
	fProcesses = new G4ProcessVector();

	ResolveParameters();

}


TsFilterByInteractionProcess::~TsFilterByInteractionProcess()
{;}


void TsFilterByInteractionProcess::ResolveParameters() {
	fProcesses->clear();

	G4String* processNames = fPm->GetStringVector(GetFullParmName(GetName()));
	G4int length = fPm->GetVectorLength(GetFullParmName(GetName()));

	for (G4int i = 0; i < length; i++) {
		G4ProcessTable* theProcessTable = G4ProcessTable::GetProcessTable();
		G4ProcessVector* oneVector = theProcessTable->FindProcesses(processNames[i]);

		for (unsigned int j=0; j < oneVector->size(); j++)
			fProcesses->insert((*oneVector)[j]);

		if (fProcesses->size()==0) {
			G4cerr << "Topas is exiting due to a serious error in scoring setup." << G4endl;
			G4cerr << GetName() << " = " << processNames[i] << " refers to an unknown process." << G4endl;
			fPm->AbortSession(1);
		}
	}

	TsVFilter::ResolveParameters();
}


G4bool TsFilterByInteractionProcess::Accept(const G4Step* aStep) const {
	if (fParentFilter && !fParentFilter->Accept(aStep)) return false;

	for (unsigned int i=0; i < fProcesses->size(); i++)
		if (aStep->GetPostStepPoint()->GetProcessDefinedStep() == (*fProcesses)[i]) {
			if (fInvert) return false;
			else return true;
	    }

	if (fInvert) return true;
	else return false;
}


G4bool TsFilterByInteractionProcess::AcceptTrack(const G4Track*) const {
	G4cerr << "Topas is exiting due to a serious error in source setup." << G4endl;
	G4cerr << "Sources cannot be filtered by " << GetName() << G4endl;
	fPm->AbortSession(1);
	return false;
}
