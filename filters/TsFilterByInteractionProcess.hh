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

#ifndef TsFilterByInteractionProcess_hh
#define TsFilterByInteractionProcess_hh

#include "TsVFilter.hh"

class G4ProcessVector;

class TsFilterByInteractionProcess : public TsVFilter
{
public:
	TsFilterByInteractionProcess(G4String name, TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM,
							 TsFilterManager* fM, TsVGenerator* generator, TsVScorer* scorer, TsVFilter* parentFilter);
	virtual ~TsFilterByInteractionProcess();

	void ResolveParameters();

	virtual G4bool Accept(const G4Step*) const;
	virtual G4bool AcceptTrack(const G4Track*) const;

private:
	G4ProcessVector* fProcesses;
};
#endif
