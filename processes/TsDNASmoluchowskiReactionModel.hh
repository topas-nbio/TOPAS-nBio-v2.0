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
// $Id: TsDNASmoluchowskiReactionModel.hh 100802 2016-11-02 14:55:27Z gcosmo $
//
// Author: Mathieu Karamitros
//
// The code is developed in the framework of the ESA AO7146
//
// We would be very happy hearing from you, send us your feedback! :)
//
// In order for Geant4-DNA to be maintained and still open-source,
// article citations are crucial. 
// If you use Geant4-DNA chemistry and you publish papers about your software, 
// in addition to the general paper on Geant4-DNA:
//
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we would be very happy if you could please also cite the following
// reference papers on chemistry:
//
// J. Comput. Phys. 274 (2014) 841-882
// Prog. Nucl. Sci. Tec. 2 (2011) 503-508 


#ifndef TsDNASmoluchowskiReactionModel_hh
#define TsDNASmoluchowskiReactionModel_hh

#include "AddClone_def.hh"
#include "G4VDNAReactionModel.hh"
#include <vector>

class G4DNAMolecularReactionData;

/**
  * TsDNASmoluchowskiReactionModel should be used
  * for very fast reactions (high reaction rate) : the reactions between
  * reactants occuring at encounter.
  * When the time step is constrained this model
  * uses brownian bridge : "Absorbing (Smoluchowski) boundary condition"
  * Reference : On the simulation of diffusion processes close to boundaries,
  * N. J. B. Green, Molecular Physics, 65: 6, 1399 — 1408(1988)
  */

class TsDNASmoluchowskiReactionModel : public G4VDNAReactionModel
{
public :
    TsDNASmoluchowskiReactionModel();
    virtual ~TsDNASmoluchowskiReactionModel();

    TsDNASmoluchowskiReactionModel(const TsDNASmoluchowskiReactionModel&);

    G4IT_ADD_CLONE(G4VDNAReactionModel, TsDNASmoluchowskiReactionModel)

    virtual void Initialise(G4MolecularConfiguration*, const G4Track&) ;
    virtual void InitialiseToPrint(G4MolecularConfiguration*) ;
    virtual G4double GetReactionRadius(G4MolecularConfiguration*,
                                       G4MolecularConfiguration*);
    virtual G4double GetReactionRadius(const G4int);

    virtual G4bool FindReaction(const G4Track&,
                                const G4Track&,
                                const G4double /*reactionRadius*/,
                                G4double& /*separationDistance*/,
                                const G4bool /*alongStepInteraction*/) ;

private :
    const std::vector<const G4DNAMolecularReactionData*>* fReactionData ;
    TsDNASmoluchowskiReactionModel& operator=(const TsDNASmoluchowskiReactionModel&);
};

#endif
