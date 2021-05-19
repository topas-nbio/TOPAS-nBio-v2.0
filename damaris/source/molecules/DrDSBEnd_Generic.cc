// Extra Class for use by DaMaRiS
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
// Created by John-William Warmenhoven on 27/10/18.
// DaMaRiS is developed at the University of Manchester.
// See README for references.
//

#include "DrDSBEnd_Generic.hh"
#include <G4PhysicalConstants.hh>
#include <G4SystemOfUnits.hh>
#include <G4ParticleTable.hh>

DrDSBEnd_Generic* DrDSBEnd_Generic::theInstance = nullptr;
std::map<G4String,DrDSBEnd_Generic*> DrDSBEnd_Generic::theInstanceMap;

DrDSBEnd_Generic* DrDSBEnd_Generic::Definition(const G4String name)
{
    G4double DNAMass    = ((10*607.4)+157.9)*1.6605e-27*kg* c_squared;
    G4double DNALength  = 10.0*0.332*nm;
    G4int DNAnAtoms  = int(10.0*63.5+0.5)+4;


    //@@@@ If the definition exist then just return the pointer
    //@@@@ to that definition, else create the definition and
    //@@@@ store it in the map
    if(theInstanceMap.find(name) != theInstanceMap.end()){
        //@@@@ dereference the static pointer to the map and
        //@@@@ then use name as key to return the pointer to
        //@@@@ the appropriate molecule defeinition(?)
        return theInstanceMap[name];
    }

    G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* anInstance = pTable->FindParticle(name);
    if (anInstance == nullptr){
        G4double diffCoef = 1.4 * (nm*nm/s);
        anInstance = new G4MoleculeDefinition(name, DNAMass, diffCoef, 0,-1, DNALength, DNAnAtoms, 0);
    }

    //@@@@ The definition has just been created and so needs to
    //@@@@ stored in the map.
    theInstanceMap[name] = reinterpret_cast<DrDSBEnd_Generic*>(anInstance);
    return theInstanceMap[name];
}