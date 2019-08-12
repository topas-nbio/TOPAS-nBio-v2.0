// Extra Class for TsEmDNAChemistryExtended 
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
 #include "TsO2.hh"
 #include "G4PhysicalConstants.hh"
 #include "G4SystemOfUnits.hh"
 #include "G4ParticleTable.hh"
  
 // ######################################################################
 // ###                         Peroxyde                               ###
 // ######################################################################
 TsO2* TsO2::theInstance = 0;
 
 TsO2* TsO2::Definition()
 {
   if (theInstance != 0) return theInstance;
 
   const G4String name = "TsO2";
 
   // search in particle table]
   G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
   G4ParticleDefinition* anInstance = pTable->FindParticle(name);
   if (anInstance == 0)
   {
     const G4String formatedName = "O_{2}";
 
     // create molecule
     //
     //      G4MoleculeDefinition(const G4String& name,
     //          G4double mass,
     //          G4double diffCoeff,
     //          G4int    charge = 0,
     //          G4int    electronicLevels = 0,
     //          G4double radius = -1,
     //          G4int    atomsNumber = -1,
     //          G4double lifetime = -1,
    //          G4String aType = "",
     //          G4FakeParticleID ID = G4FakeParticleID::Create()
    //      );
 
     G4double mass = 31.998 * g / Avogadro * c_squared;
     anInstance = new G4MoleculeDefinition(name, mass, 2.1e-9 * (m * m / s), 0, // charge
                                           2, // number of occupancies
                                           1.52 * angstrom, // radius
                                           2 // number of atoms
                                           );
 
     ((G4MoleculeDefinition*) anInstance)->SetLevelOccupation(0);
     ((G4MoleculeDefinition*) anInstance)->SetFormatedName(formatedName);
   }
   theInstance = reinterpret_cast<TsO2*>(anInstance);
   return theInstance;
 }

