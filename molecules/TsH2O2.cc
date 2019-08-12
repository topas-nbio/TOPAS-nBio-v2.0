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
 #include "TsH2O2.hh"
 #include "G4PhysicalConstants.hh"
 #include "G4SystemOfUnits.hh"
 #include "G4ParticleTable.hh"
  
 // ######################################################################
 // ###                         Peroxyde                               ###
 // ######################################################################
 TsH2O2* TsH2O2::theInstance = 0;
 
 TsH2O2* TsH2O2::Definition()
 {
   if (theInstance != 0) return theInstance;
 
   const G4String name = "TsH2O2";
 
   // search in particle table]
   G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
   G4ParticleDefinition* anInstance = pTable->FindParticle(name);
   if (anInstance == 0)
   {
     const G4String formatedName = "H_{2}O_{2}";
 
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
 
     G4double mass = 34.01468 * g / Avogadro * c_squared;
     anInstance = new G4MoleculeDefinition(name, mass, 1.4e-9 * (m * m / s), 0, // charge
                                           8, // number of occupancies
                                           3 * angstrom, // radius
                                            4 // number of atoms
                                           );
 
     ((G4MoleculeDefinition*) anInstance)->SetLevelOccupation(0);
     ((G4MoleculeDefinition*) anInstance)->SetLevelOccupation(1);
     ((G4MoleculeDefinition*) anInstance)->SetLevelOccupation(2);
     ((G4MoleculeDefinition*) anInstance)->SetLevelOccupation(3);
    ((G4MoleculeDefinition*) anInstance)->SetLevelOccupation(4);
     ((G4MoleculeDefinition*) anInstance)->SetLevelOccupation(5);
     ((G4MoleculeDefinition*) anInstance)->SetLevelOccupation(6);
     ((G4MoleculeDefinition*) anInstance)->SetLevelOccupation(7);
     ((G4MoleculeDefinition*) anInstance)->SetFormatedName(formatedName);
   }
   theInstance = reinterpret_cast<TsH2O2*>(anInstance);
   return theInstance;
 }

