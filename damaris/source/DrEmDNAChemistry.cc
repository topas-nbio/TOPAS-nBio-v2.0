// Physics Module for DrEmDNAChemistry
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
// Created by John-William Warmenhoven.
// DaMaRiS is developed at the University of Manchester.
// See README for references.
//
#include "DrEmDNAChemistry.hh"
#include "TsParameterManager.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//JWW_DaMaRiS
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@ Models:
#include "DrDNAMolecularStepByStepModel.hh"
#include "DrDNASmoluchowskiReactionModel.hh"
//@@@@ Processes:
#include "DrFRAP.hh"
#include "DrTransportation_Subdiffusion_CTRW.hh"
#include "DrWaitingTime.hh"
#include "DrReportSystem.hh"
#include "DrReportMSD.hh"
//@@@@ Utils:
#include "DrDefinitions.hh"
//@@@@ For the clock molecule
#include "DrClock.hh"
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// Chemistry
#include "G4DNAWaterDissociationDisplacer.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAWaterExcitationStructure.hh"
#include "G4ProcessManager.hh"
#include "G4DNAGenericIonsManager.hh"

// Physics
#include "G4DNAElectronSolvation.hh"
#include "G4DNAAttachment.hh"
#include "G4DNAVibExcitation.hh"
#include "G4DNAElastic.hh"
#include "G4DNAChampionElasticModel.hh"
#include "G4DNAScreenedRutherfordElasticModel.hh"
#include "G4DNAUeharaScreenedRutherfordElasticModel.hh"
#include "G4DNAMolecularDissociation.hh"
#include "G4DNABrownianTransportation.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4DNAMolecularStepByStepModel.hh"
#include "G4VDNAReactionModel.hh"
#include "G4DNASmoluchowskiReactionModel.hh"
#include "G4DNAElectronHoleRecombination.hh"

// particles
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4GenericIon.hh"

// Molecules
#include "G4MoleculeTable.hh"
#include "G4H2O.hh"
#include "G4H2.hh"
#include "G4Hydrogen.hh"
#include "G4OH.hh"
#include "G4H3O.hh"
#include "G4Electron_aq.hh"
#include "G4H2O2.hh"


#include "G4PhysicsListHelper.hh"
#include "G4BuilderType.hh"

/****/
#include "G4DNAMoleculeEncounterStepper.hh"
#include "G4ProcessVector.hh"
#include "G4ProcessTable.hh"
#include "G4DNASecondOrderReaction.hh"
#include "G4MolecularConfiguration.hh"
/****/

// factory
#include "G4PhysicsConstructorFactory.hh"

G4_DECLARE_PHYSCONSTR_FACTORY(DrEmDNAChemistry);

#include "G4Threading.hh"

#include <map>

DrEmDNAChemistry::DrEmDNAChemistry()
        :G4VUserChemistryList()
{
    G4DNAChemistryManager::Instance()->SetChemistryList(this);
}

DrEmDNAChemistry::DrEmDNAChemistry(TsParameterManager* pM)
        :G4VUserChemistryList(), fPm(pM)
{
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    //@@@@ First call to set up instance and pass param manager
    DrDefinitions::Instance()->Initialise(fPm);
    //@@@@ Sets general parameters used in the simulation
    DrDefinitions::Instance()->InitialiseParameters();
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    DefineParameters();
    fPm->SetNeedsChemistry();
    G4DNAChemistryManager::Instance()->SetChemistryList(this);
}

DrEmDNAChemistry::~DrEmDNAChemistry() {}

void DrEmDNAChemistry::DefineParameters()
{
    fExistingMolecules["solvatedelectron"] = "e_aq";
    fExistingMolecules["hydroxyl"] = "OH";
    fExistingMolecules["hydrogen"] = "H";
    fExistingMolecules["hydronium"] = "H3Op";
    fExistingMolecules["dyhydrogen"] = "H2";
    fExistingMolecules["hydroxide"] = "OHm";
    fExistingMolecules["hydrogenperoxide"] = "H2O2";

    fDiffusionCoefficients["e_aq"] = 4.9e-9;
    fDiffusionCoefficients["OH"] = 2.8e-9;
    fDiffusionCoefficients["H"] = 7.0e-9;
    fDiffusionCoefficients["H3Op"] = 9.0e-9;
    fDiffusionCoefficients["H2"] = 4.8e-9;
    fDiffusionCoefficients["OHm"] = 5.0e-9;
    fDiffusionCoefficients["H2O2"] = 2.3e-9;

    fName = "Default";
    if ( fPm->ParameterExists("Ch/ChemistryName") )
        fName = fPm->GetStringParameter("Ch/ChemistryName");

    // Re-set diffusion coefficients
    if ( fPm->ParameterExists(GetFullParmName("DiffusionCoefficients/Molecules")) ) {
        G4String* molecules = fPm->GetStringVector(GetFullParmName("DiffusionCoefficients/Molecules"));
        G4double* diffusionCoefficients = fPm->GetUnitlessVector(GetFullParmName("DiffusionCoefficients/Values"));
        G4int nbOfMolecules = fPm->GetVectorLength(GetFullParmName("DiffusionCoefficients/Molecules"));

        if ( nbOfMolecules != fPm->GetVectorLength(GetFullParmName("DiffusionCoefficients/Values")) )
            Quit(GetFullParmName("DiffusionCoefficients/Molecules"),
                 "must match with " + GetFullParmName("DiffusionCoefficients/Values"));

        for ( int i = 0; i < nbOfMolecules; i++ ) {
            molecules[i].toLower();
            if ( !MoleculeExists(molecules[i]) ) {
                Quit(GetFullParmName("DiffusionCoefficients/Molecules"), ". Molecule: " + molecules[i] + " was not found in database");
            } else if ( molecules[i] == "water" ) {
                Quit(GetFullParmName("DiffusionCoefficients/Molecules"), ". Molecule: " + molecules[i] + " is not allowed to has a diffusion coefficient.");
            } else {
                std::map<G4String, G4String>::const_iterator it = fExistingMolecules.find(molecules[i]);
                fDiffusionCoefficients[it->second] = diffusionCoefficients[i];
            }
        }

    }

    // Re-set reaction rates
    std::vector<G4String>* reactionNames = new std::vector<G4String>;
    fPm->GetParameterNamesBracketedBy("Ch/" + fName, "ReactionRate", reactionNames);
    G4int numberOfReactions = reactionNames->size();
    G4int prefixLength = G4String("Ch/" + fName + "/BinaryReaction/").length();

    if ( numberOfReactions > 0 ) {
        for ( int i = 0; i < numberOfReactions; i++ ) {
            G4String parName = (*reactionNames)[i];
            G4String reactions = parName.substr(prefixLength, parName.find("ReactionRate")-prefixLength-1);
            G4String reactorA = reactions.substr(0, reactions.find("/"));
            G4String reactorB = reactions.substr(reactions.find("/") + 1);
            G4String* product = fPm->GetStringVector(parName.substr(0,parName.find("ReactionRate")-1) + "/Products");
            G4int nbOfProduct = fPm->GetVectorLength(parName.substr(0,parName.find("ReactionRate")-1) + "/Products");
            G4double reactionRate = fPm->GetUnitlessParameter(parName.substr(0,parName.find("ReactionRate")-1) + "/ReactionRate");

            std::vector<G4String> reactors;
            std::vector<G4String> products;

            reactorA.toLower();
            reactorB.toLower();
            if ( !MoleculeExists(reactorA) )
                Quit(parName, ". Molecule: " + reactorA + " was not found in database");

            if ( !MoleculeExists(reactorB) )
                Quit(parName, ". Molecule: " + reactorB + " was not found in database");

            reactors.push_back(fExistingMolecules[reactorA]);
            reactors.push_back(fExistingMolecules[reactorB]);

            for ( int i = 0; i < nbOfProduct; i++ ) {
                product[i].toLower();
                if ( product[i] == "noproduct" ) // Note on comparison of strings: See comment on ConstructReactionTable()
                    products.push_back("noproduct");
                else if ( product[i] == "water") {
                    products.push_back("H2O");
                } else {
                    if ( !MoleculeExists(product[i]) )
                        Quit(parName, ". Molecule: " + product[i] + " was not found in database");
                    else {
                        products.push_back(fExistingMolecules[product[i]]);
                    }
                }
            }

            fReactionSpecies.push_back(reactors);
            fReactionProducts.push_back(products);
            fReactionRates.push_back(reactionRate);
        }
    }

    fSetWaterConfiguration = false;
    if ( fPm->ParameterExists(GetFullParmName("SetWaterMolecularConfiguration")) )
        fSetWaterConfiguration = fPm->GetBooleanParameter(GetFullParmName("SetWaterMolecularConfiguration"));

    // Branching ratios
    G4String parName = GetFullParmName("BranchingRatios/");

    fIonizationStates = 1.0;
    if ( fPm->ParameterExists(parName + "IonizationStates/DissociativeDecay") )
        fIonizationStates = fPm->GetUnitlessParameter(parName + "IonizationStates/DissociativeDecay");
    std::cout << "-- Branching ratio: All ionization states probability " << fIonizationStates << std::endl;

    fA1B1DissociativeDecay = 0.65;
    if ( fPm->ParameterExists(parName + "A1B1/DissociativeDecay"))
        fA1B1DissociativeDecay = fPm->GetUnitlessParameter(parName + "A1B1/DissociativeDecay");
    std::cout << "-- Branching ratio: A1B1 Dissociative decay probability " << fA1B1DissociativeDecay << std::endl;

    fA1B1Relaxation = 0.35;
    if ( fPm->ParameterExists(parName + "A1B1/Relaxation"))
        fA1B1Relaxation = fPm->GetUnitlessParameter(parName + "A1B1/Relaxation");
    std::cout << "-- Branching ratio: A1B1 Relaxation probability " << fA1B1Relaxation << std::endl;

    fB1A1AutoIonization = 0.55;
    if ( fPm->ParameterExists(parName + "B1A1/AutoIonization"))
        fB1A1AutoIonization = fPm->GetUnitlessParameter(parName + "B1A1/AutoIonization");
    std::cout << "-- Branching ratio: B1A1 Auto-ionization probability " << fB1A1AutoIonization << std::endl;

    fB1A1DissociativeDecay = 0.15;
    if ( fPm->ParameterExists(parName + "B1A1/DissociativeDecay"))
        fB1A1DissociativeDecay = fPm->GetUnitlessParameter(parName + "B1A1/DissociativeDecay");
    std::cout << "-- Branching ratio: B1A1 Dissociative decay probability " << fB1A1DissociativeDecay << std::endl;

    fB1A1Relaxation = 0.30;
    if ( fPm->ParameterExists(parName + "B1A1/Relaxation"))
        fB1A1Relaxation = fPm->GetUnitlessParameter(parName + "B1A1/Relaxation");
    std::cout << "-- Branching ratio: B1A1 Relaxation probability " << fB1A1Relaxation << std::endl;

    fRydDiffAutoIonization = 0.5;
    if ( fPm->ParameterExists(parName + "RydbergStatesAndDiffuseBands/AutoIoinisation"))
        fRydDiffAutoIonization = fPm->GetUnitlessParameter(parName + "RydbergStatesAndDiffuseBands/AutoIoinisation");
    std::cout << "-- Branching ratio: Rydberg states and diffuse bands auto-ionisation probability " << fRydDiffAutoIonization << std::endl;

    fRydDiffRelaxation = 0.5;
    if ( fPm->ParameterExists(parName + "RydbergStatesAndDiffuseBands/Relaxation"))
        fRydDiffRelaxation = fPm->GetUnitlessParameter(parName + "RydbergStatesAndDiffuseBands/Relaxation");
    std::cout << "-- Branching ratio: Rydberg states and diffuse bands relaxation probability " << fRydDiffRelaxation << std::endl;

    fDissociativeAttachment = 1.0;
    if ( fPm->ParameterExists(parName + "DissociativeAttachment"))
        fDissociativeAttachment = fPm->GetUnitlessParameter(parName + "DissociativeAttachment");
    std::cout << "-- Branching ratio: dissociative attachment probability " << fDissociativeAttachment << std::endl;
}

G4bool DrEmDNAChemistry::MoleculeExists(G4String name) {
    if ( fExistingMolecules.find(name) == fExistingMolecules.end() )
        return false;
    return true;
}

void DrEmDNAChemistry::ConstructMolecule() {
    // Create the definition
    G4Electron::Definition(); // safety
    G4H2O::Definition();
    G4Hydrogen::Definition();
    G4H3O::Definition();
    G4OH::Definition();
    G4Electron_aq::Definition();
    G4H2O2::Definition();
    G4H2::Definition();

    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    DrClock::Definition();
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    //____________________________________________________________________________

    G4MoleculeTable::Instance()->CreateConfiguration("H3Op", G4H3O::Definition());
    G4MolecularConfiguration *OHm =
            G4MoleculeTable::Instance()->CreateConfiguration(
                    "OHm", // just a tag to store and retrieve from
                    // G4MoleculeTable
                    G4OH::Definition(),
                    -1, // charge
                    5.0e-9 * (m2 / s));
    OHm->SetMass(17.0079 * g / Avogadro * c_squared);
    G4MoleculeTable::Instance()->CreateConfiguration("OH", G4OH::Definition());
    G4MoleculeTable::Instance()->CreateConfiguration("e_aq",
                                                     G4Electron_aq::Definition());
    G4MoleculeTable::Instance()->CreateConfiguration("H",
                                                     G4Hydrogen::Definition());
    G4MoleculeTable::Instance()->CreateConfiguration("H2", G4H2::Definition());
    G4MoleculeTable::Instance()->CreateConfiguration("H2O2",
                                                     G4H2O2::Definition());

    if ( fSetWaterConfiguration )
        G4MoleculeTable::Instance()->CreateConfiguration("H2O", G4H2O::Definition());

    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    G4MoleculeTable::Instance()->CreateConfiguration("Clock", DrClock::Definition());
    DrDefinitions::Instance()->InitialisePathway();
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    // Re-set diffusion coefficients
    std::map<G4String, G4double>::iterator it = fDiffusionCoefficients.begin();
    while ( it != fDiffusionCoefficients.end() ) {
        G4String molecule = it->first;
        G4double diffusionCoefficient = it->second;
        G4MoleculeTable::Instance()->GetConfiguration(molecule)->SetDiffusionCoefficient(diffusionCoefficient * (m2/s));
        std::cout << "-- Diffusion coefficient for molecule: " << molecule
                  << " used the value " << diffusionCoefficient << " (m2/s) " << std::endl;
        it++;
    }
    std::cout << "" << std::endl;
}

void DrEmDNAChemistry::ConstructDissociationChannels() {
    //-----------------------------------
    // Get the molecular configuration
    G4MolecularConfiguration* OH = G4MoleculeTable::Instance()->GetConfiguration("OH");
    G4MolecularConfiguration* OHm = G4MoleculeTable::Instance()->GetConfiguration("OHm");
    G4MolecularConfiguration* e_aq = G4MoleculeTable::Instance()->GetConfiguration("e_aq");
    G4MolecularConfiguration* H2 = G4MoleculeTable::Instance()->GetConfiguration("H2");
    G4MolecularConfiguration* H3O = G4MoleculeTable::Instance()->GetConfiguration("H3Op");
    G4MolecularConfiguration* H = G4MoleculeTable::Instance()->GetConfiguration("H");

    //-------------------------------------
    //Define the decay channels
    G4MoleculeDefinition* water = G4H2O::Definition();
    G4MolecularDissociationChannel* decCh1;
    G4MolecularDissociationChannel* decCh2;

    G4ElectronOccupancy* occ = new G4ElectronOccupancy(*(water->GetGroundStateElectronOccupancy()));

    //////////////////////////////////////////////////////////
    //            EXCITATIONS                               //
    //////////////////////////////////////////////////////////
    G4DNAWaterExcitationStructure waterExcitation;
    //--------------------------------------------------------
    //---------------Excitation on the fifth layer------------

    decCh1 = new G4MolecularDissociationChannel("A^1B_1_Relaxation");
    decCh2 = new G4MolecularDissociationChannel("A^1B_1_DissociativeDecay");
    //Decay 1 : OH + H
    decCh1->SetEnergy(waterExcitation.ExcitationEnergy(0));
    decCh1->SetProbability(fA1B1Relaxation);
    decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::NoDisplacement);

    decCh2->AddProduct(OH);
    decCh2->AddProduct(H);
    decCh2->SetProbability(fA1B1DissociativeDecay);
    decCh2->SetDisplacementType(G4DNAWaterDissociationDisplacer::A1B1_DissociationDecay);

    //  water->AddExcitedState("A^1B_1");
    occ->RemoveElectron(4, 1); // this is the transition form ground state to
    occ->AddElectron(5, 1); // the first unoccupied orbital: A^1B_1

    water->NewConfigurationWithElectronOccupancy("A^1B_1", *occ);
    water->AddDecayChannel("A^1B_1", decCh1);
    water->AddDecayChannel("A^1B_1", decCh2);

    //--------------------------------------------------------
    //---------------Excitation on the fourth layer-----------
    decCh1 = new G4MolecularDissociationChannel("B^1A_1_Relaxation_Channel");
    decCh2 = new G4MolecularDissociationChannel("B^1A_1_DissociativeDecay");
    G4MolecularDissociationChannel* decCh3 = new G4MolecularDissociationChannel("B^1A_1_AutoIonisation_Channel");

    //Decay 1 : energy
    decCh1->SetEnergy(waterExcitation.ExcitationEnergy(1));
    decCh1->SetProbability(fB1A1Relaxation);

    //Decay 2 : 2OH + H_2
    decCh2->AddProduct(H2);
    decCh2->AddProduct(OH);
    decCh2->AddProduct(OH);
    decCh2->SetProbability(fB1A1DissociativeDecay);
    decCh2->SetDisplacementType(G4DNAWaterDissociationDisplacer::B1A1_DissociationDecay);

    //Decay 3 : OH + H_3Op + e_aq
    decCh3->AddProduct(OH);
    decCh3->AddProduct(H3O);
    decCh3->AddProduct(e_aq);
    decCh3->SetProbability(fB1A1AutoIonization);
    decCh3->SetDisplacementType(G4DNAWaterDissociationDisplacer::AutoIonisation);

    *occ = *(water->GetGroundStateElectronOccupancy());
    occ->RemoveElectron(3); // this is the transition form ground state to
    occ->AddElectron(5, 1); // the first unoccupied orbital: B^1A_1

    water->NewConfigurationWithElectronOccupancy("B^1A_1", *occ);
    water->AddDecayChannel("B^1A_1", decCh1);
    water->AddDecayChannel("B^1A_1", decCh2);
    water->AddDecayChannel("B^1A_1", decCh3);

    //-------------------------------------------------------
    //-------------------Excitation of 3rd layer-----------------
    decCh1 = new G4MolecularDissociationChannel("Excitation3rdLayer_AutoIonisation_Channel");
    decCh2 = new G4MolecularDissociationChannel("Excitation3rdLayer_Relaxation_Channel");

    //Decay channel 1 : : OH + H_3Op + e_aq
    decCh1->AddProduct(OH);
    decCh1->AddProduct(H3O);
    decCh1->AddProduct(e_aq);

    decCh1->SetProbability(fRydDiffAutoIonization);
    decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::AutoIonisation);

    //Decay channel 2 : energy
    decCh2->SetEnergy(waterExcitation.ExcitationEnergy(2));
    decCh2->SetProbability(fRydDiffRelaxation);

    //Electronic configuration of this decay
    *occ = *(water->GetGroundStateElectronOccupancy());
    occ->RemoveElectron(2, 1);
    occ->AddElectron(5, 1);

    //Configure the water molecule
    water->NewConfigurationWithElectronOccupancy("Excitation3rdLayer", *occ);
    water->AddDecayChannel("Excitation3rdLayer", decCh1);
    water->AddDecayChannel("Excitation3rdLayer", decCh2);

    //-------------------------------------------------------
    //-------------------Excitation of 2nd layer-----------------
    decCh1 = new G4MolecularDissociationChannel("Excitation2ndLayer_AutoIonisation_Channel");
    decCh2 = new G4MolecularDissociationChannel("Excitation2ndLayer_Relaxation_Channel");

    //Decay Channel 1 : : OH + H_3Op + e_aq
    decCh1->AddProduct(OH);
    decCh1->AddProduct(H3O);
    decCh1->AddProduct(e_aq);

    decCh1->SetProbability(fRydDiffAutoIonization);
    decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::AutoIonisation);

    //Decay channel 2 : energy
    decCh2->SetEnergy(waterExcitation.ExcitationEnergy(3));
    decCh2->SetProbability(fRydDiffRelaxation);

    *occ = *(water->GetGroundStateElectronOccupancy());
    occ->RemoveElectron(1, 1);
    occ->AddElectron(5, 1);

    water->NewConfigurationWithElectronOccupancy("Excitation2ndLayer", *occ);
    water->AddDecayChannel("Excitation2ndLayer", decCh1);
    water->AddDecayChannel("Excitation2ndLayer", decCh2);

    //-------------------------------------------------------
    //-------------------Excitation of 1st layer-----------------
    decCh1 = new G4MolecularDissociationChannel("Excitation1stLayer_AutoIonisation_Channel");
    decCh2 = new G4MolecularDissociationChannel("Excitation1stLayer_Relaxation_Channel");

    *occ = *(water->GetGroundStateElectronOccupancy());
    occ->RemoveElectron(0, 1);
    occ->AddElectron(5, 1);

    //Decay Channel 1 : : OH + H_3Op + e_aq
    decCh1->AddProduct(OH);
    decCh1->AddProduct(H3O);
    decCh1->AddProduct(e_aq);
    decCh1->SetProbability(fRydDiffAutoIonization);
    decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::AutoIonisation);

    //Decay channel 2 : energy
    decCh2->SetEnergy(waterExcitation.ExcitationEnergy(4));
    decCh2->SetProbability(fRydDiffRelaxation);

    water->NewConfigurationWithElectronOccupancy("Excitation1stLayer", *occ);
    water->AddDecayChannel("Excitation1stLayer", decCh1);
    water->AddDecayChannel("Excitation1stLayer", decCh2);

    /////////////////////////////////////////////////////////
    //                  IONISATION                         //
    /////////////////////////////////////////////////////////
    //--------------------------------------------------------
    //------------------- Ionisation -------------------------

    decCh1 = new G4MolecularDissociationChannel("Ionisation_Channel");

    //Decay Channel 1 : : OH + H_3Op
    decCh1->AddProduct(H3O);
    decCh1->AddProduct(OH);
    decCh1->SetProbability(fIonizationStates);
    decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::Ionisation_DissociationDecay);

    *occ = *(water->GetGroundStateElectronOccupancy());
    occ->RemoveElectron(4, 1);
    // this is a ionized h2O with a hole in its last orbital
    water->NewConfigurationWithElectronOccupancy("Ionisation5", *occ);
    water->AddDecayChannel("Ionisation5", decCh1);

    *occ = *(water->GetGroundStateElectronOccupancy());
    occ->RemoveElectron(3, 1);
    water->NewConfigurationWithElectronOccupancy("Ionisation4", *occ);
    water->AddDecayChannel("Ionisation4", new G4MolecularDissociationChannel(*decCh1));

    *occ = *(water->GetGroundStateElectronOccupancy());
    occ->RemoveElectron(2, 1);
    water->NewConfigurationWithElectronOccupancy("Ionisation3", *occ);
    water->AddDecayChannel("Ionisation3", new G4MolecularDissociationChannel(*decCh1));

    *occ = *(water->GetGroundStateElectronOccupancy());
    occ->RemoveElectron(1, 1);
    water->NewConfigurationWithElectronOccupancy("Ionisation2", *occ);
    water->AddDecayChannel("Ionisation2", new G4MolecularDissociationChannel(*decCh1));

    *occ = *(water->GetGroundStateElectronOccupancy());
    occ->RemoveElectron(0, 1);
    water->NewConfigurationWithElectronOccupancy("Ionisation1", *occ);
    water->AddDecayChannel("Ionisation1", new G4MolecularDissociationChannel(*decCh1));

    //////////////////////////////////////////////////////////
    //            Dissociative Attachment                   //
    //////////////////////////////////////////////////////////
    decCh1 = new G4MolecularDissociationChannel("DissociativeAttachment");

    //Decay 1 : 2OH + H_2
    decCh1->AddProduct(H2);
    decCh1->AddProduct(OHm);
    decCh1->AddProduct(OH);
    decCh1->SetProbability(fDissociativeAttachment);
    decCh1->SetDisplacementType(G4DNAWaterDissociationDisplacer::DissociativeAttachment);

    *occ = *(water->GetGroundStateElectronOccupancy());
    occ->AddElectron(5, 1); // H_2O^-
    water->NewConfigurationWithElectronOccupancy("DissociativeAttachment", *occ);
    water->AddDecayChannel("DissociativeAttachment", decCh1);

    delete occ;
}

void DrEmDNAChemistry::ConstructReactionTable(
        G4DNAMolecularReactionTable *theReactionTable) {
    //-----------------------------------
    //Get the molecular configuration
    G4MolecularConfiguration* OH = G4MoleculeTable::Instance()->GetConfiguration("OH");
    G4MolecularConfiguration* OHm = G4MoleculeTable::Instance()->GetConfiguration("OHm");
    G4MolecularConfiguration* e_aq = G4MoleculeTable::Instance()->GetConfiguration("e_aq");
    G4MolecularConfiguration* H2 = G4MoleculeTable::Instance()->GetConfiguration("H2");
    G4MolecularConfiguration* H3Op = G4MoleculeTable::Instance()->GetConfiguration("H3Op");
    G4MolecularConfiguration* H = G4MoleculeTable::Instance()->GetConfiguration("H");
    G4MolecularConfiguration* H2O2 = G4MoleculeTable::Instance()->GetConfiguration("H2O2");

    std::map<G4String, G4MolecularConfiguration*> reactions;
    reactions["OH"] = OH;
    reactions["OHm"] = OHm;
    reactions["e_aq"] = e_aq;
    reactions["H2"] = H2;
    reactions["H3Op"] = H3Op;
    reactions["H"] = H;
    reactions["H2O2"] = H2O2;

    if ( fSetWaterConfiguration ) {
        G4MolecularConfiguration* H2O = G4MoleculeTable::Instance()->GetConfiguration("H2O");
        reactions["H2O"] = H2O;
    }

    G4DNAMolecularReactionData* reactionData;
    for ( size_t t = 0; t < fReactionSpecies.size(); t++ ) {

        reactionData = new G4DNAMolecularReactionData(fReactionRates[t] * (1e-3 * m3 / (mole * s)), reactions[fReactionSpecies[t][0]], reactions[fReactionSpecies[t][1]]);

        for ( size_t u = 0; u < fReactionProducts[t].size(); u++ ) {
            if ( "noproduct" != fReactionProducts[t][u] ) // This comparison crashes if the order is fReactionProducts[t][u] != "noproduct"
                reactionData->AddProduct(reactions[ fReactionProducts[t][u] ] );
        }

        theReactionTable->SetReaction(reactionData);
    }

    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    std::vector<std::pair<G4double, std::vector<G4String>>> tempStore = DrDefinitions::Instance()->InitialiseReactions();

    for (auto reaction : tempStore) {
        G4MolecularConfiguration *reactant1 = G4MoleculeTable::Instance()->GetConfiguration(reaction.second[0]);
        G4MolecularConfiguration *reactant2 = G4MoleculeTable::Instance()->GetConfiguration(reaction.second[1]);
        G4MolecularConfiguration *product = G4MoleculeTable::Instance()->GetConfiguration(reaction.second[2]);
        G4double reactionRadius = reaction.first;

        G4double reactionRate{0.0};
        if (reactant1->GetName() != reactant2->GetName()) {
            reactionRate = 2 * 4 * pi * Avogadro * reactionRadius * reactant1->GetDiffusionCoefficient();
        } else {
            reactionRate = 4 * pi * Avogadro * reactionRadius * reactant1->GetDiffusionCoefficient();
        }
        reactionData = new G4DNAMolecularReactionData(reactionRate, reactant1, reactant2);
        reactionData->AddProduct(product);
        theReactionTable->SetReaction(reactionData);
    }
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
}

void DrEmDNAChemistry::ConstructProcess() {
    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

    //===============================================================
    // Extend vibrational to low energy
    // Anyway, solvation of electrons is taken into account from 7.4 eV
    // So below this threshold, for now, no accurate modeling is done
    //
    G4VProcess* process =
            G4ProcessTable::GetProcessTable()->
                    FindProcess("e-_G4DNAVibExcitation", "e-");

    if (process)
    {
        G4DNAVibExcitation* vibExcitation = (G4DNAVibExcitation*) process;
        G4VEmModel* model = vibExcitation->EmModel();
        G4DNASancheExcitationModel* sancheExcitationMod =
                dynamic_cast<G4DNASancheExcitationModel*>(model);
        if(sancheExcitationMod)
        {
            sancheExcitationMod->ExtendLowEnergyLimit(0.025 * eV);
        }
    }

    //===============================================================
    // Define processes for molecules
    //
    G4MoleculeTable *theMoleculeTable = G4MoleculeTable::Instance();
    G4MoleculeDefinitionIterator iterator =
            theMoleculeTable->GetDefintionIterator();
    iterator.reset();

    while (iterator()) {
        G4MoleculeDefinition *moleculeDef = iterator.value();
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        if (moleculeDef->GetName().substr(0, 3) == "DSB") {
            //-----------------------Sub-Diffusion-----------------------
            if (DrDefinitions::Instance()->GetIsSubDiffusion()) {
                DrWaitingTime *SubDiff = new DrWaitingTime();
                moleculeDef->GetProcessManager()->AddDiscreteProcess(SubDiff, 1);
                DrTransportation_Subdiffusion_CTRW *subdiff = new DrTransportation_Subdiffusion_CTRW();
                ph->RegisterProcess(subdiff, moleculeDef);
            }
                //-------------------Brownian-Diffusion-----------------------
            else {
                G4DNABrownianTransportation *brown = new G4DNABrownianTransportation();
                ph->RegisterProcess(brown, moleculeDef);
            }
        }
            //-----------Brownian-Diffusion-For-Other-Particles---------------
        else {
            //---------------------Below from M.K.------------------------
            if (moleculeDef != G4H2O::Definition()) {
                G4DNABrownianTransportation *brown = new G4DNABrownianTransportation();
                ph->RegisterProcess(brown, moleculeDef);
            } else {
                moleculeDef->GetProcessManager()->AddRestProcess(
                        new G4DNAElectronHoleRecombination(), 2);
                G4DNAMolecularDissociation *dissociationProcess =
                        new G4DNAMolecularDissociation("H2O_DNAMolecularDecay");
                dissociationProcess->SetDisplacer(
                        moleculeDef, new G4DNAWaterDissociationDisplacer);
                dissociationProcess->SetVerboseLevel(1);
                moleculeDef->GetProcessManager()->AddRestProcess(dissociationProcess,
                                                                 1);
            }
            //--------------------Above from M.K.------------------------
        }
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        /*
         * Warning : end of particles and processes are needed by
         * EM Physics builders
         */
    }

    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    //@@@@ DaMaRiS Processes
    DrClock::Definition()->GetProcessManager()->AddDiscreteProcess(new DrReportMSD(), 10);
    DrClock::Definition()->GetProcessManager()->AddDiscreteProcess(new DrReportSystem(), 8);
    DrClock::Definition()->GetProcessManager()->AddDiscreteProcess(new DrFRAP(), 6);
    for(auto read: DrDefinitions::Instance()->GetProcList()){
        read.first->GetProcessManager()->AddDiscreteProcess(read.second);
    }
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    G4DNAChemistryManager::Instance()->Initialize();
}

void DrEmDNAChemistry::ConstructTimeStepModel(
        G4DNAMolecularReactionTable *reactionTable) {
    //=========================================
    // Diffusion controlled reaction model
    //=========================================
    //
    // The reaction model defines how to compute the reaction range between
    // molecules
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    G4VDNAReactionModel *reactionRadiusComputer = new DrDNASmoluchowskiReactionModel();
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    reactionTable->PrintTable(reactionRadiusComputer);

    // The StepByStep model tells the step manager how to behave before and
    // after each step, how to compute the time steps.
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    DrDNAMolecularStepByStepModel *stepByStep = new DrDNAMolecularStepByStepModel();
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    stepByStep->SetReactionModel(reactionRadiusComputer);

    RegisterTimeStepModel(stepByStep, 0);
}

G4bool DrEmDNAChemistry::IsWaterMolecularConfigurationActivated() {
    return fSetWaterConfiguration;
}


G4String DrEmDNAChemistry::GetFullParmName(G4String suffix ) {
    return "Ch/" + fName + "/" + suffix;
}


void DrEmDNAChemistry::Quit(G4String parName, G4String message) {
    std::cerr << "TOPAS is exiting due to an error in Chemistry configuration." << std::endl;
    std::cerr << "Parameter: " << parName << " " << message << std::endl;
    exit(1);
}