// Physics Module for TsEmDNAChemistry
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

#include "TsEmDNAChemistry.hh"
#include "TsParameterManager.hh"
#include "G4DNAOneStepThermalizationModel.hh"
#include "TsDNARuddIonisationExtendedModel.hh"
#include "G4DNAMolecularStepByStepModel.hh"
#include "G4DNASmoluchowskiReactionModel.hh"

#include "TsDNAFirstOrderReaction.hh"
#include "TsDNARemoveInMaterial.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"

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
#include "G4VDNAReactionModel.hh"
#include "G4DNAElectronHoleRecombination.hh"
#include "G4DNAIonisation.hh"

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
#include "TsScavengerProduct.hh"

#include "G4PhysicsListHelper.hh"
#include "G4BuilderType.hh"
#include "G4DNAMoleculeEncounterStepper.hh"
#include "G4ProcessVector.hh"
#include "G4ProcessTable.hh"


#include "G4MolecularConfiguration.hh"
#include "G4PhysicsConstructorFactory.hh"
#include "G4Threading.hh"

#include <map>

G4_DECLARE_PHYSCONSTR_FACTORY(TsEmDNAChemistry);

TsEmDNAChemistry::TsEmDNAChemistry()
:G4VUserChemistryList()
{
    G4DNAChemistryManager::Instance()->SetChemistryList(this);
}

TsEmDNAChemistry::TsEmDNAChemistry(TsParameterManager* pM)
:G4VUserChemistryList(), fPm(pM)
{
    DefineParameters();
    fPm->SetNeedsChemistry();
    G4DNAChemistryManager::Instance()->SetChemistryList(this);
}


TsEmDNAChemistry::~TsEmDNAChemistry()
{
}

void TsEmDNAChemistry::DefineParameters()
{
    fExistingMolecules["solvatedelectron"] = "e_aq";
    fExistingMolecules["hydroxyl"] = "OH";
    fExistingMolecules["hydrogen"] = "H";
    fExistingMolecules["hydronium"] = "H3Op";
    fExistingMolecules["dyhydrogen"] = "H2";
    fExistingMolecules["hydroxide"] = "OHm";
    fExistingMolecules["hydrogenperoxide"] = "H2O2";
    
    fDiffusionCoefficients["e_aq"] = 4.9e-9*m2/s;
    fDiffusionCoefficients["OH"] = 2.8e-9*m2/s;
    fDiffusionCoefficients["H"] = 7.0e-9*m2/s;
    fDiffusionCoefficients["H3Op"] = 9.0e-9*m2/s;
    fDiffusionCoefficients["H2"] = 4.8e-9*m2/s;
    fDiffusionCoefficients["OHm"] = 5.0e-9*m2/s;
    fDiffusionCoefficients["H2O2"] = 2.3e-9*m2/s;
    fDiffusionCoefficients["Product"] = 0.0*m2/s;
    
    fName = "Default";
    if ( fPm->ParameterExists("Ch/ChemistryName") )
        fName = fPm->GetStringParameter("Ch/ChemistryName");
    
    // Re-set diffusion coefficients
    if ( fPm->ParameterExists(GetFullParmName("DiffusionCoefficients/Molecules")) ) {
        G4String* molecules = fPm->GetStringVector(GetFullParmName("DiffusionCoefficients/Molecules"));
        G4double* diffusionCoefficients = fPm->GetDoubleVector(GetFullParmName("DiffusionCoefficients/Values"),"surface perTime");
        G4int nbOfMolecules = fPm->GetVectorLength(GetFullParmName("DiffusionCoefficients/Molecules"));
        
        if ( nbOfMolecules != fPm->GetVectorLength(GetFullParmName("DiffusionCoefficients/Values")) )
            Quit(GetFullParmName("DiffusionCoefficients/Molecules"),
                 "must match with " + GetFullParmName("DiffusionCoefficients/Values"));
        
        for ( int i = 0; i < nbOfMolecules; i++ ) {
            molecules[i].toLower();
            if ( !MoleculeExists(molecules[i]) ) {
                Quit(GetFullParmName("DiffusionCoefficients/Molecules"), "\n--- Molecule: " + molecules[i] + " was not found in database");
            } else if ( molecules[i] == "water" ) {
                Quit(GetFullParmName("DiffusionCoefficients/Molecules"), "\n--- Molecule: " + molecules[i] + " is not allowed to has a diffusion coefficient.");
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
    G4int prefixLength = G4String("Ch/" + fName + "/Reaction/").length();
    
    if ( numberOfReactions > 0 ) {
        for ( int i = 0; i < numberOfReactions; i++ ) {
            G4String parName = (*reactionNames)[i];
			if (fPm->ParameterExists(parName.substr(0,parName.find("ReactionRate")-1) + "/CompatibleWithStepByStep") &&
				!fPm->GetBooleanParameter(parName.substr(0,parName.find("ReactionRate")-1) + "/CompatibleWithStepByStep"))
				continue;
			
            G4String reactions = parName.substr(prefixLength, parName.find("ReactionRate")-prefixLength-1);
            G4String reactorA = reactions.substr(0, reactions.find("/"));
            G4String reactorB = reactions.substr(reactions.find("/") + 1);
            G4String* product = fPm->GetStringVector(parName.substr(0,parName.find("ReactionRate")-1) + "/Products");
            G4int nbOfProduct = fPm->GetVectorLength(parName.substr(0,parName.find("ReactionRate")-1) + "/Products");
            G4double reactionRate = fPm->GetDoubleParameter(parName.substr(0,parName.find("ReactionRate")-1) + "/ReactionRate","perMolarConcentration perTime");
            
            std::vector<G4String> reactors;
            std::vector<G4String> products;
            
            reactorA.toLower();
            reactorB.toLower();
            if ( !MoleculeExists(reactorA) )
                Quit(parName, "\n--- Molecule: " + reactorA + " was not found in database");
            
            if ( !MoleculeExists(reactorB) )
                Quit(parName, "\n--- Molecule: " + reactorB + " was not found in database");
            
            reactors.push_back(fExistingMolecules[reactorA]);
            reactors.push_back(fExistingMolecules[reactorB]);
            
            for ( int j = 0; j < nbOfProduct; j++ ) {
                product[j].toLower();
                if ( product[j] == "none" ) // Note on comparison of strings: See comment on ConstructReactionTable()
                    products.push_back("none");
                else if ( product[j] == "water") {
                    products.push_back("H2O");
                } else {
                    if ( !MoleculeExists(product[j]) )
                        Quit(parName, "\n--- Molecule: " + product[j] + " was not found in database");
                    else {
                        products.push_back(fExistingMolecules[product[j]]);
                    }
                }
            }
            
            fReactionSpecies.push_back(reactors);
            fReactionProducts.push_back(products);
            fReactionRates.push_back(reactionRate);
        }
    }

    else {
        Quit("Ch/"+fName +"/Reaction/","No Reactions were declared for chemistry!");
    } 

    // First order reactions (Scavengers)
    if ( fPm->ParameterExists(GetFullParmName("Scavenger/Molecules")) ) {
        G4String* scavengedMolecules = fPm->GetStringVector(GetFullParmName("Scavenger/Molecules"));
        G4double* scavengerConcentrations = fPm->GetDoubleVector(GetFullParmName("Scavenger/Concentrations"), "molar concentration");
        G4double* scavengerReactionRates = fPm->GetDoubleVector(GetFullParmName("Scavenger/ReactionRates"), "perMolarConcentration perTime");
        G4bool* hasProduct = fPm->GetBooleanVector(GetFullParmName("Scavenger/HasProducts"));
        G4int nbOfScavenged = fPm->GetVectorLength(GetFullParmName("Scavenger/Molecules"));
        for ( int i = 0; i < nbOfScavenged; i++ ) {
            G4String aName = scavengedMolecules[i];
            aName.toLower();
            if ( !MoleculeExists(aName) )  {
                Quit(GetFullParmName("Scavenger/Molecules"), "\n--- Molecule name " + aName + " was not found");
            }
            fScavengedMolecules.push_back( fExistingMolecules[aName] );
            fScavengedCapacities.push_back( scavengerConcentrations[i] * scavengerReactionRates[i] );
            fScavengerHasProducts.push_back( hasProduct[i] );
        }
    }
	// Remove molecules in material named (at contact, no following a power low as First order reactions)
	if ( fPm->ParameterExists(GetFullParmName("RemoveInMaterial"))) {
		fRemoveInMaterial = fPm->GetStringParameter(GetFullParmName("RemoveInMaterial"));
		G4String* theMolecules = fPm->GetStringVector(GetFullParmName("RemoveInMaterialTheseMolecules"));
		G4int nTheMaterials = fPm->GetVectorLength(GetFullParmName("RemoveInMaterialTheseMolecules"));
		for ( int i = 0; i < nTheMaterials; i++ ) {
			if ( !MoleculeExists(theMolecules[i]) )
				Quit(GetFullParmName("RemoveInMaterialTheseMolecules"), "\n--- Molecule: " + theMolecules[i] + " was not found in database");
			fRemoveInMaterialTheseMolecules.push_back(fExistingMolecules[theMolecules[i]]);
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
    if ( fPm->ParameterExists(parName + "RydbergStatesAndDiffuseBands/AutoIoinization"))
        fRydDiffAutoIonization = fPm->GetUnitlessParameter(parName + "RydbergStatesAndDiffuseBands/AutoIoinization");
    std::cout << "-- Branching ratio: Rydberg states and diffuse bands auto-ionization probability " << fRydDiffAutoIonization << std::endl;
    
    fRydDiffRelaxation = 0.5;
    if ( fPm->ParameterExists(parName + "RydbergStatesAndDiffuseBands/Relaxation"))
        fRydDiffRelaxation = fPm->GetUnitlessParameter(parName + "RydbergStatesAndDiffuseBands/Relaxation");
    std::cout << "-- Branching ratio: Rydberg states and diffuse bands relaxation probability " << fRydDiffRelaxation << std::endl;
    
    fDissociativeAttachment = 1.0;
    if ( fPm->ParameterExists(parName + "DissociativeAttachment"))
        fDissociativeAttachment = fPm->GetUnitlessParameter(parName + "DissociativeAttachment");
    std::cout << "-- Branching ratio: dissociative attachment probability " << fDissociativeAttachment << std::endl;
    
}


G4bool TsEmDNAChemistry::MoleculeExists(G4String name) {
    if ( fExistingMolecules.find(name) == fExistingMolecules.end() )
        return false;
    return true;
}


void TsEmDNAChemistry::ConstructMolecule()
{
    G4Electron::Definition();
    G4H2O::Definition();
    G4Hydrogen::Definition();
    G4H3O::Definition();
    G4OH::Definition();
    G4Electron_aq::Definition();
    G4H2O2::Definition();
    G4H2::Definition();
    TsScavengerProduct::Definition();
    
    G4MoleculeTable::Instance()->CreateConfiguration("H", G4Hydrogen::Definition());
    G4MoleculeTable::Instance()->CreateConfiguration("OH", G4OH::Definition());
    G4MoleculeTable::Instance()->CreateConfiguration("H2O2", G4H2O2::Definition());
    G4MoleculeTable::Instance()->CreateConfiguration("H2", G4H2::Definition());
    G4MoleculeTable::Instance()->CreateConfiguration("e_aq", G4Electron_aq::Definition());
    G4MoleculeTable::Instance()->CreateConfiguration("H3Op", G4H3O::Definition());
    
    G4MolecularConfiguration* OHm =
    G4MoleculeTable::Instance()-> CreateConfiguration("OHm", G4OH::Definition(), -1, 5.0e-9 * (m2 / s));
    OHm->SetMass(17.0079 * g / Avogadro * c_squared);
    
    G4MoleculeTable::Instance()->CreateConfiguration("Product", TsScavengerProduct::Definition());
    
    if ( fSetWaterConfiguration )
        G4MoleculeTable::Instance()->CreateConfiguration("H2O", G4H2O::Definition());
    
    // Re-set diffusion coefficients
    std::map<G4String, G4double>::iterator it = fDiffusionCoefficients.begin();
    while ( it != fDiffusionCoefficients.end() ) {
        G4String molecule = it->first;
        G4double diffusionCoefficient = it->second;
        G4MoleculeTable::Instance()->GetConfiguration(molecule)->SetDiffusionCoefficient(diffusionCoefficient);
        std::cout << "-- Diffusion coefficient for molecule: " << molecule
        << " used the value " << diffusionCoefficient/(m2/s) << " (m2/s) " << std::endl;
        it++;
    }
    std::cout << "" << std::endl;
}


void TsEmDNAChemistry::ConstructDissociationChannels()
{
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

void TsEmDNAChemistry::ConstructReactionTable(G4DNAMolecularReactionTable*
                                              theReactionTable)
{
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
        reactionData = new G4DNAMolecularReactionData(fReactionRates[t], reactions[fReactionSpecies[t][0]], reactions[fReactionSpecies[t][1]]);
        G4int reactionType = 1;
        reactionData->SetReactionID(reactionType);
        
        for ( size_t u = 0; u < fReactionProducts[t].size(); u++ ) {
            if ( "none" != fReactionProducts[t][u] ) // This comparison crashes if the order is fReactionProducts[t][u] != "none"
                reactionData->AddProduct(reactions[ fReactionProducts[t][u] ] );
        }
        std::cout << " Re-set reaction kobs to : " << fReactionRates[t]/(1e-3*m3/(mole*s)) << "/M/s" << std::endl;
        theReactionTable->SetReaction(reactionData);
    }
}


#include "G4NistManager.hh"

void TsEmDNAChemistry::ConstructProcess()
{
    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
    G4VProcess* process = G4ProcessTable::GetProcessTable()->FindProcess("e-_G4DNAVibExcitation", "e-");
    
    if (process != nullptr) {
        G4DNAVibExcitation* vibExcitation = (G4DNAVibExcitation*) process;
        G4VEmModel* model = vibExcitation->EmModel();
        G4DNASancheExcitationModel* sancheExcitationMod =
        dynamic_cast<G4DNASancheExcitationModel*>(model);
        if(sancheExcitationMod != nullptr) {
            sancheExcitationMod->ExtendLowEnergyLimit(0.025 * eV);
        }
    }
    
    process = G4ProcessTable::GetProcessTable()->FindProcess("e-_G4DNAElectronSolvation", "e-");
    if ( process == nullptr ) {
        G4DNAElectronSolvation* solvation = (G4DNAElectronSolvation*)process;
        G4VEmModel* solvationModel = new G4DNAOneStepThermalizationModel();
        solvation->SetEmModel(solvationModel, 0);
    }
    
    G4bool useExtendedProtonLimit = false;
    if ( fPm->ParameterExists(GetFullParmName("UseExtendedProtonLimit") ) )
        useExtendedProtonLimit = fPm->GetBooleanParameter(GetFullParmName("UseExtendedProtonLimit"));
    
    if (useExtendedProtonLimit ) {
        process = G4ProcessTable::GetProcessTable()->FindProcess("proton_G4DNAIonisation", "proton");
        if ( process ) {
            std::cout << "********************************* Use with caution ************************************** " << std::endl;
            std::cout << "************** Extending the energy range for proton ionisation interactions ************ " << std::endl;
            G4DNAIonisation* protonIonisationProcess = (G4DNAIonisation*)process;
            
            G4VEmModel* mod1 = new TsDNARuddIonisationExtendedModel();
            mod1->SetLowEnergyLimit(0*eV);
            mod1->SetHighEnergyLimit(500*keV);
            
            G4VEmModel* mod2;
            mod2= new G4DNABornIonisationModel();
            mod2->SetLowEnergyLimit(500*keV);
            mod2->SetHighEnergyLimit(100*MeV);
            
            G4VEmModel* mod3 = new TsDNARuddIonisationExtendedModel();
            mod3->SetLowEnergyLimit(100*MeV);
            mod3->SetHighEnergyLimit(1000*MeV);
            
            //protonIonisationProcess->SetEmModel(mod1,0);
            protonIonisationProcess->SetEmModel(mod2,0);
            protonIonisationProcess->SetEmModel(mod3,1);
        }
    }
    
    G4MoleculeTable* theMoleculeTable = G4MoleculeTable::Instance();
    G4MoleculeDefinitionIterator iterator =
    theMoleculeTable->GetDefintionIterator();
    iterator.reset();
    
    std::vector<G4MolecularConfiguration*> pr;
    pr.push_back(G4MoleculeTable::Instance()->GetConfiguration("Product"));
    while (iterator()) {
        G4MoleculeDefinition* moleculeDef = iterator.value();
        
        if (moleculeDef != G4H2O::Definition() && moleculeDef != TsScavengerProduct::Definition() ) {
            G4DNABrownianTransportation* brown = new G4DNABrownianTransportation();
            ph->RegisterProcess(brown, moleculeDef);
            
            for ( size_t u = 0; u < fScavengedMolecules.size(); u++ ) {
                if ( moleculeDef->GetName() == fScavengedMolecules[u] ) {
                    TsDNAFirstOrderReaction* firstOrd = new TsDNAFirstOrderReaction();
                    
                    G4MolecularConfiguration* mC = G4MoleculeTable::Instance()->GetConfiguration(fScavengedMolecules[u]);
                    if ( !fScavengerHasProducts[u] )
                        firstOrd->SetReaction(mC, fScavengedCapacities[u]);
                    else
                        firstOrd->SetReaction(mC, pr, fScavengedCapacities[u]);
                    moleculeDef->GetProcessManager()->AddDiscreteProcess(firstOrd);
                    std::cout << "-- Set scavenging capacity for molecule " << fScavengedMolecules[u] <<
					" of " << fScavengedCapacities[u]*s << " /s " << std::endl;
                }
            }
	    for ( size_t u = 0; u < fRemoveInMaterialTheseMolecules.size(); u++ ) {
				if ( moleculeDef->GetName() == G4MoleculeTable::Instance()->
					GetConfiguration(fRemoveInMaterialTheseMolecules[u])->GetDefinition()->GetName())
				{
					TsDNARemoveInMaterial* removeProcess = new TsDNARemoveInMaterial();
					G4MolecularConfiguration* mC = G4MoleculeTable::Instance()->
												   GetConfiguration(fRemoveInMaterialTheseMolecules[u]);
					if ( !G4Material::GetMaterial(fRemoveInMaterial) )
						Quit("", "\n--- Material " + fRemoveInMaterial + " was not found " );
					removeProcess->SetReaction(mC, G4Material::GetMaterial(fRemoveInMaterial));
					moleculeDef->GetProcessManager()->AddDiscreteProcess(removeProcess);
					G4cout << "-- Molecule " << mC->GetName() << " will be removed at contact of material " << fRemoveInMaterial << G4endl;
				}
	    }
            
        } else if ( moleculeDef == G4H2O::Definition() ) {
            moleculeDef->GetProcessManager()
            ->AddRestProcess(new G4DNAElectronHoleRecombination(), 2);
            G4DNAMolecularDissociation* dissociationProcess =
            new G4DNAMolecularDissociation("H2O_DNAMolecularDecay");
            dissociationProcess->SetDisplacer(moleculeDef, new G4DNAWaterDissociationDisplacer);
            dissociationProcess->SetVerboseLevel(1);
            moleculeDef->GetProcessManager()->AddRestProcess(dissociationProcess, 1);
        }
    }
    G4DNAChemistryManager::Instance()->Initialize();
}


void TsEmDNAChemistry::ConstructTimeStepModel(G4DNAMolecularReactionTable* reactionTable)
{
    G4VDNAReactionModel* reactionRadiusComputer =
    new G4DNASmoluchowskiReactionModel();
    reactionTable->PrintTable(reactionRadiusComputer);
    
    G4DNAMolecularStepByStepModel* stepByStep =
    new G4DNAMolecularStepByStepModel();
    stepByStep->SetReactionModel(reactionRadiusComputer);
    
    RegisterTimeStepModel(stepByStep, 0);
}


G4bool TsEmDNAChemistry::IsWaterMolecularConfigurationActivated() {
    return fSetWaterConfiguration;
}


G4String TsEmDNAChemistry::GetFullParmName(G4String suffix ) {
    return "Ch/" + fName + "/" + suffix;
}


void TsEmDNAChemistry::Quit(G4String parName, G4String message) {
    G4cerr << "TOPAS is exiting due to an error in Chemistry configuration." << std::endl;
    G4cerr << "--- Parameter: " << parName << " " << message << std::endl;
    fPm->AbortSession(1);
}

