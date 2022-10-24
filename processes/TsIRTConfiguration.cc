// Extra Class for TsIRT
#include "TsIRTConfiguration.hh"
#include "TsIRTUtils.hh"
#include "TsParameterManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4RandomDirection.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "G4UIcommand.hh"
#include "G4IosFlagsSaver.hh"

#include <vector>
#include <map>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <functional>
#include <sstream>

TsIRTConfiguration::TsIRTConfiguration(G4String name, TsParameterManager* pM)
: fPm(pM), fName(name), fReactionID(0), fTotalBinaryReaction(0), fLastMoleculeID(0),
fKick(false), fAllTotallyDiffusionControlled(false)
{
	G4String chemistryList = fPm->GetStringParameter("Ch/ChemistryName");
	
	fUseSimpleScavengerModel = true;
	
	fUtils = new TsIRTUtils();
	fLowerTime = 1.0e-13*s;
	fUpperTime = 1.0e-6*s;
	
	fExistingMolecules.clear();
	
	fExistingMolecules["hydrogen"] = "H^0";
	fExistingMolecules["hydroxyl"] = "OH^0";
	fExistingMolecules["hydrogenperoxide"] = "H2O2^0";
	fExistingMolecules["dyhydrogen"] = "H_2^0";
	fExistingMolecules["solvatedelectron"] = "e_aq^-1";
	fExistingMolecules["hydronium"] = "H3O^1";
	fExistingMolecules["hydroxide"] = "OH^-1";
	fExistingMolecules["oxygen"] = "O2^0";
	fExistingMolecules["superoxideanion"] = "O2^-1";
	fExistingMolecules["hydroperoxy"] = "HO2^0";
	fExistingMolecules["dioxidanide"] = "HO2^-1";
	fExistingMolecules["atomicoxygen"] = "O3P^0";
	fExistingMolecules["oxyde"] = "O^-1";
	fExistingMolecules["trioxide"] = "O3^-1";
	fExistingMolecules["ozone"] = "O3^0";
	fExistingMolecules["none"] = "None";
	
	AddMolecule("H^0",     7.0e9*nm*nm/s,   0, 0.19*nm);
	AddMolecule("OH^0",    2.2e9*nm*nm/s,   0, 0.22*nm);
	AddMolecule("H2O2^0",  2.3e9*nm*nm/s,   0, 0.21*nm);
	AddMolecule("H_2^0",   4.8e9*nm*nm/s,   0, 0.14*nm);
	AddMolecule("e_aq^-1", 4.9e9*nm*nm/s,  -1, 0.50*nm);
	AddMolecule("H3O^1",   9.46e9*nm*nm/s,  1, 0.25*nm);
	AddMolecule("OH^-1",   5.3e9*nm*nm/s,  -1, 0.33*nm);
	AddMolecule("O2^0",    2.4e9*nm*nm/s,   0, 0.17*nm);
	AddMolecule("O2^-1",   1.75e9*nm*nm/s, -1, 0.22*nm);
	AddMolecule("HO2^0",   2.3e9*nm*nm/s,  0, 0.21*nm);
	AddMolecule("HO2^-1",  1.4e9*nm*nm/s, -1, 0.25*nm);
	AddMolecule("O3P^0",   2.0e9*nm*nm/s,  0, 0.20*nm);
	AddMolecule("O^-1",    2.0e9*nm*nm/s, -1, 0.25*nm);
	AddMolecule("O3^-1",   2.0e9*nm*nm/s, -1, 0.20*nm);
	AddMolecule("O3^0",    2.0e9*nm*nm/s,  0, 0.20*nm);
	AddMolecule("None",    0.0e9*nm*nm/s,  0, 0.00*nm);
	
	// Re-set diffusion coefficient, radius and/or charge.
	std::vector<G4String>* moleculeNames = new std::vector<G4String>;
	fPm->GetParameterNamesStartingWith("Mo/", moleculeNames);
	G4int numberOfMolecules = moleculeNames->size();
	std::vector<G4String> moleculesDontExist;
	
	for ( int i = 0; i < numberOfMolecules; i++ ) {
		G4String fullName = (*moleculeNames)[i];
		fullName.toLower();
		G4bool moleculeExists = false;
		
		if ( fullName.contains("diffusioncoefficient") ) {
			G4String molName = fullName.substr(3, fullName.find("diffusioncoefficient")-4);
			if ( fExistingMolecules.find(molName) != fExistingMolecules.end() ) {
				G4int molID =  fMoleculesID[fExistingMolecules[molName]];
				fMoleculesDefinition[molID].diffusionCoefficient = fPm->GetDoubleParameter(fullName,"surface perTime");
				G4cout << molName << ": Re-set diffusion coefficient to "
				<< fMoleculesDefinition[molID].diffusionCoefficient/(nm*nm/s) << " nm2/s" << std::endl;
				moleculeExists = true;
			}
		}
		
		if (fullName.contains("radius") ) {
			G4String molName = fullName.substr(3, fullName.find("radius")-4);
			if ( fExistingMolecules.find(molName) != fExistingMolecules.end() ) {
				G4int molID =  fMoleculesID[fExistingMolecules[molName]];
				fMoleculesDefinition[molID].radius = fPm->GetDoubleParameter(fullName, "Length");
				G4cout << molName << ": Re-set radius to "
				<< fMoleculesDefinition[molID].radius/nm << " nm" << std::endl;
				moleculeExists = true;
			}
		}
		
		if (fullName.contains("charge")) {
			G4String molName = fullName.substr(3, fullName.find("charge")-4);
			if ( fExistingMolecules.find(molName) != fExistingMolecules.end() ) {
				G4int molID =  fMoleculesID[fExistingMolecules[molName]];
				fMoleculesDefinition[molID].charge = fPm->GetUnitlessParameter(fullName);
				G4cout << molName << ": Re-set charge to "
				<< fMoleculesDefinition[molID].radius/nm << " nm" << std::endl;
				moleculeExists = true;
			}
		}
		
		if (fullName.contains("symbol") && !moleculeExists) {
			moleculesDontExist.push_back(fullName.substr(3, fullName.find("symbol")-4));
		}
	}
	
	// Creates user-defined molecules.
	for ( size_t u = 0; u < moleculesDontExist.size(); u++ ) {
		moleculesDontExist[u].toLower();
		G4double charge = fPm->GetUnitlessParameter("Mo/" + moleculesDontExist[u] + "/Charge");
		G4double radius = fPm->GetDoubleParameter("Mo/" + moleculesDontExist[u] + "/Radius", "Length");
		G4double diffusionCoefficient = fPm->GetDoubleParameter("Mo/" + moleculesDontExist[u] +
																"/DiffusionCoefficient","surface perTime");
		G4String symbol = fPm->GetStringParameter("Mo/" + moleculesDontExist[u] + "/Symbol");
		
		if ( fPm->ParameterExists("Mo/" + moleculesDontExist[u] + "/AssignMoleculeID"))
			AddMolecule(symbol, fPm->GetIntegerParameter("Mo/" + moleculesDontExist[u] + "/AssignMoleculeID"),
						diffusionCoefficient, charge, radius);
		else
			AddMolecule(symbol, diffusionCoefficient, charge, radius);
		
		fExistingMolecules[moleculesDontExist[u]] = symbol;
	}
	
	G4String parName = "Ch/" + chemistryList + "/SetAllReactionsTotallyDiffusionControlled";
	if ( fPm->ParameterExists(parName) )
		fAllTotallyDiffusionControlled = fPm->GetBooleanParameter(parName);
	
	std::vector<G4String>* reactionNames = new std::vector<G4String>;
	fPm->GetParameterNamesBracketedBy("Ch/" , "Products", reactionNames);
	G4int numberOfReactions = reactionNames->size();
	G4int prefixLength = G4String("Ch/" + chemistryList + "/Reaction/").length();
	
	for ( int i = 0; i < numberOfReactions; i++ ) {
		G4String aparName = (*reactionNames)[i];
		if ( aparName.contains("BackgroundReaction"))
			continue;
		
		if ( fPm->ParameterExists(aparName.substr(0,aparName.find("Products")-1) + "/Active") &&
			!fPm->GetBooleanParameter(aparName.substr(0,aparName.find("Products")-1) + "/Active") )
			continue;
		
		G4String reactions = aparName.substr(prefixLength, aparName.find("Products")-prefixLength-1);
		G4String reactorA = reactions.substr(0, reactions.find("/"));
		G4String reactorB = reactions.substr(reactions.find("/") + 1);
		reactorA.toLower();
		reactorB.toLower();
		G4String* product = fPm->GetStringVector(aparName);
		G4int nbOfProduct = fPm->GetVectorLength(aparName);
		
		for ( int j = 0; j < nbOfProduct; j++ )
			product[j].toLower();
		
		reactorA.toLower();
		reactorB.toLower();
		
		G4int reactionType = fPm->GetIntegerParameter(aparName.substr(0,aparName.find("Products")-1) + "/ReactionType");
		
		
		G4double reactionRate = fPm->GetDoubleParameter(aparName.substr(0,aparName.find("Products")-1) +
														"/ReactionRate","perMolarConcentration perTime");
		
		if ( nbOfProduct == 1 )
			InsertReaction(reactorA, reactorB,
						   product[0],"None","None",reactionRate,reactionType);
		else if ( nbOfProduct == 2 )
			InsertReaction(reactorA, reactorB,product[0],product[1],"None",reactionRate,reactionType);
		else
			InsertReaction(reactorA, reactorB,product[0],product[1],product[2],reactionRate,reactionType);

		if (fPm->ParameterExists(aparName.substr(0,aparName.find("Products")-1) + "/ActivationRate") ||
			fPm->ParameterExists(aparName.substr(0,aparName.find("Products")-1) + "/DiffusionRate")){
			if (fReactions[fReactionID-1].reactionType == 2 || fReactions[fReactionID-1].reactionType == 4) {
				G4double activationRate = fPm->GetDoubleParameter(aparName.substr(0,aparName.find("Products")-1) +
															"/ActivationRate","perMolarConcentration perTime");
				G4double diffusionRate  = fPm->GetDoubleParameter(aparName.substr(0,aparName.find("Products")-1) +
															"/DiffusionRate","perMolarConcentration perTime");
				fReactions[fReactionID-1].kdif = diffusionRate;
				fReactions[fReactionID-1].kact = activationRate;
				G4cout << "Advance Reaction Mode for reaction" << fReactionID-1 << " : " 
				       << fExistingMolecules[reactorA] << " + " << fExistingMolecules[reactorB] << G4endl;
			}
			else {
				G4cout << "WARNING: Advance Reaction Options Only available for reations type II or IV" << G4endl;
			}
		}
	}
	
	prefixLength = G4String("Ch/" + chemistryList + "/BackgroundReaction/").length();
	for ( int i = 0; i < numberOfReactions; i++ ) {
		G4String aparName = (*reactionNames)[i];
		if ( aparName.contains("/Reaction/"))
			continue;
		
		if ( fPm->ParameterExists(aparName.substr(0,aparName.find("Products")-1) + "/Active") &&
			!fPm->GetBooleanParameter(aparName.substr(0,aparName.find("Products")-1) + "/Active") )
			continue;
		
		G4String reactions = aparName.substr(prefixLength, aparName.find("Products")-prefixLength-1);
		G4String reactorA = reactions.substr(0, reactions.find("/"));
		G4String reactorB = reactions.substr(reactions.find("/") + 1);
		reactorA.toLower();
		reactorB.toLower();
		G4String* product = fPm->GetStringVector(aparName);
		G4int nbOfProduct = fPm->GetVectorLength(aparName);
		
		for ( int j = 0; j < nbOfProduct; j++ )
			product[j].toLower();
		
		reactorA.toLower();
		reactorB.toLower();
		
		if ( fPm->ParameterExists(aparName.substr(0,aparName.find("Products")-1) + "/ScavengingCapacity") ) {
			G4double scavengingCapacity = fPm->GetDoubleParameter(aparName.substr(0,aparName.find("Products")-1) +
																  "/ScavengingCapacity","perTime");
			if ( nbOfProduct == 1 )
				InsertBackgroundReaction(reactorA, reactorB,product[0],"None","None",scavengingCapacity,false);
			else if ( nbOfProduct == 2 )
				InsertBackgroundReaction(reactorA, reactorB,product[0],product[1],"None",scavengingCapacity,false);
			else
				InsertBackgroundReaction(reactorA, reactorB,product[0],product[1],product[2],scavengingCapacity,false);
		} else {
			G4double concentration = fPm->GetDoubleParameter(aparName.substr(0,aparName.find("Products")-1) +
															 "/Concentration","molar concentration");
			G4double reactionRate = fPm->GetDoubleParameter(aparName.substr(0,aparName.find("Products")-1) +
															"/ReactionRate","perMolarConcentration perTime");
			
			G4String scavengingExponentialModel = "exponentialsinglefactor";
			if ( fPm->ParameterExists(aparName.substr(0,aparName.find("Products")-1) + "/ScavengingModel"))
				scavengingExponentialModel = fPm->GetStringParameter(aparName.substr(0,aparName.find("Products")-1) + "/ScavengingModel");
			
			scavengingExponentialModel.toLower();
			if ( scavengingExponentialModel == "exponentialsinglefactor" ) {
				if ( nbOfProduct == 1 )
					InsertBackgroundReaction(reactorA, reactorB,product[0],"None","None",reactionRate,concentration,false);//scavengingCapacity);
				else if ( nbOfProduct == 2 )
					InsertBackgroundReaction(reactorA, reactorB,product[0],product[1],"None",reactionRate,concentration,false);//,scavengingCapacity);
				else
					InsertBackgroundReaction(reactorA, reactorB,product[0],product[1],product[2],reactionRate,concentration,false);//scavengingCapacity);
				
			} else if (scavengingExponentialModel == "exponentialdoublefactor" ) {
				if ( nbOfProduct == 1 )
					InsertBackgroundReaction(reactorA, reactorB,product[0],"None","None",reactionRate,concentration,true);
				else if ( nbOfProduct == 2 )
					InsertBackgroundReaction(reactorA, reactorB,product[0],product[1],"None",reactionRate,concentration,true);
				else
					InsertBackgroundReaction(reactorA, reactorB,product[0],product[1],product[2],
											 reactionRate,concentration,true);
			} else {
				Quit(aparName.substr(0,aparName.find("Products")-1) + "/ScavengingModel",
					 "Scavenging model does not exists in TOPAS-nBio database\n Use ExponentialSingleFactor or ExponentialDoubleFactor");
			}
		}
	}
	
	// Reactions for quality check of the IRT loop.
	parName = "Ch/" + chemistryList + "/TestIRTForQualityAssurance";
	fQualityAssurance = false;
	if (fPm->ParameterExists(parName) && fPm->GetBooleanParameter(parName) ) {
		fQualityAssurance = true;
		G4double kUnit = fPm->GetUnitValue("/M/s");
		// Fake molecules for quality assurance.
		AddMolecule("A",       17, 5.0e9*nm*nm/s,  0, 1.00*nm);
		AddMolecule("B",       18, 5.0e9*nm*nm/s,  0, 1.00*nm);
		AddMolecule("AA",      19, 5.0e9*nm*nm/s,  0, 1.00*nm);
		AddMolecule("AB",      20, 5.0e9*nm*nm/s,  0, 1.00*nm);
		AddMolecule("BB",      21, 5.0e9*nm*nm/s,  0, 1.00*nm);
		AddMolecule("AAB",     22, 5.0e9*nm*nm/s,  0, 1.00*nm);
		AddMolecule("ABB",     23, 5.0e9*nm*nm/s,  0, 1.00*nm);
		// Reactions for quality assurance.
		InsertReaction("A", "A", "AA", "None", "None", 5.0e9 * kUnit, 1);
		InsertReaction("A", "B", "AB", "None", "None", 5.0e9 * kUnit, 1);
		InsertReaction("B", "B", "BB", "None", "None", 5.0e9 * kUnit, 1);
		InsertReaction("AB", "A", "AAB", "None", "None", 5.0e9 * kUnit, 1);
		InsertReaction("AB", "B", "ABB", "None", "None", 5.0e9 * kUnit, 1);
	}
	
	// Re-scale chemistry parameters based on temperature
	parName = "Ch/" + chemistryList + "/Temperature";
	fScaleForTemperature = false;
	if ( fPm->ParameterExists(parName) ) {
		fTemperature = fPm->GetUnitlessParameter(parName);
		fScaleForTemperature = true;
		
		parName = "Ch/" + chemistryList + "/ApplyCorrectionScalingForTemperature";
		fKick = false;
		if ( fPm->ParameterExists(parName) )
			fKick = fPm->GetBooleanParameter(parName);
		
		//AdjustDiffusionCoefficientyForTemperature(fTemperature);
		AdjustDiffusionCoefficientyForTemperatureArrehniusFit(fTemperature);
	}
	
	ResolveReactionRateCoefficients();
	
	if ( fScaleForTemperature )
		AdjustReactionRateForTemperature(fTemperature);
	
	CalculateContactProbabilities();
	ResolveRemainerReactionParameters();
	
	if ( fPm->ParameterExists("Ch/"+chemistryList+"/ModelAcidPropertiesFromSubstance") ) {
		fpHSolventConcentration = 0.0;
		fpHValue                = 7.1;
		
		fpHSolvent = fPm->GetStringParameter("Ch/"+chemistryList+"/ModelAcidPropertiesFromSubstance");
		fpHSolvent.toLower();
		if (fPm->ParameterExists("Ch/"+chemistryList+"/ModelAcidPropertiesWithConcentration") &&
			fPm->ParameterExists("Ch/"+chemistryList+"/ModelAcidPropertiesWithpH")) {
			G4String message = "Cannot be defined when parameter: Ch/" + chemistryList + "/ModelAcidPropertiesWithConcentration is used.";
			Quit("Ch/"+chemistryList+"/ModelAcidPropertiesWithpH",message);
		}
		
		if (fPm->ParameterExists("Ch/"+chemistryList+"/ModelAcidPropertiesWithConcentration") ) {
			fpHSolventConcentration = fPm->GetDoubleParameter("Ch/"+chemistryList+"/ModelAcidPropertiesWithConcentration","molar concentration");
			AdjustReactionRateForPH("Concentration");
		}
		
		if (fPm->ParameterExists("Ch/"+chemistryList+"/ModelAcidPropertiesWithpH") ) {
			fpHValue = fPm->GetUnitlessParameter("Ch/"+chemistryList+"/ModelAcidPropertiesWithpH");
			AdjustReactionRateForPH("PH");
		}
	}
	
	PrintReactionsInformation();
}


TsIRTConfiguration::~TsIRTConfiguration()
{
	delete fUtils;
}


void TsIRTConfiguration::AddMolecule(G4String name, G4int moleculeID, G4double diffusionCoefficient,
									 G4double charge, G4double radius) {
	TsMoleculeDefinition aMolecule;
	aMolecule.diffusionCoefficient = diffusionCoefficient;
	aMolecule.charge = charge;
	aMolecule.radius = radius;
	
	fMoleculesDefinition[moleculeID] = aMolecule;
	fMoleculesID[name] = moleculeID;
	fMoleculesName[moleculeID] = name;
	
	fLastMoleculeID = moleculeID;
	
	G4bool found = false;
	for (auto it = fExistingMolecules.begin(); it != fExistingMolecules.end(); ++it) {
		if (it->second == name) {
			found = true;
			break;
		}
	}
	
	if ( !found )
		fExistingMolecules[name] = name;
	
}


void TsIRTConfiguration::AddMolecule(G4String name, G4double diffusionCoefficient,
									 G4double charge, G4double radius) {
	fLastMoleculeID++;
	TsMoleculeDefinition aMolecule;
	aMolecule.diffusionCoefficient = diffusionCoefficient;
	aMolecule.charge = charge;
	aMolecule.radius = radius;
	
	fMoleculesDefinition[fLastMoleculeID] = aMolecule;
	fMoleculesID[name] = fLastMoleculeID;
	fMoleculesName[fLastMoleculeID] = name;
	
	G4bool found = false;
	for (auto it = fExistingMolecules.begin(); it != fExistingMolecules.end(); ++it) {
		if (it->second == name) {
			found = true;
			break;
		}
	}
	
	if ( !found )
		fExistingMolecules[name] = name;
	
}


void TsIRTConfiguration::AddMolecule(G4String name) {
	fLastMoleculeID++;
	AddMolecule(name, fLastMoleculeID, 0.0, 0.0, 0.0);
	return;
}


void TsIRTConfiguration::AdjustDiffusionCoefficientyForTemperatureArrehniusFit(G4double temperatureInCelsius) {
	G4double D = 0.0;
	G4double A = 0., B = 0., C =0.;
	for ( auto& aMol : fMoleculesDefinition ) {
		G4int molID = aMol.first;
		G4cout << "Updating diffCoeff for molecule " << std::setw(7) << fMoleculesName[molID] << ": " <<
		std::setprecision(5) << aMol.second.diffusionCoefficient/(nm*nm/s) << " -> ";
		
		if ( fMoleculesName[molID] == "OH^0" ) {
			A = -8120.0; B = 1.45065; C = 11.42870;
		} else if (fMoleculesName[molID] == "e_aq^-1" ) {
			A = -1387.9; B = 1.05709; C = 13.04977;
		} else if (fMoleculesName[molID] == "H3O^1" ) {
			A = -9560.7; B = 1.59117; C = 11.09443;
		} else if (fMoleculesName[molID] == "H2O2^0" ) {
			A = -8122.5; B = 1.45081; C = 11.44698;
		} else if (fMoleculesName[molID] == "H_2^0" ) {
			A = -8221.0; B = 1.45315; C = 11.76385;
		} else if (fMoleculesName[molID] == "H^0" ) {
			A = -8333.0; B = 1.45576; C = 11.92479;
		} else if (fMoleculesName[molID] == "OH^-1" ) {
			A = -113704.1; B = 2.00721; C = 10.95892;
		} else if (fMoleculesName[molID] == "HO2^0" ) {
			A = -16782.7; B = 1.59057; C = 11.30206;
		} else if (fMoleculesName[molID] == "HO2^-1" ) {
			A = -16648.6; B = 1.58897; C = 11.08846;
		} else if (fMoleculesName[molID] == "O2^0" ) {
			A = -16797.5; B = 1.59077; C = 11.32031;
		} else if (fMoleculesName[molID] == "O2^-1" ) {
			A = -16717.0; B = 1.58975; C = 11.18462;
		} else if (fMoleculesName[molID] == "O3P^0" ) {
			A = -16740.6; B = 1.59003; C = 11.24232;
		} else if (fMoleculesName[molID] == "O^-1" ) {
			A = -16740.4; B = 1.59003; C = 11.24233;
		} else if (fMoleculesName[molID] == "O3^-1" ) {
			A = -16740.3; B = 1.59003; C = 11.24233;
		} else if (fMoleculesName[molID] == "O3^0" ) {
			A = -16740.0; B = 1.59002; C = 11.24234;
		} else {
			continue;
		}
		
		if ( D == 0 ) {
			G4double invTemp = 1./(temperatureInCelsius+273.15);
			D = A * std::pow(invTemp, B) + C;
			D = std::pow(10., D) * nm * nm / s;
		}
		
		aMol.second.diffusionCoefficient = D;
		G4cout << std::setprecision(5) << D/(nm*nm/s) << " nm2/s " << G4endl;
		D = 0;
	}
}


void TsIRTConfiguration::AdjustDiffusionCoefficientyForTemperature(G4double temperatureInCelsius) {
	G4double DH2O_T = 1.0704e-9 + 4.7082e-11*temperatureInCelsius - 1.2633e-14 * std::pow(temperatureInCelsius,2) +
	5.0605e-15 * std::pow(temperatureInCelsius,3) - 2.3048e-17 * std::pow(temperatureInCelsius,4) +
	2.9217e-20 * std::pow(temperatureInCelsius,5) + 1.3128e-23 * std::pow(temperatureInCelsius,6);
	G4double ambientTemperature = 25.0;
	G4double DH2O_25C = 1.0704e-9 + 4.7082e-11*ambientTemperature - 1.2633e-14 * std::pow(ambientTemperature,2) +
	5.0605e-15 * std::pow(ambientTemperature,3) - 2.3048e-17 * std::pow(ambientTemperature,4) +
	2.9217e-20 * std::pow(ambientTemperature,5) + 1.3128e-23 * std::pow(ambientTemperature,6);
	
	G4double D = 0.0;
	for ( auto& aMol : fMoleculesDefinition ) {
		G4int molID = aMol.first;
		G4cout << "Updating diffCoeff for molecule " << std::setw(7) << fMoleculesName[molID] << ": " <<
		std::setprecision(5) << aMol.second.diffusionCoefficient/(nm*nm/s) << " -> ";
		
		if ( fMoleculesName[molID] == "OH^-1" ) {
			D = 2.666e-9 + 9.769e-11 * temperatureInCelsius +
			3.303e-13 * temperatureInCelsius*temperatureInCelsius -
			7.295e-16 * temperatureInCelsius*temperatureInCelsius*temperatureInCelsius;
			aMol.second.diffusionCoefficient = D * m * m / s;
			
		} else if ( fMoleculesName[molID] == "H3O^1" ) {
			D = 5.361e-9 + 1.659e-10 * temperatureInCelsius -
			7.48e-14 * temperatureInCelsius*temperatureInCelsius -
			1.0018e-16 * temperatureInCelsius*temperatureInCelsius*temperatureInCelsius;
			aMol.second.diffusionCoefficient = D * m * m / s;
			
		} else if ( fMoleculesName[molID] == "e_aq^-1" ) {
			if (temperatureInCelsius <= 90) {
				D = -4.690 + 1.8e-2 * temperatureInCelsius -
				1.253e-4 * temperatureInCelsius*temperatureInCelsius +
				5.458e-7 * temperatureInCelsius*temperatureInCelsius*temperatureInCelsius;
				D = std::pow(10, D) * 1E-4;
				aMol.second.diffusionCoefficient = D * m * m / s;
			}
			
			else {
				D = Arrhenius(1.98373E-05, temperatureInCelsius + 273.15, 20.75);
				aMol.second.diffusionCoefficient = D * m * m / s;
			}
			
		} else {
			D = aMol.second.diffusionCoefficient * DH2O_T/DH2O_25C;
			aMol.second.diffusionCoefficient = D;
		}
		G4cout << std::setprecision(5) << aMol.second.diffusionCoefficient/(nm*nm/s) << " nm2/s " << G4endl;
	}
}


G4bool TsIRTConfiguration::MoleculeExists(G4String name) {
	if ( fMoleculesID.find(name) == fMoleculesID.end() )
		return false;
	return true;
}


void TsIRTConfiguration::QuitIfMoleculeNotFound(G4String mol) {
	if (fExistingMolecules.find(mol) == fExistingMolecules.end()) {
		G4cerr << "TOPAS is exiting due to a fatal error in IRT Chemistry setup!" << G4endl;
		G4cerr << "--- Molecule " << mol << " does not exists in database" << G4endl;
		fPm->AbortSession(1);
	}
	return;
}


G4double TsIRTConfiguration::GetMoleculeRadius(G4int moleculeID) {
	return fMoleculesDefinition[moleculeID].radius;
}


G4int TsIRTConfiguration::GetMoleculeCharge(G4int moleculeID) {
	return fMoleculesDefinition[moleculeID].charge;
}


void TsIRTConfiguration::Diffuse(TsMolecule& mol, G4double dt) {
	G4double sigma, x, y, z;
	G4double diffusionCoefficient = fMoleculesDefinition[mol.id].diffusionCoefficient;
	
	sigma = std::sqrt(2.0 * diffusionCoefficient * dt);
	
	x = G4RandGauss::shoot(0., 1.0)*sigma;
	y = G4RandGauss::shoot(0., 1.0)*sigma;
	z = G4RandGauss::shoot(0., 1.0)*sigma;
	mol.position += G4ThreeVector(x, y, z);
	mol.time += dt;
}


void TsIRTConfiguration::PrintMoleculesInformation() {
	for ( auto& molecules : fMoleculesID ) {
		G4String name = molecules.first;
		G4int id = molecules.second;
		TsMoleculeDefinition aMolecule = fMoleculesDefinition[id];
		std::cout << "----" << std::endl;
		std::cout << "Molecule name           : " << name << std::endl;
		std::cout << "  Molecule ID           : " << id << std::endl;
		std::cout << "  Diffusion coefficient : " << aMolecule.diffusionCoefficient/(nm*nm/s)
		<< " nm2/s " << std::endl;
		std::cout << "  Charge                : " << aMolecule.charge << " e " << std::endl;
		std::cout << "  Radius                : " << aMolecule.radius/nm << " nm " << std::endl;
		std::cout << "----" << std::endl;
		std::cout << "" << std::endl;
	}
}


void TsIRTConfiguration::PrintReactionsInformation() {
	G4IosFlagsSaver iosfs(G4cout);
	std::map<size_t, std::vector<TsMolecularReaction> > temporal;
	for ( size_t i = 0; i < fReactions.size(); i++) {
		G4int type = fReactions[i].reactionType;
		temporal[type].push_back(fReactions[i]);
	}
	
	G4String* outputReactionParentA = new G4String[fReactions.size()];
	G4String* outputReactionParentB = new G4String[fReactions.size()];
	G4String* outputReactionProducts = new G4String[fReactions.size()];
	G4String* outputKobs = new G4String[fReactions.size()];
	G4String* outputKdif = new G4String[fReactions.size()];
	G4String* outputKact = new G4String[fReactions.size()];
	G4String* outputScav = new G4String[fReactions.size()];
	G4String* outputAlpha = new G4String[fReactions.size()];
	G4String* outputReactionRadius = new G4String[fReactions.size()];
	G4String* outputReactionRadiusEff = new G4String[fReactions.size()];
	G4String* outputReactionRadiusEff1 = new G4String[fReactions.size()];
	G4String* outputProbability = new G4String[fReactions.size()];
	G4String* outputModel = new G4String[fReactions.size()];
	
	G4cout << G4endl;
	G4int n = 0;
	for ( size_t i = 1; i <= temporal.size(); i++ ) {
		for ( auto& reactions : temporal[i] ) {
			G4int type = reactions.reactionType;
			G4int molA = reactions.reactorA;
			G4int molB = reactions.reactorB;
			std::vector<G4int> products = reactions.products;
			
			outputReactionParentA[n] = fMoleculesName[molA];
			outputReactionParentB[n] = fMoleculesName[molB];
			
			if ( products.size() > 0 ) {
				outputReactionProducts[n] = fMoleculesName[products[0]];
				for ( size_t u = 1; u < products.size(); u++ ) {
					outputReactionProducts[n] += " + " + fMoleculesName[products[u]];
				}
			} else {
				outputReactionProducts[n] = "None";
			}
			
			outputKobs[n] = G4UIcommand::ConvertToString(reactions.kobs/fPm->GetUnitValue("/M/s"));
			outputKdif[n] = G4UIcommand::ConvertToString(reactions.kdif/fPm->GetUnitValue("/M/s"));
			outputKact[n] = G4UIcommand::ConvertToString(reactions.kact/fPm->GetUnitValue("/M/s"));
			outputAlpha[n] = G4UIcommand::ConvertToString(reactions.alpha * nm);
			outputReactionRadius[n] = G4UIcommand::ConvertToString(reactions.reactionRadius/nm);
			outputReactionRadiusEff[n] = G4UIcommand::ConvertToString(reactions.effectiveReactionRadius/nm);
			outputReactionRadiusEff1[n] = G4UIcommand::ConvertToString(reactions.effectiveTildeReactionRadius/nm);
			outputProbability[n] = G4UIcommand::ConvertToString(reactions.probabilityOfReaction);
			if (type == 6 ) {
				outputScav[n] = G4UIcommand::ConvertToString(reactions.scavengingCapacity * s);
				outputModel[n] = reactions.sampleExponential ? "true" : "false";
			} else {
				outputScav[n] = "0.0";
				outputModel[n] = "none";
			}
			n++;
		}
	}
	
	G4int maxLengthOutputReactionParentA = -1;
	G4int maxLengthOutputReactionParentB = -1;
	G4int maxLengthOutputReactionProducts = -1;
	G4int maxLengthOutputKobs = -1;
	G4int maxLengthOutputKdif = -1;
	G4int maxLengthOutputKact = -1;
	G4int maxLengthOutputProbability = -1;
	G4int maxLengthOutputReactionRadius = -1;
	G4int maxLengthOutputReactionRadiusEff = -1;
	G4int maxLengthOutputReactionRadiusEff1 = -1;
	G4int maxLengthOutputAlpha = -1;
	
	for ( int i = 0; i < n; i++ ) {
		if ( maxLengthOutputReactionParentA < (G4int)outputReactionParentA[i].length())
			maxLengthOutputReactionParentA = outputReactionParentA[i].length();
		
		if ( maxLengthOutputReactionParentB < (G4int)outputReactionParentB[i].length())
			maxLengthOutputReactionParentB = outputReactionParentB[i].length();
		
		if ( maxLengthOutputReactionProducts < (G4int)outputReactionProducts[i].length())
			maxLengthOutputReactionProducts = outputReactionProducts[i].length();
		
		if ( maxLengthOutputKobs < (G4int)outputKobs[i].length())
			maxLengthOutputKobs = outputKobs[i].length();
		
		if ( maxLengthOutputKdif < (G4int)outputKdif[i].length())
			maxLengthOutputKdif = outputKdif[i].length();
		
		if ( maxLengthOutputKact < (G4int)outputKact[i].length())
			maxLengthOutputKact = outputKact[i].length();
		
		if ( maxLengthOutputAlpha < (G4int)outputAlpha[i].length())
			maxLengthOutputAlpha = outputAlpha[i].length();
		
		if ( maxLengthOutputProbability < (G4int) outputProbability[i].length())
			maxLengthOutputProbability = outputProbability[i].length();
		
		if ( maxLengthOutputReactionRadius < (G4int) outputReactionRadius[i].length())
			maxLengthOutputReactionRadius = outputReactionRadius[i].length();
		
		if ( maxLengthOutputReactionRadiusEff < (G4int) outputReactionRadiusEff[i].length())
			maxLengthOutputReactionRadiusEff = outputReactionRadiusEff[i].length();
		
		if ( maxLengthOutputReactionRadiusEff1 < (G4int) outputReactionRadiusEff1[i].length())
			maxLengthOutputReactionRadiusEff1 = outputReactionRadiusEff1[i].length();
	}
	
	/*maxLengthOutputReactionParentA += 2;
	 maxLengthOutputReactionParentB += 2;*/
	maxLengthOutputReactionProducts += 2;
	maxLengthOutputKobs += 4;
	maxLengthOutputKdif += 4;
	maxLengthOutputKact += 4;
	maxLengthOutputAlpha += 2;
	 maxLengthOutputProbability += 2;
	 maxLengthOutputReactionRadius += 2;
	 maxLengthOutputReactionRadiusEff += 2;
	 maxLengthOutputReactionRadiusEff1 += 2;
	
	G4String* title = new G4String[11];
	title[0] = "Type";
	title[1] = "Molecules";
	title[2] = "Products";
	title[3] = "Kobs (/M/s)";
	title[4] = "Kdif (/M/s)";
	title[5] = "Kact (/M/s)";
	title[6] = "r (nm)";
	title[7] = "Preact";
	title[8] = "alpha (/nm)";
	title[9] = "reff (nm)";
	title[10] = "reff1 (nm)";
	n = 0;
	for ( size_t i = 1; i <= temporal.size(); i++ ) {
		for (int k = 0; k < 120; k++ )
			G4cout << "=";
		G4cout << G4endl;
		if ( i == 1 )
			G4cout << std::setfill(' ') << std::setw(7) << title[0]
			<< std::setw(maxLengthOutputReactionParentA+maxLengthOutputReactionParentB+3) << std::left << title[1]
			<< std::setw(maxLengthOutputReactionProducts+7) << std::left << title[2]
			<< std::setw(maxLengthOutputKobs) << std::left << title[3]
			<< std::setw(maxLengthOutputReactionRadius) << std::left << title[6]
			<< std::setw(maxLengthOutputProbability) << std::left << title[7] << G4endl;
		else if ( i == 2 )
			G4cout << std::setfill(' ') << std::setw(7) << title[0]
			<< std::setw(maxLengthOutputReactionParentA+maxLengthOutputReactionParentB+3) << std::left << title[1]
			<< std::setw(maxLengthOutputReactionProducts+7) << std::left << title[2]
			<< std::setw(maxLengthOutputKobs) << std::left << title[3]
			<< std::setw(maxLengthOutputKdif) << std::left << title[4]
			<< std::setw(maxLengthOutputKact) << std::left << title[5]
			<< std::setw(maxLengthOutputReactionRadius) << std::left << title[6]
			<< std::setw(maxLengthOutputProbability) << std::left << title[7]
			<< std::setw(maxLengthOutputAlpha) << std::left << title[8] << G4endl;
		else if ( i == 3 )
			G4cout << std::setfill(' ') << std::setw(7) << title[0]
			<< std::setw(maxLengthOutputReactionParentA+maxLengthOutputReactionParentB+3) << std::left << title[1]
			<< std::setw(maxLengthOutputReactionProducts+7) << std::left << title[2]
			<< std::setw(maxLengthOutputKobs) << std::left << title[3]
			<< std::setw(maxLengthOutputReactionRadius) << std::left << title[6]
			<< std::setw(maxLengthOutputReactionRadiusEff) << std::left << title[9]
			<< std::setw(maxLengthOutputProbability) << std::left << title[7] << G4endl;
		else if ( i == 4 )
			G4cout << std::setfill(' ') << std::setw(7) << title[0]
			<< std::setw(maxLengthOutputReactionParentA+maxLengthOutputReactionParentB+3) << std::left << title[1]
			<< std::setw(maxLengthOutputReactionProducts+7) << std::left << title[2]
			<< std::setw(maxLengthOutputKobs) << std::left << title[3]
			<< std::setw(maxLengthOutputKdif) << std::left << title[4]
			<< std::setw(maxLengthOutputKact) << std::left << title[5]
			<< std::setw(maxLengthOutputReactionRadius) << std::left << title[6]
			<< std::setw(maxLengthOutputReactionRadiusEff) << std::left << title[9]
			<< std::setw(maxLengthOutputReactionRadiusEff1) << std::left << title[10]
			<< std::setw(maxLengthOutputProbability) << std::left << title[7] << G4endl;
		else if ( i == 5 )
			G4cout << std::setfill(' ') << std::setw(7) << title[0]
			<< std::setw(maxLengthOutputReactionParentA+maxLengthOutputReactionParentB+3) << std::left << title[1]
			<< std::setw(maxLengthOutputReactionProducts+7) << std::left << title[2]
			<< std::setw(maxLengthOutputKobs) << std::left << title[3]
			<< std::setw(maxLengthOutputReactionRadius) << std::left << title[6]
			<< std::setw(maxLengthOutputReactionRadiusEff) << std::left << title[9]
			<< std::setw(maxLengthOutputProbability) << std::left << title[7] << G4endl;
		else if ( i == 6 )
			G4cout << std::setfill(' ') << std::setw(7) << title[0]
			<< std::setw(maxLengthOutputReactionParentA+maxLengthOutputReactionParentB+3) << std::left << title[1]
			<< std::setw(maxLengthOutputReactionProducts+7) << std::left << title[2]
			<< std::setw(maxLengthOutputKobs) << std::left << "Scavenging capacity /s" << G4endl;
		
		for (int k = 0; k < 120; k++ )
			G4cout << "=";
		G4cout << G4endl;
	}
	
	delete[] title;
	delete[] outputProbability;
	delete[] outputReactionRadius;
	delete[] outputReactionRadiusEff;
	delete[] outputReactionRadiusEff1;
	delete[] outputKobs;
	delete[] outputKdif;
	delete[] outputKact;
	delete[] outputScav;
	delete[] outputAlpha;
	delete[] outputReactionProducts;
	delete[] outputReactionParentA;
	delete[] outputReactionParentB;
	delete[] outputModel;
	temporal.clear();
}


G4int TsIRTConfiguration::GetReactionIndex(G4int pdgA, G4int pdgB) {
	for ( size_t u = 0; u < fMoleculeCanReactWith[pdgA].size(); u++ ) {
		if ( pdgB == fMoleculeCanReactWith[pdgA][u].first )
			return fMoleculeCanReactWith[pdgA][u].second;
	}
	return -1;
}


G4double TsIRTConfiguration::GetOnsagerRadius(G4int molA, G4int molB) {
	if ( !fScaleForTemperature )
		return (fMoleculesDefinition[molA].charge * fMoleculesDefinition[molB].charge)/
		(4*CLHEP::pi*epsilon0*k_Boltzmann) / (293.15 * 80.1);
	else {
		G4double temperatureInKelvin = fTemperature + 273.15;
		G4double epsilon =((5321.0 * std::pow(temperatureInKelvin,-1)) + (233.76 * std::pow(temperatureInKelvin,+0))
						   - (0.9297 * std::pow(temperatureInKelvin,+1)) + (0.001417*std::pow(temperatureInKelvin,+2))
						   - (8.292E-7*std::pow(temperatureInKelvin,+3))) *  8.85418781762037E-12;
		
		G4double electronCharge = 1.60217662e-19 ;
		G4double rc = std::pow(electronCharge,2) / (4 * CLHEP::pi * epsilon * 1.3806488E-23 * temperatureInKelvin)
		* fMoleculesDefinition[molA].charge * fMoleculesDefinition[molB].charge;
		return rc*1E9*nm;
	}
}


void TsIRTConfiguration::ResolveRemainerReactionParameters() {
	for ( int i = 0; i < (int)fReactions.size(); i++ ) {
		G4int molA = fReactions[i].reactorA;
		G4int molB = fReactions[i].reactorB;
		G4int reactionType = fReactions[i].reactionType;
		
		G4double observedReactionRate = fReactions[i].kobs;
		G4double diffusionReactionRate = fReactions[i].kdif;
		G4double activationReactionRate = fReactions[i].kact;
		G4double probability = fReactions[i].probabilityOfReaction;
		G4double rc = fReactions[i].OnsagerRadius;
		
		G4double reactionRadius = 0;
		G4double effectiveReactionRadius = 0;
		G4double effectiveTildeReactionRadius = 0;
		G4double alpha = 0;
		
		if (reactionType == 1 || reactionType == 3 || reactionType == 5) {
			G4double sumDiffCoeff = 0;
			if (molA == molB) {
				sumDiffCoeff = fMoleculesDefinition[molA].diffusionCoefficient;
				effectiveReactionRadius = observedReactionRate / (4. * CLHEP::pi * sumDiffCoeff * CLHEP::Avogadro);
			} else {
				sumDiffCoeff = fMoleculesDefinition[molA].diffusionCoefficient
				+ fMoleculesDefinition[molB].diffusionCoefficient;
				effectiveReactionRadius = observedReactionRate / (4. * CLHEP::pi * sumDiffCoeff * CLHEP::Avogadro);
			}
			
			effectiveReactionRadius /= probability;
			
			if ( rc == 0 ) {
				reactionRadius = effectiveReactionRadius;
			} else {
				reactionRadius = rc/std::log(1 + rc/effectiveReactionRadius);
			}
			
		} else if ( reactionType == 2 || reactionType == 4 ) {
			// R = RA + RB, see Plante 2011 after Eq.17
			reactionRadius = fMoleculesDefinition[molA].radius + fMoleculesDefinition[molB].radius;
			G4double sumDiffCoeff = fMoleculesDefinition[molA].diffusionCoefficient + fMoleculesDefinition[molB].diffusionCoefficient;
			
			if ( reactionType == 2 ) {
				effectiveReactionRadius = reactionRadius;
				effectiveTildeReactionRadius = effectiveReactionRadius;
				alpha = (activationReactionRate + diffusionReactionRate)/(diffusionReactionRate*reactionRadius);
				
			} else {
				effectiveReactionRadius = -rc / (1-exp(rc / reactionRadius));
				rc /= nm;
				G4double kact = activationReactionRate/fPm->GetUnitValue("/M/s");
				G4double r = reactionRadius/nm;
				G4double reff = effectiveReactionRadius/nm;
				G4double v = 0.93*kact/(4 * CLHEP::pi * std::pow(reff,2)*6.022140857e-1);
				sumDiffCoeff /= nm*nm/s;
				effectiveTildeReactionRadius =  rc/(std::exp(rc/r) * (1 + sumDiffCoeff*rc/(std::pow(r,2)*v))-1) * nm;
				
				G4double nm3persToPerMPers = 6.022140857e-1;
				kact /= nm3persToPerMPers;
				kact *= nm*nm*nm/s;
				rc *= nm;
				sumDiffCoeff *= nm*nm/s;
				alpha = 0.25*kact/(CLHEP::pi*reactionRadius*reactionRadius) +
				rc*sumDiffCoeff/(reactionRadius*reactionRadius*(1.0 - std::exp(-rc/reactionRadius)));
			}
		} else {
			continue; // TypeVI
		}
		
		if ( fQualityAssurance ) {
			reactionRadius = 1.0 * nm;
			effectiveReactionRadius = 1.0 * nm;
			probability = 1;
			fReactions[i].reactionRadius = reactionRadius;
			fReactions[i].effectiveReactionRadius = effectiveReactionRadius;
			fReactions[i].probabilityOfReaction = probability;
		}
		
		fReactions[i].reactionRadius = reactionRadius;
		fReactions[i].effectiveReactionRadius = effectiveReactionRadius;
		fReactions[i].effectiveTildeReactionRadius = effectiveTildeReactionRadius;
		fReactions[i].alpha = alpha;
	}
}


void TsIRTConfiguration::CalculateContactProbabilities() {
	for ( int i = 0; i < (int)fReactions.size(); i++ ) {
		G4double probability = 1;
		G4int reactionType = fReactions[i].reactionType;
		G4int molA = fReactions[i].reactorA;
		G4int molB = fReactions[i].reactorB;
		if ( reactionType == 1 || reactionType == 3 || reactionType == 5 ) {
			if ( reactionType == 5 )
				probability = 0.25;
			else
				probability = 1.0;
			
		} else if ( reactionType == 2 || reactionType == 4 ){
			G4double Rs = 0.29 * nm;
			G4double diffusionReactionRate = fReactions[i].kdif;
			G4double activationReactionRate = fReactions[i].kact;
			G4double rc = GetOnsagerRadius(molA, molB);
			
			G4double effectiveReactionRadius = -rc / (1-exp(rc /
															(fMoleculesDefinition[molA].radius + fMoleculesDefinition[molB].radius)));
			
			if ( reactionType == 2 )
				probability =  Rs / (Rs + (diffusionReactionRate / activationReactionRate) *
									 ((fMoleculesDefinition[molA].radius + fMoleculesDefinition[molB].radius) + Rs));
			else
				probability = Rs / (Rs + (diffusionReactionRate / activationReactionRate) * (effectiveReactionRadius + Rs));
			
		} else {
			continue;
		}
		
		fReactions[i].probabilityOfReaction = probability;
	}
}


void TsIRTConfiguration::ResolveReactionRateCoefficients() {
	for ( int i = 0; i < (int)fReactions.size(); i++ ) {
		G4int molA = fReactions[i].reactorA;
		G4int molB = fReactions[i].reactorB;
		G4int reactionType = fReactions[i].reactionType;
		G4double kobs = fReactions[i].kobs;
		
		G4double rc = GetOnsagerRadius(molA, molB);
		G4double diffusionReactionRate = 0;
		G4double activationReactionRate = 0;
		fReactions[i].OnsagerRadius = rc;

		if (fReactions[i].kdif != -1 && fReactions[i].kact != -1) {
			continue;
		}
		
		if (reactionType == 1 || reactionType == 3 || reactionType == 5) {
			diffusionReactionRate = kobs;
			
		} else if ( reactionType == 2 || reactionType == 4 ) {
			
			// R = RA + RB, see Plante 2011 after Eq.17
			fReactionRadius = fMoleculesDefinition[molA].radius + fMoleculesDefinition[molB].radius;
			G4double sumDiffCoeff = fMoleculesDefinition[molA].diffusionCoefficient + fMoleculesDefinition[molB].diffusionCoefficient;
			
			if ( reactionType == 2 ) {
				diffusionReactionRate = 4 * pi * sumDiffCoeff * (fMoleculesDefinition[molA].radius +
																 fMoleculesDefinition[molB].radius) * Avogadro;
				if (molA == molB)
					diffusionReactionRate/=2;
				
				activationReactionRate = diffusionReactionRate * kobs / (diffusionReactionRate - kobs);
				
			} else {
				G4double effectiveReactionRadius = -rc / (1-exp(rc /
																(fMoleculesDefinition[molA].radius +
																 fMoleculesDefinition[molB].radius)));
				diffusionReactionRate = 4 * pi * sumDiffCoeff * effectiveReactionRadius * Avogadro;
				
				if (molA == molB) diffusionReactionRate/=2;
				
				activationReactionRate = diffusionReactionRate * kobs / (diffusionReactionRate - kobs);
			}
		} else {
			continue; //Type VI;
		}
		
		fReactions[i].kdif = diffusionReactionRate;
		fReactions[i].kact = activationReactionRate;
	}
}


void TsIRTConfiguration::InsertReaction(G4int molA, G4int molB, std::vector<G4int> products,
										G4double kobs, G4int reactionType)
{
	G4int index = fReactionID;
	fReactionID++;
	fTotalBinaryReaction++;
	
	TsMolecularReaction aMolecularReaction;
	aMolecularReaction.reactorA = molA;
	aMolecularReaction.reactorB = molB;
	aMolecularReaction.products = products;
	aMolecularReaction.index = index;
	aMolecularReaction.kact  = -1;
	aMolecularReaction.kdif  = -1;
	
	if (fAllTotallyDiffusionControlled) {
		if ( reactionType == 2 ) {
			reactionType = 1;
		} else if ( reactionType == 4 ) {
			reactionType = 3;
		} else if ( reactionType == 5 ) {
			if ( fMoleculesDefinition[molA].charge * fMoleculesDefinition[molB].charge != 0 )
				reactionType = 3;
			else
				reactionType = 1;
		}
	}
	
	aMolecularReaction.reactionType = reactionType;
	aMolecularReaction.kobs = kobs;
	aMolecularReaction.scavengingCapacity = 0.0;
	
	G4int pdgA = molA;
	G4int pdgB = molB;
	if ( fMoleculeCanReactWith.find(pdgA) == fMoleculeCanReactWith.end() ) { // key not found
		fMoleculeCanReactWith[pdgA].push_back(std::make_pair(pdgB, index));
	} else {
		G4bool found = false;
		for ( size_t u = 0; u < fMoleculeCanReactWith[pdgA].size(); u++ ) {
			if ( pdgB == fMoleculeCanReactWith[pdgA][u].first) {
				found = true;
				break;
			}
		}
		if ( !found ) {
			fMoleculeCanReactWith[pdgA].push_back(std::make_pair(pdgB,index));
		}
	}
	
	if ( fMoleculeCanReactWith.find(pdgB) == fMoleculeCanReactWith.end() ) { // key not found
		fMoleculeCanReactWith[pdgB].push_back(std::make_pair(pdgA,index));
	} else {
		G4bool found = false;
		for ( size_t u = 0; u < fMoleculeCanReactWith[pdgB].size(); u++ ) {
			if ( pdgA == fMoleculeCanReactWith[pdgB][u].first ) {
				found = true;
				break;
			}
		}
		if ( !found ) {
			fMoleculeCanReactWith[pdgB].push_back(std::make_pair(pdgA,index));
		}
	}
	
	fReactions[index] = aMolecularReaction;
}


void TsIRTConfiguration::InsertReaction(G4String A, G4String B, G4String p1, G4String p2, G4String p3,
										G4double kobs, G4int reactionType)
{
	QuitIfMoleculeNotFound(A);
	QuitIfMoleculeNotFound(B);
	G4String molNameA  = fExistingMolecules[A];
	G4String molNameB  = fExistingMolecules[B];

	std::vector<G4int> products;
	G4String molNameP1 = p1;
	if (p1 != "None") {
		QuitIfMoleculeNotFound(p1);
		molNameP1 = fExistingMolecules[p1];
		products.push_back(fMoleculesID[molNameP1]);
	}

	G4String molNameP2 = p2;
	if (p2 != "None") {
		QuitIfMoleculeNotFound(p2);
		molNameP2 = fExistingMolecules[p2];
		products.push_back(fMoleculesID[molNameP2]);
	}

	G4String molNameP3 = p3;
	if (p3 != "None") {
		QuitIfMoleculeNotFound(p3);
		molNameP3 = fExistingMolecules[p3];
		products.push_back(fMoleculesID[molNameP3]);
	}

	G4int molA = fMoleculesID[molNameA];
	G4int molB = fMoleculesID[molNameB];
	
	if ( GetReactionIndex(molA, molB) > -1 ) {
		G4cerr << "TOPAS is exiting due to a fatal error in IRT Chemistry setup!" << G4endl;
		G4cerr << "--- Reaction: " << molNameA << " + " << molNameB << " -> " 
		                           << molNameP1 << " + " << molNameP2 << " + " << molNameP3 << " already exists." << G4endl;
		fPm->AbortSession(1);
		return;
	}
	
	InsertReaction(molA, molB, products, kobs, reactionType);
	
}


void TsIRTConfiguration::InsertBackgroundReaction(G4String A, G4String B, G4String p1,
												  G4String p2, G4String p3,
												  G4double scavengingCapacity, G4bool sampleExponential)
{
	QuitIfMoleculeNotFound(A);
	QuitIfMoleculeNotFound(B);
	G4String molNameA = fExistingMolecules[A];
	G4String molNameB = fExistingMolecules[B];
	
	G4int index = fReactionID;
	fReactionID++;
	G4int molA = fMoleculesID[molNameA];
	G4int molB = fMoleculesID[molNameB];
	// No need to test to confirm if this reaction already exists, because it could exists as a second order reaction

	std::vector<G4int> products;
	G4String molNameP1 = p1;
	if (p1 != "None") {
		QuitIfMoleculeNotFound(p1);
		molNameP1 = fExistingMolecules[p1];
		products.push_back(fMoleculesID[molNameP1]);
	}

	G4String molNameP2 = p2;
	if (p2 != "None") {
		QuitIfMoleculeNotFound(p2);
		molNameP2 = fExistingMolecules[p2];
		products.push_back(fMoleculesID[molNameP2]);
	}

	G4String molNameP3 = p3;
	if (p3 != "None") {
		QuitIfMoleculeNotFound(p3);
		molNameP3 = fExistingMolecules[p3];
		products.push_back(fMoleculesID[molNameP3]);
	}
	
	TsMolecularReaction aMolecularReaction;
	aMolecularReaction.reactorA = molA;
	aMolecularReaction.reactorB = molB;
	aMolecularReaction.products = products;
	aMolecularReaction.scavengingCapacity = scavengingCapacity;
	aMolecularReaction.kobs = 0.0;
	aMolecularReaction.index = index;
	aMolecularReaction.reactionType = 6;
	aMolecularReaction.sampleExponential = sampleExponential;
	
	fReactions[index] = aMolecularReaction;
}


void TsIRTConfiguration::InsertBackgroundReaction(G4String A, G4String B, G4String p1,
												  G4String p2, G4String p3, G4double kobs,
												  G4double concentration, G4bool sampleExponential)
{
	QuitIfMoleculeNotFound(A);
	QuitIfMoleculeNotFound(B);
	G4String molNameA = fExistingMolecules[A];
	G4String molNameB = fExistingMolecules[B];
	
	G4int index = fReactionID;
	fReactionID++;
	G4int molA = fMoleculesID[molNameA];
	G4int molB = fMoleculesID[molNameB];
	// No test to confirm if this reaction already exists, it could be a first order reaction
	std::vector<G4int> products;
	G4String molNameP1 = p1;
	if (p1 != "None") {
		QuitIfMoleculeNotFound(p1);
		molNameP1 = fExistingMolecules[p1];
		products.push_back(fMoleculesID[molNameP1]);
	}

	G4String molNameP2 = p2;
	if (p2 != "None") {
		QuitIfMoleculeNotFound(p2);
		molNameP2 = fExistingMolecules[p2];
		products.push_back(fMoleculesID[molNameP2]);
	}

	G4String molNameP3 = p3;
	if (p3 != "None") {
		QuitIfMoleculeNotFound(p3);
		molNameP3 = fExistingMolecules[p3];
		products.push_back(fMoleculesID[molNameP3]);
	}
	
	TsMolecularReaction aMolecularReaction;
	aMolecularReaction.reactorA = molA;
	aMolecularReaction.reactorB = molB;
	aMolecularReaction.products = products;
	aMolecularReaction.kobs = kobs;
	G4double sumDiffCoeff = fMoleculesDefinition[molA].diffusionCoefficient
	+ fMoleculesDefinition[molB].diffusionCoefficient;
	aMolecularReaction.effectiveReactionRadius = kobs / (4. * CLHEP::pi * sumDiffCoeff * CLHEP::Avogadro);
	aMolecularReaction.concentration = concentration;
	aMolecularReaction.scavengingCapacity = kobs * concentration;
	aMolecularReaction.index = index;
	aMolecularReaction.reactionType = 6;
	aMolecularReaction.sampleExponential = sampleExponential;// true;
	fReactions[index] = aMolecularReaction;
}


G4double TsIRTConfiguration::Arrhenius(G4double A, G4double temperatureInKelvin, G4double E) {
	G4double den = 0.0083146*temperatureInKelvin;
	return A * std::exp(-E/den);
}


G4double TsIRTConfiguration::Noyes(G4double kobs, G4double kact, G4double kdiff) {
	if (kobs == 0)
		return (kdiff * kact) / (kact + kdiff);
	
	else if (kact == 0)
		return  (kobs * kdiff) / (kdiff - kobs);
	
	else if (kdiff == 0 )
		return (kobs * kact ) / (kact - kobs);
	
	else
		return 0;
}


G4double TsIRTConfiguration::Smoluchowski(G4double Beta, TsMolecularReaction Reaction) {
	G4double r = fMoleculesDefinition[Reaction.reactorA].radius +
	fMoleculesDefinition[Reaction.reactorB].radius;
	G4double D = fMoleculesDefinition[Reaction.reactorA].diffusionCoefficient +
	fMoleculesDefinition[Reaction.reactorB].diffusionCoefficient;
	
	G4double diffusionReactionRate = 4 * pi * Beta * Avogadro * D * r;
	if (Reaction.reactorA == Reaction.reactorB) diffusionReactionRate/=2;
	
	return diffusionReactionRate/fPm->GetUnitValue("/M/s");
}


G4double TsIRTConfiguration::Debye(G4double Beta, TsMolecularReaction Reaction, G4double temperatureInCelsius) {
	G4int molA = Reaction.reactorA;
	G4int molB = Reaction.reactorB;
	
	G4double sumDiffCoeff = fMoleculesDefinition[molA].diffusionCoefficient +
	fMoleculesDefinition[molB].diffusionCoefficient;
	G4double effectiveReactionRadius = fMoleculesDefinition[molA].radius +
	fMoleculesDefinition[molB].radius;
	
	G4double fD = DebyeFactor(temperatureInCelsius, molA, molB, fMoleculesDefinition[molA].radius +
							  fMoleculesDefinition[molB].radius);
	
	G4double diffusionReactionRate = 4 * pi * sumDiffCoeff * Beta * effectiveReactionRadius * Avogadro * fD;
	
	if (molA == molB) diffusionReactionRate/=2;
	
	return diffusionReactionRate/fPm->GetUnitValue("/M/s");
}


G4double TsIRTConfiguration::DebyeFactor(G4double , G4int molA, G4int molB, G4double r) {
	G4double rc = GetOnsagerRadius(molA, molB);
	G4double d = rc/r;
	return d/(std::exp(d) - 1);
}


G4double TsIRTConfiguration::ProbabilityOfReactionT(TsMolecularReaction reaction, G4double , G4double ) {
	G4double rc = GetOnsagerRadius(reaction.reactorA, reaction.reactorB);
	G4double Rs = 0.29 * nm;
	G4double kdiv = 1;
	if ( reaction.reactionType == 2 || reaction.reactionType == 4) {
		G4double r = fMoleculesDefinition[reaction.reactorA].radius +
		fMoleculesDefinition[reaction.reactorB].radius;
		
		if ( reaction.reactionType == 4 )
			r = -rc / (1-exp(rc / r));
		
		kdiv = reaction.kdif / reaction.kact;
		return Rs/(Rs + kdiv*(r + Rs));
		
	} else {
		return reaction.probabilityOfReaction;
	}
}


G4double TsIRTConfiguration::OnsagerRadius(G4double temperatureInKelvin) {// <- erase
	G4double epsilon =((5321.0 * pow(temperatureInKelvin,-1))
					   + (233.76 * pow(temperatureInKelvin,+0))
					   - (0.9297 * pow(temperatureInKelvin,+1))
					   + (0.001417*pow(temperatureInKelvin,+2))
					   - (8.292E-7*pow(temperatureInKelvin,+3)))
	*  8.85418781762037E-12;
	G4double electronCharge = 1.60217662E-19;
	G4double RC = pow(electronCharge,2) / (4 * pi * epsilon * 1.3806488E-23 * temperatureInKelvin);
	
	return RC * m;
}


void TsIRTConfiguration::AdjustReactionRateForTemperature(G4double temperatureInCelsius) {
	G4double temperatureInKelvin = 273.15 + temperatureInCelsius;
	
	std::vector<std::vector<G4String>> reactions;
	reactions.push_back({"R1",   "e_aq^-1","e_aq^-1"});
	reactions.push_back({"R2",   "e_aq^-1","H3O^1"});
	reactions.push_back({"R3",   "e_aq^-1","H^0"});
	reactions.push_back({"R4",   "e_aq^-1","OH^0"});
	reactions.push_back({"R5",   "e_aq^-1","H2O2^0"});
	reactions.push_back({"R6",   "H3O^1",  "OH^-1"});
	reactions.push_back({"R7",   "H^0",    "H^0"});
	reactions.push_back({"R8",   "H^0",    "OH^0"});
	reactions.push_back({"R9",   "H^0",    "H2O2^0"});
	reactions.push_back({"R10",  "OH^0",   "OH^0"});
	for (size_t i = 0 ; i < fReactions.size() ; i++) {
		G4String ReactA = fMoleculesName[fReactions[i].reactorA];
		G4String ReactB = fMoleculesName[fReactions[i].reactorB];
		//G4bool isReaction = false;
		
		for (size_t j = 0 ; j < reactions.size() ; j++) {
			if (((ReactA == reactions[j][1]) && (ReactB == reactions[j][2]))
				||  ((ReactA == reactions[j][2]) && (ReactB == reactions[j][1]))) {
				//isReaction = true;
				fReactions[i].kobs /= fPm->GetUnitValue("/M/s");
				fReactions[i].kact /= fPm->GetUnitValue("/M/s");
				fReactions[i].kdif /= fPm->GetUnitValue("/M/s");
				fReactions[i].scavengingCapacity /= 1.0/s;
				
				if (reactions[j][0] == "R1") {
					/*fReactions[i].kobs = 2.0 * pow(10,12.281-
							         3.768e2/temperatureInKelvin-
 							         6.673e4/pow(temperatureInKelvin,2)-
							         1.075e7/pow(temperatureInKelvin,3));*/
					if (fKick)
						fReactions[i].kobs = 3.431e+13 * exp(-2296.32/temperatureInKelvin)*0.81;
					else
						fReactions[i].kobs = 3.431e+13 * exp(-2296.32/temperatureInKelvin);
                                } 
                                else if ( reactions[j][0] == "R2" ) {
					/*fReactions[i].kobs = pow(10,39.127-
							         3.888e4/temperatureInKelvin+
								 2.054e7/pow(temperatureInKelvin,2)-
								 4.899e9/pow(temperatureInKelvin,3)+
								 4.376e11/pow(temperatureInKelvin,4));*/
				
					if (fKick)
						fReactions[i].kobs = 1.693e+12 * exp(-1281.58/temperatureInKelvin)*1.08;
					else
						fReactions[i].kobs = 1.693e+12 * exp(-1281.58/temperatureInKelvin);
				
				}
				else if ( reactions[j][0] == "R3" ) {
					if (fKick)
						fReactions[i].kobs = 7.837e+12 * exp(-1681.53/temperatureInKelvin)*0.99;
					else
						fReactions[i].kobs = 7.837e+12 * exp(-1681.53/temperatureInKelvin);
				}
				else if ( reactions[j][0] == "R4" ) {
					/*fReactions[i].kobs = pow(10,13.123-
								 1.023e3/temperatureInKelvin+
								 7.634e4/pow(temperatureInKelvin,2));*/
					if (fKick)
						fReactions[i].kobs = 2.728e+12 * exp(-1300.38/temperatureInKelvin)*0.93;
					else
						fReactions[i].kobs = 2.728e+12 * exp(-1300.38/temperatureInKelvin);
				}
				else if ( reactions[j][0] == "R5" ) {
					if (fKick)
						fReactions[i].kobs = 8.316e+12 * exp(-1922.07/temperatureInKelvin)*0.93;
					else
						fReactions[i].kobs = 8.316e+12 * exp(-1922.07/temperatureInKelvin);
				}
				else if ( reactions[j][0] == "R6" ) {
					/*fReactions[i].kobs = pow(10,20.934-
								 1.236e4/temperatureInKelvin+
								 6.364e6/pow(temperatureInKelvin,2)-
								 1.475e9/pow(temperatureInKelvin,3)+
								 1.237e11/pow(temperatureInKelvin,4));*/
					if (fKick)
						fReactions[i].kobs = 1.555e+13 * exp(-1429.76/temperatureInKelvin)*1.21;
					else
						fReactions[i].kobs = 1.555e+13 * exp(-1429.76/temperatureInKelvin);
				}
				else if ( reactions[j][0] == "R7" ) {
					if (fKick)
						fReactions[i].kobs = 4.531e+12 * exp(-1816.37/temperatureInKelvin)*1.69;
					else
						fReactions[i].kobs = 4.531e+12 * exp(-1816.37/temperatureInKelvin);
				}
				else if ( reactions[j][0] == "R8" ) {
					if (fKick)
						fReactions[i].kobs = 4.531e+11 * exp(-1109.71/temperatureInKelvin)*1.51;
					else
						fReactions[i].kobs = 4.531e+11 * exp(-1109.71/temperatureInKelvin);
				}
				else if ( reactions[j][0] == "R9" ) {
					if (fKick)
						fReactions[i].kobs = 1.968e+11 * exp(-2559.09/temperatureInKelvin)*2.84;
					else
						fReactions[i].kobs = 1.968e+11 * exp(-2559.09/temperatureInKelvin);
				}
				else if ( reactions[j][0] == "R10" ) {
					if (fKick)
						fReactions[i].kobs = 8.728e+10 * exp(-858.30/temperatureInKelvin)*1.18;
					else
						fReactions[i].kobs = 8.728e+10 * exp(-858.30/temperatureInKelvin);
					/*fReactions[i].kobs = pow(10, 8.054 + 
								 2.193e3/temperatureInKelvin - 
								 7.395e5/pow(temperatureInKelvin,2) + 
								 6.870e7/pow(temperatureInKelvin,3));*/
				}
				fReactions[i].kobs *= fPm->GetUnitValue("/M/s");
				fReactions[i].kact *= fPm->GetUnitValue("/M/s");
				fReactions[i].kdif *= fPm->GetUnitValue("/M/s");
				fReactions[i].scavengingCapacity *= 1.0/s;
			}
		}	
	}
}
/*	
void TsIRTConfiguration::AdjustReactionRateForTemperature(G4double temperatureInCelsius) {
	G4double temperatureInKelvin = 273.15 + temperatureInCelsius;
	
	//May change with Temperature or H2SO4 Concentration
	G4double Kw      = K_water(temperatureInCelsius);
	G4double Hp_Con  = sqrt(Kw);
	G4double OHm_Con = Kw / Hp_Con;
	
	G4double Hp_25  = 1.003E-7;
	G4double OHm_25 = 1.003E-7;
	
	G4double H2O_Con = lH2Ol(temperatureInCelsius);
	G4double H2O_25  = lH2Ol(25.);
	
	std::vector<std::vector<G4String>> HerveReactions;
	HerveReactions.push_back({"R1",   "e_aq^-1","OH^0"});
	HerveReactions.push_back({"R2-1", "e_aq^-1","H3O^1"});
	HerveReactions.push_back({"R2-2", "H^0",    "None"});
	HerveReactions.push_back({"R3",   "e_aq^-1","e_aq^-1"});
	HerveReactions.push_back({"R4",   "e_aq^-1","H^0"});
	HerveReactions.push_back({"R5",   "e_aq^-1","H2O2^0"});
	HerveReactions.push_back({"R6-1", "e_aq^-1","None"});
	HerveReactions.push_back({"R6-2", "H^0",    "OH^-1"});
	HerveReactions.push_back({"R7",   "OH^0",   "OH^0"});
	HerveReactions.push_back({"R8",   "H^0",    "OH^0"});
	HerveReactions.push_back({"R9-1", "OH^0",   "OH^-1"});
	HerveReactions.push_back({"R9-2", "O^-1",   "None"});
	HerveReactions.push_back({"R10",  "H^0",    "H^0"});
	HerveReactions.push_back({"R11",  "H3O^1",  "HO2^-1"});
	HerveReactions.push_back({"R12",  "H3O^1",  "OH^-1"});
	HerveReactions.push_back({"R13",  "H3O^1",  "O^-1"});
	HerveReactions.push_back({"R14",  "H^0",    "O3P^0"});
	HerveReactions.push_back({"R15",  "H^0",    "H2O2^0"});
	
	// DC: kobs fitted by a Smoluchowski or Debye equation (3 or 4 in Herve et al, 2000)
	// DC*, A and A*: Arrhenius fits to kobs.
	// P: polynomial fit to kobs.
	// N: Arrhenius fit to kact.
	//
	// DC and DC*: kobs = kdif
	// A and A*: get kact from Noyes equation
	// P: kobs from fitting data of Elliot and kact from Noyes equation
	// N: after having kact, then get kobs from Noyes eq.
	
	for (size_t i = 0 ; i < fReactions.size() ; i++) {
		G4String ReactA = fMoleculesName[fReactions[i].reactorA];
		G4String ReactB = fMoleculesName[fReactions[i].reactorB];
		G4bool isHerveReaction = false;
		
		for (size_t j = 0 ; j < HerveReactions.size() ; j++) {
			if (((ReactA == HerveReactions[j][1]) && (ReactB == HerveReactions[j][2]))
				||  ((ReactA == HerveReactions[j][2]) && (ReactB == HerveReactions[j][1]))) {
				isHerveReaction = true;
				fReactions[i].kobs /= fPm->GetUnitValue("/M/s");
				fReactions[i].kact /= fPm->GetUnitValue("/M/s");
				fReactions[i].kdif /= fPm->GetUnitValue("/M/s");
				fReactions[i].scavengingCapacity /= 1.0/s;
				
				if (HerveReactions[j][0] == "R1") {
					//if ( fKick )
					//	fReactions[i].kobs = 0.844 * pow(10, 13.123 - 1.023e3/ temperatureInKelvin + 7.634e4/(pow(temperatureInKelvin,2)));
					//else
					//	fReactions[i].kobs = pow(10, 13.123 - 1.023e3/ temperatureInKelvin + 7.634e4/(pow(temperatureInKelvin,2)));
					
					fReactions[i].kact = Arrhenius(0.0304e12, temperatureInKelvin,-3.5);//Noyes(fReactions[i].kobs, 0 , fReactions[i].kdif);
					if (fKick)
						fReactions[i].kact *= 1.071;
					
					//fReactions[i].kdif = Noyes(fReactions[i].kobs, fReactions[i].kact, 0);
					fReactions[i].kobs = Noyes(0, fReactions[i].kact, fReactions[i].kdif);
				}
				
				else if (HerveReactions[j][0] == "R2-1") {
					if (fReactions[i].reactionType != 6) {
						//if ( fKick )
						//	fReactions[i].kobs = 1.09 * pow(10, 39.127 - 3.888e4 / temperatureInKelvin + 2.054e7/(pow(temperatureInKelvin,2)) - 4.899e9/(pow(temperatureInKelvin,3)) + 4.376e11/(pow(temperatureInKelvin,4)));
						//else
						//fReactions[i].kobs = pow(10, 39.127 - 3.888e4 / temperatureInKelvin + 2.054e7/(pow(temperatureInKelvin,2)) - 4.899e9/(pow(temperatureInKelvin,3)) + 4.376e11/(pow(temperatureInKelvin,4)));
						fReactions[i].kobs = Arrhenius(1.24e12, temperatureInKelvin, 10.1 );
						if ( fKick )
							fReactions[i].kobs *= 1.095;
						
						//fReactions[i].kdif = Debye(1.0, fReactions[i], temperatureInKelvin);
						//if ( fKick)
						//	fReactions[i].kdif *= 0.98;
						
						fReactions[i].kact = Noyes(fReactions[i].kobs, 0 , fReactions[i].kdif);
					}
					else {
						fReactions[i].scavengingCapacity = 1.0004*Arrhenius(1.24E12, temperatureInKelvin, 10.1) * Hp_Con;
					}
				}
				
				else if (HerveReactions[j][0] == "R2-2") {
					fReactions[i].kobs = k_26(temperatureInCelsius);
				}
				
				else if (HerveReactions[j][0] == "R3") {
					if (temperatureInCelsius <= 150 ) {//2.33e13,20.3) From Elliot 1994
						//if ( fKick )
						//	fReactions[i].kobs = 2.0 * 0.76 * pow(10 , 12.281 - 3.768e2 / temperatureInKelvin - 6.673e4/(pow(temperatureInKelvin, 2)) - 1.075e7/(pow(temperatureInKelvin,3)));
						//else
						//	fReactions[i].kobs = 2.0 * pow(10 , 12.281 - 3.768e2 / temperatureInKelvin - 6.673e4/(pow(temperatureInKelvin, 2)) - 1.075e7/(pow(temperatureInKelvin,3)));
						fReactions[i].kobs = 2.0 * Arrhenius(2.33e13, temperatureInKelvin, 20.3);
						if ( fKick )
							fReactions[i].kobs *= 8.50E-01;
						
						fReactions[i].kdif = fReactions[i].kobs;
					}
					else {
						fReactions[i].kobs = pow(10 , - 47.532 + 4.920e4 / temperatureInKelvin - 1.036e7/(pow(temperatureInKelvin,2)));
						fReactions[i].kdif = fReactions[i].kobs;
					}
				}
				
				else if (HerveReactions[j][0] == "R4") {
					if ( fKick )
						fReactions[i].kobs = 0.90 * Arrhenius(1.14e13, temperatureInKelvin, 14.93052722) ;
					else
						fReactions[i].kobs =  Arrhenius(1.14e13, temperatureInKelvin, 14.93052722);
					fReactions[i].kdif = fReactions[i].kobs;
				}
				
				else if (HerveReactions[j][0] == "R5") {
					if ( fKick ) {
						//fReactions[i].kact = 0.6961*Arrhenius(6.26E12, temperatureInKelvin, 14.0);
						fReactions[i].kobs = 8.08E-01 * Arrhenius(7.70e12, temperatureInKelvin, 15.71126816);
					}
					else {
						//fReactions[i].kact = Arrhenius(6.26E12, temperatureInKelvin, 14.0);//0.697;
						fReactions[i].kobs = Arrhenius(7.70e12, temperatureInKelvin, 15.71126816);
					}
					//fReactions[i].kobs = Noyes(0, fReactions[i].kact, fReactions[i].kdif);
					//fReactions[i].kdif = Noyes(fReactions[i].kobs, fReactions[i].kact, 0);
					fReactions[i].kact = Noyes(fReactions[i].kobs, 0 , fReactions[i].kdif);
				}
				
				else if (HerveReactions[j][0] == "R6-1") {
					fReactions[i].kobs = Arrhenius(133E12, temperatureInKelvin, 38.38) * K_water(temperatureInCelsius) /
					(lH2Ol(temperatureInCelsius) * K_H(temperatureInCelsius));
				}
				
				else if (HerveReactions[j][0] == "R6-2") {
					if (fReactions[i].reactionType != 6) {
						if ( fKick )
							fReactions[i].kobs = 0.9992*Arrhenius(133.0E12, temperatureInKelvin, 38.38); // Could Be Eo = 38.8
						else
							fReactions[i].kobs = Arrhenius(133.0E12, temperatureInKelvin, 38.38); // Could Be Eo = 38.8
						fReactions[i].kact = Noyes(fReactions[i].kobs, 0 , fReactions[i].kdif);
					}
					
					else {
						fReactions[i].scavengingCapacity = 0.9992*Arrhenius(133.0E12, temperatureInKelvin, 38.38) * OHm_Con;
					}
				}
				
				else if (HerveReactions[j][0] == "R7") {
					//if ( fKick )
					//	fReactions[i].kobs = 1.14 * pow(10 , 8.054 + 2.193E3 / temperatureInKelvin - 7.395e5/(pow(temperatureInKelvin,2)) + 6.870e7/(pow(temperatureInKelvin,3)));
					//else
					//	fReactions[i].kobs = pow(10 , 8.054 + 2.193E3 / temperatureInKelvin - 7.395e5/(pow(temperatureInKelvin,2)) + 6.870e7/(pow(temperatureInKelvin,3)));
					//fReactions[i].kobs = Arrhenius(1.04e11, temperatureInKelvin, 7.65);
					
					//fReactions[i].kdif = Smoluchowski(1.0, fReactions[i]);
					fReactions[i].kact = 2.0 * Arrhenius(0.0369e12, temperatureInKelvin, 3.0);
					if ( fKick )
						fReactions[i].kdif *= 1.01E+00;
					
					fReactions[i].kobs = Noyes(0, fReactions[i].kact, fReactions[i].kdif);
				}
				
				else if (HerveReactions[j][0] == "R8") {
					if ( fKick )
						fReactions[i].kobs = 1.83 * Arrhenius(4.26e11, temperatureInKelvin, 9.07871174);
					else
						fReactions[i].kobs = Arrhenius(4.26e11, temperatureInKelvin, 9.07871174);
					
					fReactions[i].kact = Arrhenius(0.178e12, temperatureInKelvin, 4.5);
					if (fKick)
						fReactions[i].kact *= 2.305;
					//fReactions[i].kact = Noyes(fReactions[i].kobs, 0 ,fReactions[i].kdif);
				}
				
				else if (HerveReactions[j][0] == "R9-1") {
					if (fReactions[i].reactionType != 6) {
						G4double invTemp = 1./(temperatureInCelsius+273.15);
						G4double kobs = -28232.6 * std::pow(invTemp, 1.72145) +  11.66062;
						G4double kact = -19846.6 * std::pow(invTemp, 1.65751) +  12.06992;
						kobs = std::pow(10., kobs);
						kact = std::pow(10., kact);
						if ( fKick ) {
							fReactions[i].kobs = 0.4914197*kobs;
							fReactions[i].kact = 0.2504124*kact;
						}
						else {
							fReactions[i].kobs = kobs;
							fReactions[i].kact = kact;
						}
						fReactions[i].kdif = Noyes(kobs, kact, 0);
					}
					
					else {
						G4double invTemp = 1./(temperatureInCelsius+273.15);
						G4double kobs = -28232.6 * std::pow(invTemp, 1.72145) +  11.66062;
						kobs = std::pow(10., kobs);
						if ( fKick )
							fReactions[i].scavengingCapacity = 0.4914197 * kobs * OHm_Con;
						else
							fReactions[i].scavengingCapacity = kobs * OHm_Con;
					}
				}
				
				else if (HerveReactions[j][0] == "R9-2") {
					fReactions[i].kobs = k29_inverse(temperatureInCelsius);
				}
				
				else if (HerveReactions[j][0] == "R10") {
					if ( fKick )
						fReactions[i].kobs = 2 * 1.51 * Arrhenius(2.70e12, temperatureInKelvin, 15.5275155);
					else
						fReactions[i].kobs = 2 * Arrhenius(2.70e12, temperatureInKelvin, 15.5275155);
					fReactions[i].kdif = fReactions[i].kobs;
				}
				
				else if (HerveReactions[j][0] == "R11") {
					if (fReactions[i].reactionType != 6) {
						G4double invTemp = 1./(temperatureInCelsius+273.15);
						G4double kobs = -319.8 * std::pow(invTemp, 0.85502) +  13.14262;
						G4double kact = -651.2 * std::pow(invTemp, 0.98174) +  13.53771;
						kobs = std::pow(10., kobs);
						kact = std::pow(10., kact);
						fReactions[i].kobs = kobs;
						fReactions[i].kact = kact;
						fReactions[i].kdif = Noyes(kobs, kact, 0);
					}
					else {
						G4double invTemp = 1./(temperatureInCelsius+273.15);
						G4double kobs = -319.8 * std::pow(invTemp, 0.85502) +  13.14262;
						kobs = std::pow(10., kobs);
						fReactions[i].scavengingCapacity = kobs * Hp_Con;
					}
				}
				
				else if (HerveReactions[j][0] == "R12") {
					if (fReactions[i].reactionType != 6) {
						//if ( fKick )
						//	fReactions[i].kobs = 1.21 * pow(10,20.934 - 1.236e4/temperatureInKelvin + 6.364e6/(pow(temperatureInKelvin,2)) - 1.475e9/(pow(temperatureInKelvin,3)) + 1.237e11/(pow(temperatureInKelvin,4)));
						//else
						//	fReactions[i].kobs = pow(10,20.934 - 1.236e4/temperatureInKelvin + 6.364e6/(pow(temperatureInKelvin,2)) - 1.475e9/(pow(temperatureInKelvin,3)) + 1.237e11/(pow(temperatureInKelvin,4)));
						fReactions[i].kobs = Debye(1, fReactions[i], temperatureInCelsius);
						if ( fKick)
							fReactions[i].kobs *= 1.24E+00;
						fReactions[i].kdif = fReactions[i].kobs;
					}
					
					else {
						if ( fKick )
							fReactions[i].scavengingCapacity = 0.9766 * Debye(1, fReactions[i], temperatureInCelsius) * OHm_Con;
						else
							fReactions[i].scavengingCapacity = Debye(1, fReactions[i], temperatureInCelsius) * OHm_Con;
					}
				}
				
				else if (HerveReactions[j][0] == "R13") {
					if (fReactions[i].reactionType != 6) {
						G4double invTemp = 1./(temperatureInCelsius+273.15);
						G4double kobs = -319.8 * std::pow(invTemp, 0.85502) + 13.14262;
						G4double kact = -762.0 * std::pow(invTemp, 1.02866) + 13.25234;
						kobs = std::pow(10., kobs);
						kact = std::pow(10., kact);
						fReactions[i].kobs = kobs;
						fReactions[i].kact = kact;
						fReactions[i].kdif = Noyes(kobs, kact, 0);
					}
					
					else {
						G4double invTemp = 1./(temperatureInCelsius+273.15);
						G4double kobs = -319.8 * std::pow(invTemp, 0.85502) + 13.14262;
						kobs = std::pow(10., kobs);
						fReactions[i].scavengingCapacity = kobs * Hp_Con;
					}
				}
				
				else if (HerveReactions[j][0] == "R14") {
					//R14: DC* (13)
					fReactions[i].kobs = Arrhenius(10.9E12, temperatureInKelvin, 15.6);
					fReactions[i].kdif = fReactions[i].kobs;
				}
				else if (HerveReactions[j][0] == "R15") {
					
					if ( fKick )
						fReactions[i].kobs = 2.47 * Arrhenius(1.79e11, temperatureInKelvin, 21.06587056);
					else
						fReactions[i].kobs = Arrhenius(1.79e11, temperatureInKelvin, 21.06587056);
					//fReactions[i].kdif = Smoluchowski(1, fReactions[i]);
					fReactions[i].kact = Noyes(fReactions[i].kobs, 0 , fReactions[i].kdif);
				}
				fReactions[i].kobs *= fPm->GetUnitValue("/M/s");
				fReactions[i].kact *= fPm->GetUnitValue("/M/s");
				fReactions[i].kdif *= fPm->GetUnitValue("/M/s");
				fReactions[i].scavengingCapacity *= 1.0/s;
			}
		}
		
		if (!isHerveReaction) {
			if (fReactions[i].reactionType == 6) {
				if (ReactB == "H3O^1") {
					fReactions[i].scavengingCapacity = fReactions[i].scavengingCapacity * Hp_Con / Hp_25;
				}
				else if (ReactB == "OH^-1") {
					fReactions[i].scavengingCapacity = fReactions[i].scavengingCapacity * OHm_Con / OHm_25;
				}
				else if (ReactB == "None") {
					fReactions[i].scavengingCapacity = fReactions[i].scavengingCapacity * H2O_Con / H2O_25;
				}
			}
		}
	}
}
*/

void TsIRTConfiguration::AdjustReactionRateForPH(G4String pHOrConcentration) {
	std::vector<G4double> AcidComponents;
	G4double Ionic   = 0.0;
	G4double HCon    = 0.0;
	G4double HCon25  = 1.00E-7;
	G4double OHCon   = 0.0;
	G4double OHCon25 = 1.00E-7;
	G4double HSO4Con = 0.0;
	
	if (pHOrConcentration == "PH" && fpHSolvent == "h2so4") {
		AcidComponents = GetH2SO4ComponentsConcentrationPH(fpHValue);
		Ionic          = GetIonicStrength(AcidComponents);
		HCon           = AcidComponents[0];
		HSO4Con        = AcidComponents[1];
		OHCon          = AcidComponents[3];
		G4cout << "-- Adjust for PH of H2SO4 " << G4endl;
		
	}
	
	else if (pHOrConcentration == "Concentration" && fpHSolvent == "h2so4") {
		AcidComponents = GetH2SO4ComponentsConcentrationP(fpHSolventConcentration/fPm->GetUnitValue("M"));
		Ionic          = GetIonicStrength(AcidComponents);
		HCon           = AcidComponents[0];
		HSO4Con        = AcidComponents[1];
		OHCon          = AcidComponents[3];
		G4cout << "-- Adjust for Concentration of H2SO4 " << G4endl;
	}
	
	else if (fpHSolvent == "generic") {
		HCon  = pow(10,-fpHValue);
		OHCon = 1E-14 / HCon;
		AcidComponents = {HCon, 0.0, 0.0, OHCon, 0.0, 0.0};
		Ionic          = GetIonicStrength(AcidComponents);
		G4cout << "-- Adjust for a generic substance " << G4endl;
		
	} else {
		G4cout << "-- Is doing nothing " << pHOrConcentration << " " << fpHSolvent << G4endl;
	}
	/*if ( fpHSolvent == "h2so4" ) {
	 if ( pHOrConcentration == "PH" ) {
	 AcidComponents = GetH2SO4ComponentsConcentrationPH(fpHValue);
	 } else {
	 AcidComponents = GetH2SO4ComponentsConcentrationP(fpHSolventConcentration);
	 }
	 
	 Ionic          = GetIonicStrength(AcidComponents);
	 HCon           = AcidComponents[0];
	 HSO4Con        = AcidComponents[1];
	 OHCon          = AcidComponents[3];
	 } else {
	 HCon  = pow(10,-fpHValue);
	 OHCon = 1E-14 / HCon;
	 AcidComponents = {HCon, 0.0, 0.0, OHCon, 0.0, 0.0};
	 Ionic          = GetIonicStrength(AcidComponents);
	 }*/
	
	G4cout << G4endl;
	G4cout << " ###-------- pH Scaling Starts ---------###" << G4endl;
	
	for(size_t i = 0; i < fReactions.size(); i++) {
		G4int chargeA   = fMoleculesDefinition[fReactions[i].reactorA].charge;
		G4int chargeB   = fMoleculesDefinition[fReactions[i].reactorB].charge;
		G4String ReactA = fMoleculesName[fReactions[i].reactorA];
		G4String ReactB = fMoleculesName[fReactions[i].reactorB];
		std::vector<G4String> products;
		G4double k_Before = 0;
		
		for (size_t vsize = 0; vsize < fReactions[i].products.size(); vsize++){
			products.push_back(fMoleculesName[fReactions[i].products[vsize]]);
		}
		
		fReactions[i].kobs /= fPm->GetUnitValue("/M/s");
		fReactions[i].scavengingCapacity /= 1/s;
		
		if ((chargeA != 0 && chargeB != 0) && (!(ReactA == "e_aq^-1" && ReactB == "e_aq^-1"))) {
			if (fReactions[i].reactionType != 6) {
				k_Before = fReactions[i].kobs;
				fReactions[i].kobs = IonicRate(Ionic, fReactions[i]);
			}
			
			else if (fReactions[i].reactionType == 6 && ReactB == "H3O^1") {
				k_Before = fReactions[i].scavengingCapacity;
				fReactions[i].scavengingCapacity = IonicRate(Ionic, fReactions[i]) * HCon;
			}
			
			else if (fReactions[i].reactionType == 6 && ReactB == "OH^-1") {
				k_Before = fReactions[i].scavengingCapacity;
				fReactions[i].scavengingCapacity = IonicRate(Ionic, fReactions[i]) * OHCon;
			}
			
			else {
				G4cout << "========= "<< fReactions[i].reactionType << G4endl;
				k_Before = fReactions[i].kobs;
				fReactions[i].kobs = fReactions[i].kobs * 1E-7;
				fReactions[i].kobs = IonicRate(Ionic, fReactions[i]);
			}
		}
		
		else if ((fReactions[i].reactionType == 6) && (ReactB == "H3O^1")) {
			k_Before = fReactions[i].scavengingCapacity;
			fReactions[i].scavengingCapacity = fReactions[i].scavengingCapacity * HCon / HCon25;
		}
		
		else if ((fReactions[i].reactionType == 6) && (ReactB == "OH^-1")) {
			k_Before = fReactions[i].scavengingCapacity;
			fReactions[i].scavengingCapacity = fReactions[i].scavengingCapacity * OHCon / OHCon25;
		}
		
		else if ((fReactions[i].reactionType == 6) && (ReactB == "HSO4^-1")) {
			k_Before = fReactions[i].scavengingCapacity;
			fReactions[i].scavengingCapacity = fReactions[i].scavengingCapacity * HSO4Con;
		}
		
		if (k_Before > 0) {
			G4cout << " Reaction Type: " << fReactions[i].reactionType << " | Reaction: " << ReactA << " + " << ReactB;
			for (int prod = 0; prod < (int)products.size(); prod++) {
				if (prod == 0)
					G4cout << " -> ";
				G4cout << products[prod];
				if (prod < (int)products.size() - 1) {
					G4cout << " + ";
				}
			}
			if ( fReactions[i].reactionType == 6 )
				G4cout << G4endl << " ---- scav: " << k_Before << " ---> "  << fReactions[i].scavengingCapacity << G4endl << G4endl;
			else
				G4cout << G4endl << " ---- kobs: " << k_Before << " ---> "  << fReactions[i].kobs << G4endl << G4endl;
		}
		fReactions[i].kobs *= fPm->GetUnitValue("/M/s");
		fReactions[i].scavengingCapacity *= 1/s;
	}
	
	G4cout << " ###-------- pH Scaling Ends ---------###" << G4endl;
	G4cout << G4endl;
}


void TsIRTConfiguration::ResampleReactantsPosition(TsMolecule& molA, TsMolecule& molB, G4int index, G4double time) {
	// Position approach
	G4double D1 = fMoleculesDefinition[molA.id].diffusionCoefficient;
	G4double D2 = fMoleculesDefinition[molB.id].diffusionCoefficient;
	
	G4ThreeVector r1 = molA.position;
	G4ThreeVector r2 = molB.position;
	
	if ( D1 == 0 ) {
		molB.position = r1;
		return;
	} else if ( D2 == 0 ) {
		molA.position = r2;
		return;
	}
	
	G4ThreeVector S1 = r1-r2;
	G4double r0 = S1.mag();
	G4double effectiveReactionRadius = fReactions[index].effectiveReactionRadius;
	if ( fReactions[index].reactionType == 4 )
		effectiveReactionRadius = fReactions[index].effectiveTildeReactionRadius;
	
	S1.setMag(effectiveReactionRadius);
	G4double dt = std::fabs(time-molA.time);
	
	G4double s12 = 2.0 * D1 * dt;
	G4double s22 = 2.0 * D2 * dt;
	G4double alpha = effectiveReactionRadius * r0/(2*(D1+D2)*dt);
	
	G4ThreeVector S2 = (r1 + (s12/s22)*r2) + G4ThreeVector(G4RandGauss::shoot(0.0, s12 + s22*s22/s12),
														   G4RandGauss::shoot(0.0, s12 + s22*s22/s12),
														   G4RandGauss::shoot(0.0, s12 + s22*s22/s12));
	
	S1.setPhi(rad*G4UniformRand()*2.0*CLHEP::pi);
	S1.setTheta(rad*std::acos(1.0 + 1./alpha * std::log(1.0 - G4UniformRand()*(1.-std::exp(-2.0*alpha)))));
	
	G4ThreeVector R1 = (D1 * S1 + D2*S2)/(D1+D2);
	G4ThreeVector R2 = D2*(S2-S1)/(D1+D2);
	
	molA.position = R1;
	molA.time = time;
	molB.position = R2;
	molB.time = time;
	
}


std::vector<G4ThreeVector> TsIRTConfiguration::GetPositionOfProducts(TsMolecule molA, TsMolecule molB, G4int index) {
	G4double D1 = fMoleculesDefinition[molA.id].diffusionCoefficient;
	G4double D2 = fMoleculesDefinition[molB.id].diffusionCoefficient;
	G4ThreeVector r1 = molA.position;
	G4ThreeVector r2 = molB.position;
	
	std::vector<G4ThreeVector> result;
	// weighted-position
	G4ThreeVector position = r1*std::sqrt(D2)/(std::sqrt(D1) + std::sqrt(D2)) + r2*std::sqrt(D1)/(std::sqrt(D1) + std::sqrt(D2));
	if ( fReactions[index].products.size() == 1 ) {
		// At weighted position.
		//result.push_back(position);
		if ( G4UniformRand() > 0.5 )
			result.push_back(r1);
		else
			result.push_back(r2);
		
	} else if ( fReactions[index].products.size() == 3 ) {
		// at weighted and at parents position
		result.push_back(position);
		result.push_back(r1);
		result.push_back(r2);
	} else if ( fReactions[index].products.size() == 2 ) {
		// at parents position
		result.push_back(r1);
		result.push_back(r2);
	}
	
	return result;
	
}

G4double TsIRTConfiguration::GetRCutOff(G4double tCutOff) {
	G4double probabilityOfReaction = 0.01;
	G4double maximumReactionRadius = 1.45*nm;
	G4double maximumRelativeDiffusionCoefficient = 2.0*9.46e9 *nm*nm/s;
	G4double erfcInv = fUtils->erfcInv(probabilityOfReaction);
	return maximumReactionRadius + 2.0 * std::sqrt(maximumRelativeDiffusionCoefficient * tCutOff) * erfcInv;
}


G4double TsIRTConfiguration::GetIndependentReactionTime(TsMolecule molA, TsMolecule molB, G4int indexOfReaction) {
	G4int typeOfReaction = fReactions[indexOfReaction].reactionType;
	G4double result = -1.0*ps;
	
	if ( typeOfReaction == 1 || typeOfReaction == 3 || typeOfReaction == 5) {
		result = SampleIRTTotallyDiffusionControlled(molA, molB, indexOfReaction);
	} else if ( typeOfReaction == 2 || typeOfReaction == 4 ) {
		result = SampleIRTPartiallyDiffusionControlled(molA, molB, indexOfReaction);
	}
	
	return result;
}


G4double TsIRTConfiguration::SampleIRTTotallyDiffusionControlled(TsMolecule molA, TsMolecule molB, G4int indexOfReaction) {
	G4int typeOfReaction = fReactions[indexOfReaction].reactionType;
	G4double probFactor = 1.0;
	G4double reff = (molA.position - molB.position).mag();
	
	G4double Reff = fReactions[indexOfReaction].effectiveReactionRadius;
	if ( typeOfReaction == 4 ) // Only if considering  all reactions totally diffusion controlled (for testing purposes)
		Reff = fReactions[indexOfReaction].effectiveTildeReactionRadius;
	
	if ( fReactions[indexOfReaction].OnsagerRadius != 0 ) {
		G4double rc = fReactions[indexOfReaction].OnsagerRadius;
		reff = -rc/(1 - std::exp(rc/reff));
	}
	
	G4double Winf = probFactor * Reff/reff;
	G4double D = fMoleculesDefinition[molA.id].diffusionCoefficient + fMoleculesDefinition[molB.id].diffusionCoefficient;
	
	G4double W = G4UniformRand();
	G4double irt = -1.0 * ps;
	
	if ( W < Winf ) {
		irt = (0.25/D) * std::pow( (reff-Reff)/fUtils->erfcInv(reff*W/Reff), 2 );
	}
	return irt;
}


void TsIRTConfiguration::TestSampling(G4int indexOfReaction, G4int nHistories) {
	G4double sigma = 1.0 * nm;
	G4int molA = -1, molB = -1;
	for ( auto& reactions : fReactions ) {
		if ( reactions.second.index == indexOfReaction ) {
			molA = reactions.second.reactorA;
			molB = reactions.second.reactorA;
			break;
		}
	}
	G4ThreeVector posA;
	G4ThreeVector posB;
	std::ofstream out("TestSampling.csv");
	out << "# Reaction ID " << indexOfReaction << std::endl;
	for ( int i = 0; i < nHistories; i++ ) {
		TsMolecule A;
		A.id = molA;
		A.time = 0.01*ps;
		A.reacted = false;
		A.trackID = 0;
		A.spin = -1;
		posA = G4ThreeVector(
							 G4RandGauss::shoot(0., 1.0)*sigma,
							 G4RandGauss::shoot(0., 1.0)*sigma,
							 G4RandGauss::shoot(0., 1.0)*sigma);
		A.position = posA;
		
		TsMolecule B;
		B.id = molB;
		B.time = 0.01*ps;
		B.reacted = false;
		B.trackID = 0;
		B.spin = -1;
		posB = G4ThreeVector(
							 G4RandGauss::shoot(0., 1.0)*sigma,
							 G4RandGauss::shoot(0., 1.0)*sigma,
							 G4RandGauss::shoot(0., 1.0)*sigma);
		B.position = posB;
		
		G4double irt = GetIndependentReactionTime(A, B, indexOfReaction);
		if ( irt > 0.0 )
			out << (irt)/ps << std::endl;
	}
	out.close();
}


G4double TsIRTConfiguration::SampleIRTPartiallyDiffusionControlled(TsMolecule molA, TsMolecule molB, G4int indexOfReaction) {
	G4double r0 = (molA.position-molB.position).mag();
	
	G4double irt = 0.0;
	//G4double prob = 0.0;
	if ( fReactions[indexOfReaction].reactionType == 2 ) {
		// From paper Plante 2017, Considerations for ...
		G4double alpha = -(fReactions[indexOfReaction]).alpha;
		G4double sigma = (fReactions[indexOfReaction]).effectiveReactionRadius;
		G4double D = fMoleculesDefinition[molA.id].diffusionCoefficient + fMoleculesDefinition[molB.id].diffusionCoefficient;
		irt = fUtils->SampleTypeII(alpha, sigma, r0, D);
		
	} else {
		G4double rc = (fReactions[indexOfReaction]).OnsagerRadius;
		G4double reff = -rc/(1 - std::exp(rc/r0));
		G4double sigmaEffEff = (fReactions[indexOfReaction]).effectiveTildeReactionRadius;
		
		G4double Winf = sigmaEffEff/reff;
		G4double W = G4UniformRand();
		if ( W < Winf ) {
			irt = brents_fun(molA, molB, indexOfReaction, -W);
			//prob = W;
		} else {
			irt = -1.0*ps;
		}
	}
	return irt;
}


G4double TsIRTConfiguration::SolveTime(TsMolecule, G4int indexOfReaction, G4double offset) {
	G4double Reff = (fReactions[indexOfReaction]).effectiveReactionRadius;
	G4double CsB = (fReactions[indexOfReaction]).probabilityOfReaction/s;
	G4double nm3persToPerMPers = 6.022140857e-1;
	G4double k = (fReactions[indexOfReaction]).kobs/nm3persToPerMPers * (nm*nm*nm/s);
	
	G4double sqrtArgument = std::pow(4.0*Reff*std::sqrt(Reff/k),2) - 4.0/CsB * std::log(1.0 - offset);
	if (sqrtArgument < 0 )
		return -1.0*ps;
	G4double t1 = std::pow(-2.0*Reff*sqrt(Reff/k) + 0.5 * sqrtArgument, 2);
	G4double t2 = std::pow(-2.0*Reff*sqrt(Reff/k) - 0.5 * sqrtArgument, 2);
	
	if ( t1 <= t2 )
		return t1;
	else
		return t2;
}


G4int  TsIRTConfiguration::ContactFirstOrderAndBackgroundReactions(TsMolecule molA) {
	G4int pdgA = molA.id;
	std::vector<size_t> index;
	for ( size_t u = fTotalBinaryReaction; u < fReactions.size(); u++ )
		index.push_back(u);
	
	size_t sizeIndex = index.size();
	if ( 0 < sizeIndex )
		std::random_shuffle(index.begin(), index.end());
	
	for ( size_t v = 0; v < sizeIndex; v++ ) {
		size_t u = index[v];
		if (pdgA == fReactions[u].reactorA) {
			G4double R3 = fMoleculesDefinition[pdgA].radius*fMoleculesDefinition[pdgA].radius*fMoleculesDefinition[pdgA].radius;
			G4double Cs = fReactions[u].concentration/fPm->GetUnitValue("M");
			// nm3 to M multiply by 10^-24 Nav = 6.02214076×10^23x10^-24 = 0.602214076
			Cs /= 6.022140857e-1 * (nm*nm*nm);
			G4double prob = std::exp(-4.0*CLHEP::pi*R3*Cs/3.);
			if ( G4UniformRand() < 1. - prob ) {
				return (int)u;
			}
		}
	}
	
	return -1;
}


G4double TsIRTConfiguration::SampleExponentialTime(G4int pdgA, G4int pdgB, G4int indexOfReaction) {
	G4double D = fMoleculesDefinition[pdgA].diffusionCoefficient + fMoleculesDefinition[pdgB].diffusionCoefficient;
	D /= nm*nm/s;
	G4double Rreact = fReactions[indexOfReaction].effectiveReactionRadius;
	Rreact /= nm;
	G4double Cs = fReactions[indexOfReaction].concentration/0.60221/fPm->GetUnitValue("M");
	G4double A = -Rreact/std::sqrt(CLHEP::pi * D);
	G4double B = 0.5/std::sqrt(CLHEP::pi * D);
	G4double C = 4.0 * Rreact - std::log(1.0 - G4UniformRand())/(Rreact*Cs);
	G4double timeA = A + B * std::sqrt(C);
	G4double timeB = A - B * std::sqrt(C);
	if ( timeA < timeB )
		return timeA*timeA*s;
	else
		return timeB*timeB*s;
}


std::vector<std::pair<G4int, G4double>> TsIRTConfiguration::SampleAllIRTFirstOrderAndBackgroundReactions(TsMolecule molA ) {
	std::vector<std::pair<G4int, G4double>> result;
	G4int pdgA = molA.id;
	G4double scavengingCapacity, prob, time;
	
	for ( size_t u = fTotalBinaryReaction; u < fReactions.size(); u++ ) { // TODO-> Review this algorithm: use t or dt???
		if ( fReactions[u].reactionType == 6 && pdgA == fReactions[u].reactorA) {
			scavengingCapacity = fReactions[u].scavengingCapacity;
			prob = G4UniformRand();
			if ( fReactions[u].sampleExponential )
				time = SampleExponentialTime(pdgA,fReactions[u].reactorB, G4int(u));
			else
				time = -(std::log(1.0 - prob)/scavengingCapacity);
			
			if ( time > 0.0 )
				result.push_back(std::make_pair((int)u, time));
		}
	}
	return result;
}


std::pair<G4int, G4double> TsIRTConfiguration::SampleIRTFirstOrderAndBackgroundReactions(TsMolecule molA ) {
	G4int pdgA = molA.id;
	// Search for the IRT with the minium value
	G4double irt = 1000.0*s;
	
	G4int index = -1;
	G4bool found = false;
	
	G4double scavengingCapacity, prob, time;
	
	for ( size_t u = fTotalBinaryReaction; u < fReactions.size(); u++ ) { // TODO-> Review this algorithm: use t or dt???
		if ( fReactions[u].reactionType == 6 && pdgA == fReactions[u].reactorA) {
			scavengingCapacity = fReactions[u].scavengingCapacity;
			prob = G4UniformRand();
			if ( fReactions[u].sampleExponential )
				time = SampleExponentialTime(pdgA,fReactions[u].reactorB, G4int(u));
			else
				time = -(std::log(1.0 - prob)/scavengingCapacity);
			
			if ( time < irt && time > 0.0) {
				irt = time;
				index = (int)u;
			}
			found = true;
		}
	}
	
	if ( !found )
		return std::make_pair(-1,-1*ps);
	
	return std::make_pair(index, irt);
}


G4double TsIRTConfiguration::CalculateProbabilityPartiallyDiffusionControlled(TsMolecule molA, TsMolecule molB, G4int indexOfReaction, G4double t) {
	G4double DA = fMoleculesDefinition[molA.id].diffusionCoefficient;
	G4double DB = fMoleculesDefinition[molB.id].diffusionCoefficient;
	G4double D = DA + DB;
	G4double r0 = (molA.position-molB.position).mag();
	
	if ( fReactions[indexOfReaction].reactionType == 2 ) { // between neutral particles
		// From paper A Monte Carlo step-by-step
		G4double sigma = (fReactions[indexOfReaction]).effectiveReactionRadius;
		G4double kact = (fReactions[indexOfReaction]).kact;
		
		G4double alpha =(fReactions[indexOfReaction]).alpha;
		G4double factor = kact/(4.0*CLHEP::pi*sigma*D*r0*alpha);
		G4double x = (r0-sigma)/std::sqrt(4.0*D*t);
		G4double y = alpha*std::sqrt(D*t);
		
		return factor * (fUtils->erfc(x) - std::exp(-x*x) * fUtils->erfcx(x+y));
		
	} else {
		G4double rc = (fReactions[indexOfReaction]).OnsagerRadius;
		G4double reff = -rc/(1 - std::exp(rc/r0));
		G4double sigma = (fReactions[indexOfReaction]).reactionRadius;
		G4double sigmaEffEff = (fReactions[indexOfReaction]).effectiveTildeReactionRadius;
		
		G4double alpha = (fReactions[indexOfReaction]).alpha;
		G4double a = (4.0*sigma*sigma*alpha/(rc*rc))*std::sqrt(t/D)*std::sinh(0.5*rc/sigma)*std::sinh(0.5*rc/sigma);
		G4double b = (0.25*rc/std::sqrt(D*t))*(std::cosh(0.5*rc/r0)/std::sinh(0.5*rc/r0) - std::cosh(0.5*rc/sigma)/std::sinh(0.5*rc/sigma));
		
		G4double factor = sigmaEffEff/reff;
		
		return factor * (fUtils->erfc(b) - std::exp(-b*b) * fUtils->erfcx(a+b));
	}
}


void TsIRTConfiguration::Quit(const G4String& name, G4String message) {
	G4cerr << G4endl;
	G4cerr << "Topas is exiting due to a serious error in Chemistry IRT setup." << G4endl;
	G4cerr << "--- Parameter name: " << name << G4endl;
	G4cerr << "--- " << message << G4endl;
	G4cerr << G4endl;
	fPm->AbortSession(1);
}


void TsIRTConfiguration::SetTimeLimits(G4double lower, G4double upper) {
	fUpperTime = upper;
	fLowerTime = lower;
}


G4int TsIRTConfiguration::GetNumberOfReactions() {
	return (G4int)fReactions.size();
}


std::pair<G4String, G4String> TsIRTConfiguration::GetReactants(G4int ReactionIndex) {
	std::pair<G4String, G4String> Reactants;
	
	Reactants.first  = fMoleculesName[fReactions[ReactionIndex].reactorA];
	Reactants.second = fMoleculesName[fReactions[ReactionIndex].reactorB];
	
	return Reactants;
}


std::vector<G4String> TsIRTConfiguration::GetProducts(G4int ReactionIndex) {
	std::vector<G4String> Products;
	
	for ( size_t i = 0; i < (fReactions[ReactionIndex].products).size(); i++) {
		Products.push_back(fMoleculesName[(fReactions[ReactionIndex].products)[i]]);
	}
	
	return Products;
}


G4double TsIRTConfiguration::brents_fun(TsMolecule molA, TsMolecule molB, G4int indexOfReaction, G4double offset) {
	G4double lower = fLowerTime;
	G4double upper = fUpperTime;
	G4double tol = 0.001*ps;
	unsigned int max_iter = 1000;
	
	G4double a = lower;
	G4double b = upper;
	G4double fa = CalculateProbabilityPartiallyDiffusionControlled(molA, molB, indexOfReaction, a)+offset;
	G4double fb = CalculateProbabilityPartiallyDiffusionControlled(molA, molB, indexOfReaction, b)+offset;
	G4double fs = 0;
	
	if (!(fa * fb < 0)) {
		return -11;
	}
	
	if (std::abs(fa) < std::abs(b)) {
		std::swap(a,b);
		std::swap(fa,fb);
	}
	
	G4double c = a;
	G4double fc = fa;
	G4bool mflag = true;
	G4double ss = 0;
	G4double d = 0;
	
	for (unsigned int iter = 1; iter < max_iter; ++iter) {
		// stop if converged on root or error is less than tolerance
		if (std::abs(b-a) < tol) {
			return ss;
		}
		
		if (fa != fc && fb != fc) { // use inverse quadratic interopolation
			ss =      ( a * fb * fc / ((fa - fb) * (fa - fc)) )
			+ ( b * fa * fc / ((fb - fa) * (fb - fc)) )
			+ ( c * fa * fb / ((fc - fa) * (fc - fb)) );
		}
		else{ // secant method
			ss = b - fb * (b - a) / (fb - fa);
		}
		
		// checks to see whether we can use the faster converging quadratic && secant methods or if we need to use bisection
		if (    ( (ss < (3 * a + b) * 0.25) || (ss > b) ) ||
			( mflag && (std::abs(ss-b) >= (std::abs(b-c) * 0.5)) ) ||
			( !mflag && (std::abs(ss-b) >= (std::abs(c-d) * 0.5)) ) ||
			( mflag && (std::abs(b-c) < tol) ) ||
			( !mflag && (std::abs(c-d) < tol))    ) {
			// bisection method
			ss = (a+b)*0.5;
			mflag = true;
		}
		else {
			mflag = false;
		}
		
		fs = CalculateProbabilityPartiallyDiffusionControlled(molA, molB, indexOfReaction, ss)+offset;
		d = c;
		c = b;
		fc = fb;
		
		if ( fa * fs < 0) {
			b = ss;
			fb = fs;
		} else {
			a = ss;
			fa = fs;
		}
		
		if (std::abs(fa) < std::abs(fb)) {
			std::swap(a,b);
			std::swap(fa,fb);
		}
	}
	return -1.0*ps;
}


G4double TsIRTConfiguration::CalculateProbabilityOfScavenger(TsMolecule, G4int indexOfReaction, G4double t) {
	// Vars: Reff, k, [B]=Cs
	G4double Reff = (fReactions[indexOfReaction]).effectiveReactionRadius;
	G4double CsB = (fReactions[indexOfReaction]).probabilityOfReaction/s;
	G4double nm3persToPerMPers = 6.022140857e-1;
	G4double k = (fReactions[indexOfReaction]).kobs/nm3persToPerMPers * (nm*nm*nm/s);;
	
	G4double Wscav = 1.0 - std::exp(-CsB * (t + 4.0 * Reff * std::sqrt(t * Reff / k )));
	return Wscav;
}


G4double TsIRTConfiguration::brents_fun_scav(TsMolecule molA, G4int indexOfReaction, G4double offset) {
	G4double lower = fLowerTime;
	G4double upper = fUpperTime;
	G4double tol = 0.001*ps;
	unsigned int max_iter = 100000;
	
	G4double a = lower;
	G4double b = upper;
	G4double fa = CalculateProbabilityOfScavenger(molA, indexOfReaction, a)+offset;
	G4double fb = CalculateProbabilityOfScavenger(molA, indexOfReaction, b)+offset;
	G4double fs = 0;
	
	if (!(fa * fb < 0)) {
		return -11;
	}
	
	if (std::abs(fa) < std::abs(b)) {
		std::swap(a,b);
		std::swap(fa,fb);
	}
	
	G4double c = a;
	G4double fc = fa;
	G4bool mflag = true;
	G4double ss = 0;
	G4double d = 0;
	
	for (unsigned int iter = 1; iter < max_iter; ++iter) {
		// stop if converged on root or error is less than tolerance
		if (std::abs(b-a) < tol) {
			return ss;
		}
		
		if (fa != fc && fb != fc) { // use inverse quadratic interopolation
			ss =      ( a * fb * fc / ((fa - fb) * (fa - fc)) )
			+ ( b * fa * fc / ((fb - fa) * (fb - fc)) )
			+ ( c * fa * fb / ((fc - fa) * (fc - fb)) );
		}
		else{ // secant method
			ss = b - fb * (b - a) / (fb - fa);
		}
		
		// checks to see whether we can use the faster converging quadratic && secant methods or if we need to use bisection
		if (    ( (ss < (3 * a + b) * 0.25) || (ss > b) ) ||
			( mflag && (std::abs(ss-b) >= (std::abs(b-c) * 0.5)) ) ||
			( !mflag && (std::abs(ss-b) >= (std::abs(c-d) * 0.5)) ) ||
			( mflag && (std::abs(b-c) < tol) ) ||
			( !mflag && (std::abs(c-d) < tol))    ) {
			// bisection method
			ss = (a+b)*0.5;
			mflag = true;
		}
		else {
			mflag = false;
		}
		
		fs = CalculateProbabilityOfScavenger(molA, indexOfReaction, ss)+offset;
		d = c;
		c = b;
		fc = fb;
		
		if ( fa * fs < 0) {
			b = ss;
			fb = fs;
		} else {
			a = ss;
			fa = fs;
		}
		
		if (std::abs(fa) < std::abs(fb)) {
			std::swap(a,b);
			std::swap(fa,fb);
		}
	}
	return -1.0*ps;
}


G4bool TsIRTConfiguration::MakeReaction(std::vector<TsMolecule> &initialSpecies,
										std::map<G4int, std::map<G4int, std::map<G4int,
										std::vector<G4int>>>> &spaceBinned,
										G4int NX, G4int NY, G4int NZ, G4double XMin, G4double XMax, G4double YMin,
										G4double YMax, G4double ZMin, G4double ZMax,
										std::map<G4int, std::map<G4int, G4int>> &theGvalue,
										std::vector<G4double> timeSteps,
										G4int iM, G4int indexOfReaction, G4double irt,
										std::vector<G4bool> &used) {
	
	G4ThreeVector positions = initialSpecies[iM].position;
	std::vector<G4int> products = (GetReaction(indexOfReaction)).products;
	G4int tBin = fUtils->FindBin(irt, timeSteps);
	if ( tBin < 0 ) return false;
	
	for ( size_t u = 0; u < products.size(); u++ ) {
		TsMolecule aProd;
		aProd.id = products[u];
		aProd.position = positions;
		aProd.time = irt;
		aProd.trackID = 0;
		aProd.reacted = false;
		aProd.isDNA = false;
		if ( products[u] == 1 || products[u] == 5)
			aProd.spin = G4UniformRand() > 0.5 ? 1 : 0;
		else
			aProd.spin = -1;
		
		G4int i = fUtils->FindBin(NX, XMin, XMax, positions.x());
		G4int j = fUtils->FindBin(NY, YMin, YMax, positions.y());
		G4int k = fUtils->FindBin(NZ, ZMin, ZMax, positions.z());
		
		spaceBinned[i][j][k].push_back( int(initialSpecies.size()) );
		
		initialSpecies.push_back(aProd);
		
		for ( int ti = tBin; ti < (int)timeSteps.size(); ti++ ) {
			theGvalue[aProd.id][ti]++;
		}
		used.push_back(false);
	}
	for ( int ti = tBin; ti < (int)timeSteps.size(); ti++ )
		theGvalue[initialSpecies[iM].id][ti]--;
	
	initialSpecies[iM].reacted = true;
	return true;
}


G4bool TsIRTConfiguration::MakeReaction(std::vector<TsMolecule> &initialSpecies,
										std::map<G4int, std::map<G4int, std::map<G4int,
										std::vector<G4int>>>> &spaceBinned,
										G4int NX, G4int NY, G4int NZ, G4double XMin, G4double XMax, G4double YMin,
										G4double YMax, G4double ZMin, G4double ZMax,
										std::map<G4int, std::map<G4int, G4int>> &theGvalue,
										std::vector<G4double> timeSteps,
										G4int iM, G4int jM, G4int indexOfReaction, G4double irt,
										G4double probabilityOfReaction, std::vector<G4bool> &used) {
	
	if ( G4UniformRand() < probabilityOfReaction ) {
		std::vector<G4ThreeVector> positions = GetPositionOfProducts(initialSpecies[iM], initialSpecies[jM], indexOfReaction);
		std::vector<G4int> products = (GetReaction(indexOfReaction)).products;
		G4int tBin = fUtils->FindBin(irt, timeSteps);
		if ( tBin < 0 ) return false;
		
		for ( size_t u = 0; u < products.size(); u++ ) {
			TsMolecule aProd;
			aProd.id = products[u];
			aProd.position = positions[u];
			aProd.time = irt;
			aProd.trackID = 0;
			aProd.reacted = false;
			aProd.isDNA = false;
			if ( products[u] == 1 || products[u] == 5)
				aProd.spin = G4UniformRand() > 0.5 ? 1 : 0;
			else
				aProd.spin = -1;
			
			G4int i = fUtils->FindBin(NX, XMin, XMax, positions[u].x());
			G4int j = fUtils->FindBin(NY, YMin, YMax, positions[u].y());
			G4int k = fUtils->FindBin(NZ, ZMin, ZMax, positions[u].z());
			
			spaceBinned[i][j][k].push_back( int(initialSpecies.size()) );
			
			initialSpecies.push_back(aProd);
			
			for ( int ti = tBin; ti < (int)timeSteps.size(); ti++ ) {
				theGvalue[aProd.id][ti]++;
			}
			used.push_back(false);
		}
		for ( int ti = tBin; ti < (int)timeSteps.size(); ti++ ) {
			theGvalue[initialSpecies[iM].id][ti]--;
			theGvalue[initialSpecies[jM].id][ti]--;
		}
		initialSpecies[iM].reacted = true;
		initialSpecies[jM].reacted = true;
		return true;
	} else {
		return false;
	}
}


G4bool TsIRTConfiguration::Inside(G4ThreeVector p ) {
	G4double delta = 0.5 * 1e-3 * nm;
	G4double fDx = 0.5*um;
	G4double fDy = 0.5*um;
	G4double fDz = 0.5*um;
	G4double dist = std::max(std::max(
									  std::abs(p.x())-fDx,
									  std::abs(p.y())-fDy),
							 std::abs(p.z())-fDz);
	if (dist > delta) return false;
	return true;
}


G4bool TsIRTConfiguration::MakeReaction(std::vector<TsMolecule> &initialSpecies,
										std::map<G4int, std::map<G4int, std::map<G4int,
										std::vector<G4int>>>> &spaceBinned,
										G4int NX, G4int NY, G4int NZ, G4double XMin, G4double XMax, G4double YMin,
										G4double YMax, G4double ZMin, G4double ZMax,
										std::map<G4int, std::map<G4int, G4int>> &theGvalue,
										std::map<G4int, std::map<G4int, G4int>> &theGvalueInVolume,
										std::vector<G4double> timeSteps,
										G4int iM, G4int jM, G4int indexOfReaction, G4double irt,
										G4double probabilityOfReaction, std::vector<G4bool> &used) {
	
	if ( G4UniformRand() < probabilityOfReaction ) {
		std::vector<G4ThreeVector> positions = GetPositionOfProducts(initialSpecies[iM], initialSpecies[jM], indexOfReaction);
		std::vector<G4int> products = (GetReaction(indexOfReaction)).products;
		G4int tBin = fUtils->FindBin(irt, timeSteps);
		if ( tBin < 0 ) return false;
		
		G4bool inVolume = false;
		if ( Inside(initialSpecies[iM].position) || Inside(initialSpecies[jM].position))  // at least one specie is at scoring region
			inVolume = true;
		
		for ( size_t u = 0; u < products.size(); u++ ) {
			TsMolecule aProd;
			aProd.id = products[u];
			aProd.position = positions[u];
			aProd.time = irt;
			aProd.trackID = 0;
			aProd.reacted = false;
			aProd.isDNA = false;
			if ( products[u] == 1 || products[u] == 5)
				aProd.spin = G4UniformRand() > 0.5 ? 1 : 0;
			else
				aProd.spin = -1;
			
			G4int i = fUtils->FindBin(NX, XMin, XMax, positions[u].x());
			G4int j = fUtils->FindBin(NY, YMin, YMax, positions[u].y());
			G4int k = fUtils->FindBin(NZ, ZMin, ZMax, positions[u].z());
			
			spaceBinned[i][j][k].push_back( int(initialSpecies.size()) );
			
			initialSpecies.push_back(aProd);
			
			for ( int ti = tBin; ti < (int)timeSteps.size(); ti++ ) {
				theGvalue[aProd.id][ti]++;
				if ( inVolume )
					theGvalueInVolume[aProd.id][ti]++;
			}
			used.push_back(false);
		}
		for ( int ti = tBin; ti < (int)timeSteps.size(); ti++ ) {
			theGvalue[initialSpecies[iM].id][ti]--;
			theGvalue[initialSpecies[jM].id][ti]--;
			if ( inVolume ) {
				theGvalueInVolume[initialSpecies[iM].id][ti]--;
				theGvalueInVolume[initialSpecies[jM].id][ti]--;
			}
		}
		initialSpecies[iM].reacted = true;
		initialSpecies[jM].reacted = true;
		return true;
	} else {
		return false;
	}
}



void TsIRTConfiguration::ScoreGvalue(std::vector<TsMolecule> &initialSpecies,
									 std::map<G4int, std::map<G4int, G4int>> &theGvalueInVolume,
									 std::vector<G4double> timeSteps,
									 G4int iM, G4int jM, G4int indexOfReaction, G4double irt) {
	
	std::vector<G4int> products = (GetReaction(indexOfReaction)).products;
	G4int tBin = fUtils->FindBin(irt, timeSteps);
	if ( tBin < 0 ) return;
	
	for ( size_t u = 0; u < products.size(); u++ ) {
		for ( int ti = tBin; ti < (int)timeSteps.size(); ti++ ) {
			theGvalueInVolume[products[u]][ti]++;
		}
	}
	
	for ( int ti = tBin; ti < (int)timeSteps.size(); ti++ ) {
		theGvalueInVolume[initialSpecies[iM].id][ti]--;
		theGvalueInVolume[initialSpecies[jM].id][ti]--;
	}
	
	return;
}


std::vector<G4double> TsIRTConfiguration::GetH2SO4ComponentsConcentrationPH(G4double pH) {
	G4double Ka1 = pow(10,3);
	G4double Ka2 = pow(10,-1.987);
	G4double _Kw = 1E-14;
	
	G4double _H_pos = pow(10,-pH);
	G4double _SO_4    = (_H_pos - (_Kw / _H_pos)) / ((_H_pos / Ka2) + 2);
	G4double _HSO_4   = _SO_4 * _H_pos / Ka2;
	G4double _OH_me   = _Kw / _H_pos;
	G4double _H_2SO_4 = _HSO_4 * _H_pos / Ka1;
	G4double H_Poly = 0;
	
	std::vector<G4double> Results = {_H_pos, _HSO_4, _SO_4, _OH_me, _H_2SO_4, H_Poly};
	
	return Results;
}


std::vector<G4double> TsIRTConfiguration::GetH2SO4ComponentsConcentrationP(G4double Concentration) {
	G4double Ka1 = pow(10,3);
	G4double Ka2 = pow(10,-1.987);
	G4double _C = Concentration;
	G4double _Kw = 1E-14;
	
	G4double _p4 = 1;
	G4double _p3 = Ka1;
	G4double _p2 = (Ka1 * Ka2) - (_C * Ka1) - _Kw;
	G4double _p1 = - ((Ka1 * _Kw) + (2 * _C * Ka1 * Ka2));
	G4double _p0 = - (_Kw * Ka1 * Ka2);
	
	std::vector<G4double> P = {_p4, _p3, _p2, _p1, _p0};
	
	std::vector<G4double> Result = GetRoots(4, P);
	
	G4double _H_pos = -10;
	
	for (size_t i = 0; i < Result.size(); i++) {
		if ((Result[i] > 0 and Result[i] < 3.0*_C) and i % 2 != 1) {
			if (Result[i + 1] == 0) {
				_H_pos = Result[i];
			}
		}
	}
	
	if (_H_pos == -10)
		return {0, 0, 0, 0, 0, 0};
	
	G4double _SO_4    = (_H_pos - (_Kw / _H_pos)) / ((_H_pos / Ka2) + 2);
	G4double _HSO_4   = _SO_4 * _H_pos / Ka2;
	G4double _OH_me   = _Kw / _H_pos;
	G4double _H_2SO_4 = _HSO_4 * _H_pos / Ka1;
	
	G4double H_Poly = (P[0] * pow(_H_pos,4)) + (P[1] * pow(_H_pos,3)) +
	(P[2] * pow(_H_pos,2)) + (P[3] * _H_pos) + (P[4]);
	
	std::vector<G4double> Results = {_H_pos, _HSO_4, _SO_4, _OH_me, _H_2SO_4, H_Poly};
	
	return Results;
}


G4double TsIRTConfiguration::GetIonicStrength(std::vector<G4double> Components) {
	G4double I = .5 * ((Components[0]) + (Components[1]) + (Components[2] * 4) + (Components[3]));
	return I;
}


G4double TsIRTConfiguration::IonicRate(G4double IonicStrength, G4double Rate, G4int Charge1, G4int Charge2) {
	G4double a = 0.15;
	G4double b = a * Charge1 * Charge2;
	G4double scaleFactor = pow(10,((1.02 * Charge1 * Charge2 * ((sqrt(IonicStrength)) / (1 + sqrt(IonicStrength))))  -  (2 * b * IonicStrength)));
	G4double K = Rate * scaleFactor;
	return K;
}


G4double TsIRTConfiguration::IonicRate(G4double IonicStrength, TsMolecularReaction Reaction) {
	G4double a = 0.15;
	G4int Charge1 = fMoleculesDefinition[Reaction.reactorA].charge;
	G4int Charge2 = fMoleculesDefinition[Reaction.reactorB].charge;
	G4int BackGround = Reaction.reactionType;
	G4double Rate = 0.0;
	if ( BackGround == 6 ) {
		Rate = Reaction.scavengingCapacity;
		Rate /= 1E-7;
	} else {
		Rate = Reaction.kobs;
	}
	
	G4double b = a * Charge1 * Charge2;
	G4double scaleFactor = pow(10,((1.02 * Charge1 * Charge2 * ((sqrt(IonicStrength)) / (1 + sqrt(IonicStrength))))  -  (2 * b * IonicStrength)));
	G4double K = Rate * scaleFactor;
	return K;
}


void TsIRTConfiguration::roots(double *a,int n,double *wr,double *wi) {
	double sq,b2,c,disc;
	int i,numroots;
	
	i = n;
	numroots = 0;
	while (i > 1) {
		b2 = -0.5*a[i-2];
		c = a[i-1];
		disc = b2*b2-c;
		if (disc < 0.0) {                   // complex roots
			sq = sqrt(-disc);
			wr[i-2] = b2;
			wi[i-2] = sq;
			wr[i-1] = b2;
			wi[i-1] = -sq;
			numroots+=2;
		}
		else {                              // real roots
			sq = sqrt(disc);
			wr[i-2] = fabs(b2)+sq;
			if (b2 < 0.0) wr[i-2] = -wr[i-2];
			if (wr[i-2] == 0)
				wr[i-1] = 0;
			else {
				wr[i-1] = c/wr[i-2];
				numroots+=2;
			}
			wi[i-2] = 0.0;
			wi[i-1] = 0.0;
		}
		i -= 2;
	}
	if (i == 1) {
		wr[0] = -a[0];
		wi[0] = 0.0;
		numroots++;
	}
}


void TsIRTConfiguration::deflate(double *a,int n,double *b,double *quad,double *err) {
	double r,z;
	int i;
	
	r = quad[1];
	z = quad[0];
	
	b[1] = a[1] - r;
	
	for (i=2;i<=n;i++){
		b[i] = a[i] - r * b[i-1] - z * b[i-2];
	}
	*err = fabs(b[n])+fabs(b[n-1]);
}


void TsIRTConfiguration::find_quad(double *a,int n,double *b,double *quad,double *err, int *iter) {
	double maxiter = 500;
	double *c,dn,dr,ds,drn,dsn,eps,r,z;
	int i;
	
	c = new double [n+1];
	c[0] = 1.0;
	r = quad[1];
	z = quad[0];
	eps = 1e-15;
	*iter = 1;
	
	do {
		if (*iter > maxiter) break;
		if (((*iter) % 200) == 0) {
			eps *= 10.0;
		}
		b[1] = a[1] - r;
		c[1] = b[1] - r;
		
		for (i=2;i<=n;i++){
			b[i] = a[i] - r * b[i-1] - z * b[i-2];
			c[i] = b[i] - r * c[i-1] - z * c[i-2];
		}
		dn=c[n-1] * c[n-3] - c[n-2] * c[n-2];
		drn=b[n] * c[n-3] - b[n-1] * c[n-2];
		dsn=b[n-1] * c[n-1] - b[n] * c[n-2];
		
		if (fabs(dn) < 1e-10) {
			if (dn < 0.0) dn = -1e-8;
			else dn = 1e-8;
		}
		dr = drn / dn;
		ds = dsn / dn;
		r += dr;
		z += ds;
		(*iter)++;
	} while ((fabs(dr)+fabs(ds)) > eps);
	quad[0] = z;
	quad[1] = r;
	*err = fabs(ds)+fabs(dr);
	delete [] c;
}


void TsIRTConfiguration::diff_poly(double *a,int n,double *b) {
	double coef;
	int i;
	
	coef = (double)n;
	b[0] = 1.0;
	for (i=1;i<n;i++) {
		b[i] = a[i]*((double)(n-i))/coef;
	}
}


void TsIRTConfiguration::recurse(double *a,int n,double *b,int p,double *quad,
								 double *err,int *iter) {
	double *c,*x,rs[2],tst;
	
	if (fabs(b[p]) < 1e-16) p--;    // this bypasses roots at zero
	if (p == 2) {
		quad[0] = b[2];
		quad[1] = b[1];
		*err = 0;
		*iter = 0;
		return;
	}
	c = new double [p+1];
	x = new double [n+1];
	c[0] = x[0] = 1.0;
	rs[0] = quad[0];
	rs[1] = quad[1];
	*iter = 0;
	find_quad(b,p,c,rs,err,iter);
	tst = fabs(rs[0]-quad[0])+fabs(rs[1]-quad[1]);
	if (*err < 1e-12) {
		quad[0] = rs[0];
		quad[1] = rs[1];
	}
	
	if (((*iter > 5) && (tst < 1e-4)) || ((*iter > 20) && (tst < 1e-1))) {
		diff_poly(b,p,c);
		recurse(a,n,c,p-1,rs,err,iter);
		quad[0] = rs[0];
		quad[1] = rs[1];
	}
	delete [] x;
	delete [] c;
}


void TsIRTConfiguration::get_quads(double *a,int n,double *quad,double *x) {
	double maxiter = 500;
	double *b,*z,err,tmp;
	int iter,i,p;
	
	if ((tmp = a[0]) != 1.0) {
		a[0] = 1.0;
		for (i=1;i<=n;i++) {
			a[i] /= tmp;
		}
	}
	if (n == 2) {
		x[0] = a[1];
		x[1] = a[2];
		return;
	}
	else if (n == 1) {
		x[0] = a[1];
		return;
	}
	p = n;
	b = new double [n+1];
	z = new double [n+1];
	b[0] = 1.0;
	for (i=0;i<=n;i++) {
		z[i] = a[i];
		x[i] = 0.0;
	}
	do {
		if (n > p) {
			quad[0] = 3.14159e-1;
			quad[1] = 2.78127e-1;
		}
		do {
			for (i=0;i<5;i++) {
				find_quad(z,p,b,quad,&err,&iter);
				if ((err > 1e-7) || (iter > maxiter)) {
					diff_poly(z,p,b);
					recurse(z,p,b,p-1,quad,&err,&iter);
				}
				deflate(z,p,b,quad,&err);
				if (err < 0.001) break;
				//quad[0] = random(8) - 4.0;
				//quad[1] = random(8) - 4.0;
				quad[0] = (rand() % 8) - 4;
				quad[i] = (rand() % 8) - 4;
			}
			if (err > 0.01) {
				std::cout << "Error! Convergence failure in quadratic x^2 + r*x + s." << std::endl;
				std::cout << "Enter new trial value for 'r': ";
				std::cin >> quad[1];
				std::cout << "Enter new trial value for 's' ( 0 to exit): ";
				std::cin >> quad[0];
				if (quad[0] == 0) {
					G4cerr << "TOPAS is exiting due to fatal error en IRT Chemistry setup!" << G4endl;
					G4cerr << "--- No valid solution to the pH scaling system found!" << G4endl;
					fPm->AbortSession(1);
				}
			}
		} while (err > 0.01);
		x[p-2] = quad[1];
		x[p-1] = quad[0];
		p -= 2;
		for (i=0;i<=p;i++) {
			z[i] = b[i];
		}
	} while (p > 2);
	if (p == 2) {
		x[0] = b[1];
		x[1] = b[2];
	}
	else x[0] = b[1];
	delete [] z;
	delete [] b;
}

std::vector<double> TsIRTConfiguration::GetRoots(int Order, std::vector<double> Exponents) {
	
	double a[21],x[21],wr[21],wi[21],quad[2];
	int n,i;
	
	n = Order;
	if ((n < 1) || (n > 20)) {
		std::cout << "Error! Invalid order: n = " << n << std::endl;
		return {0};
	}
	
	for (i = 0; i <= n; i++) {
		a[i] = Exponents[i];
		if (a[0] == 0) {
			std::cout << "Error! Highest coefficient cannot be 0." << std::endl;
			return {0};
		}
	}
	if (a[n] == 0) {
		std::cout << "Error! Lowest coefficient (constant term) cannot be 0." << std::endl;
		return {0};
	}
	
	quad[0] = 2.71828e-1;
	quad[1] = 3.14159e-1;
	
	get_quads(a,n,quad,x);
	roots(x,n,wr,wi);
	std::vector<double> Results;
	
	for (i=0;i<n;i++) {
		if ((wr[i] != 0.0) || (wi[i] != 0.0)){
			Results.push_back(wr[i]);
			Results.push_back(wi[i]);
		}
	}
	return Results;
}


G4double TsIRTConfiguration::k23_inverse(G4double t) {
	return 6.62 * pow(10,10) + (1.48 * pow(10,9) * t) + (1.28 * pow(10,7) * pow(t,2)) - (6.03 * pow(10,4) * pow(t,3)) + (1.28 * pow(10,2) * pow(t,4));
}


G4double TsIRTConfiguration::k26_inverse(G4double t) {
	G4double T = 273.15 + t;
	G4double p = 18.61 - (6.94 * pow(10,3) * pow(T,-1)) + (2.12 * pow(10,6) * pow(T,-2)) - (2.34 * pow(10,8) * pow(T,-3));
	return pow(10,p);
}


G4double TsIRTConfiguration::k27_inverse(G4double t) {
	G4double T = 273.15 + t;
	return pow(10,14) * 1.33 * std::exp(-38.38/(0.0083146*T));
}


G4double TsIRTConfiguration::k29_inverse(G4double t) {
	return (k_29(t) * K_water(t)) / (lH2Ol(t) * K_OH(t));
}


G4double TsIRTConfiguration::k30_inverse(G4double t) {
	return 3.41 * pow(10,10) + (2.75 * pow(10, 8) * t) + (1.24 * pow(10,7) * pow(t,2)) - (6.23 * pow(10,4) * pow(t,3)) + (1.31 * pow(10,2) * pow(t,4));
}


G4double TsIRTConfiguration::K_H(G4double t) {
	G4double p = 10.49 - (4.103 * pow(10,-2) * t) + (1.443 * pow(10,-4) * pow(t,2)) - (2.325 * pow(10,-7) * pow(t,3)) + (2.065 * pow(10,-10) * pow(t,4));
	return pow(10,-p);
}


G4double TsIRTConfiguration::K_OH(G4double t) {
	G4double p = 12.50 - (3.317 * pow(10,-2) * t) + (1.964 * pow(10,-4) * pow(t,2)) - (6.198 * pow(10,-7) * pow(t,3)) + (8.244 * pow(10,-10) * pow(t,4));
	return pow(10,-p);
}


G4double TsIRTConfiguration::K_H2O2(G4double t) {
	G4double p = 12.50 - 3.317 * pow(10,-2) * t + 1.964 * pow(10,-4) * pow(t,2) - 6.198 * pow(10,-7) * pow(t, 3) + 8.244 * pow(10,-10) * pow(t,4);
	return pow(10,(-p));
}


G4double TsIRTConfiguration::k_26(G4double t) {
	return K_H(t)*k26_inverse(t);
}


G4double TsIRTConfiguration::lH2Ol(G4double t) {
	return 55.55 * (1.007 - (4.449 * pow(10,-4) * t) - (3.155 * pow(10,-8) * pow(t,2)) - (6.088 * pow(10,-9) * pow(t,3)));
}


G4double TsIRTConfiguration::K_water(G4double t) {
	G4double p = 14.93 - (4.131 * pow(10,-2) * t) + (1.903 * pow(10,-4) * pow(t,2)) - (4.705 * pow(10,-7) * pow(t,3)) + (5.724 * pow(10,-10) * pow(t,4));
	return pow(10,-p);
}


G4double TsIRTConfiguration::k_28(G4double t) {
	return k30_inverse(t) * K_OH(t);
}


G4double TsIRTConfiguration::k_29(G4double t) {
	return 7.22 * pow(10,9) + (1.62 * pow(10,8) * t)  + (2.40 * pow(10,6) * pow(t,2)) - (7.81 * pow(10,3) * pow(t,3)) + (1.06 * 10 * pow(t,4));
}


G4double TsIRTConfiguration::k_30(G4double t) {
	return k30_inverse(t)*K_H2O2(t);
}

