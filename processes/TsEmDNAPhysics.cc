// Physics Module for TsEmDNAPhysics
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

#include "TsEmDNAPhysics.hh"
#include "TsParameterManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4DNAGenericIonsManager.hh"

#include "G4DNAElectronSolvation.hh"
#include "G4DNAElastic.hh"
#include "G4DNAChampionElasticModel.hh"
#include "G4DNAScreenedRutherfordElasticModel.hh"
#include "G4DNAELSEPAElasticModel.hh"

#include "G4DNAExcitation.hh"
#include "G4DNAAttachment.hh"
#include "G4DNAVibExcitation.hh"
#include "G4DNAIonisation.hh"
#include "G4DNAChargeDecrease.hh"
#include "G4DNAChargeIncrease.hh"

#include "G4hMultipleScattering.hh"
#include "G4LowEWentzelVIModel.hh"

#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4GenericIon.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"

#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4PhotoElectricEffect.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4ComptonScattering.hh"
#include "G4LivermoreComptonModel.hh"
#include "G4GammaConversion.hh"
#include "G4LivermoreGammaConversionModel.hh"
#include "G4RayleighScattering.hh"
#include "G4LivermoreRayleighModel.hh"
#include "G4EmParameters.hh"

#include "G4LossTableManager.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4PhysicsListHelper.hh"
#include "G4BuilderType.hh"

#include "G4PhysicsConstructorFactory.hh"

#include "G4DNAOneStepThermalizationModel.hh"
#include "G4DNAUeharaScreenedRutherfordElasticModel.hh"
#include "G4DNACPA100ElasticModel.hh"
#include "G4DNACPA100ExcitationModel.hh"
#include "G4DNACPA100IonisationModel.hh"
#include "G4DNAEmfietzoglouIonisationModel.hh"
#include "G4DNAEmfietzoglouExcitationModel.hh"
#include "G4DNABornIonisationModel.hh"
#include "G4DNABornExcitationModel.hh"
#include "G4DNAPTBElasticModel.hh"
#include "G4DNAPTBExcitationModel.hh"
#include "G4DNAPTBIonisationModel.hh"
#include "G4DNAModelInterface.hh"
#include "G4DNAVacuumModel.hh"

#include "G4DNARuddIonisationModel.hh"
#include "TsDNARuddIonisationExtendedModel.hh"
#include "TsSplitProcessG4DNA.hh"
#include "TsDNARuddIonisationExtendedRITRACKSModel.hh"
#include "G4ProcessManager.hh"

G4_DECLARE_PHYSCONSTR_FACTORY(TsEmDNAPhysics);

TsEmDNAPhysics::TsEmDNAPhysics(G4int ver, const G4String&)
: G4VPhysicsConstructor("TsEmDNAPhysics"), verbose(ver)
{
    G4EmParameters* param = G4EmParameters::Instance();
    param->SetDefaults();
    param->SetFluo(true);
    param->SetAuger(true);
    param->SetAugerCascade(true);
    param->SetDeexcitationIgnoreCut(true);
    param->ActivateDNA();
    SetPhysicsType(bElectromagnetic);
}


TsEmDNAPhysics::TsEmDNAPhysics(TsParameterManager* pM)
: G4VPhysicsConstructor("TsEmDNAPhysics"), fPm(pM), verbose(2)
{
    G4EmParameters* param = G4EmParameters::Instance();
    param->SetDefaults();
    param->SetFluo(true);
    param->SetAuger(true);
    param->SetAugerCascade(true);
    param->SetDeexcitationIgnoreCut(true);
    param->ActivateDNA();
    SetPhysicsType(bElectromagnetic);

    fName = "Default";
    if ( fPm->ParameterExists("Ph/ListName") )
        fName = fPm->GetStringParameter("Ph/ListName");
    G4cout << fName << G4endl;
}


TsEmDNAPhysics::~TsEmDNAPhysics()
{}


void TsEmDNAPhysics::ConstructParticle()
{
    G4Gamma::Gamma();
    G4Electron::Electron();
    G4Positron::Positron();
    G4Proton::Proton();
    G4GenericIon::GenericIonDefinition();
    
    G4DNAGenericIonsManager * genericIonsManager;
    genericIonsManager=G4DNAGenericIonsManager::Instance();
    genericIonsManager->GetIon("alpha++");
    genericIonsManager->GetIon("alpha+");
    genericIonsManager->GetIon("helium");
    genericIonsManager->GetIon("hydrogen");
}


void TsEmDNAPhysics::ConstructProcess()
{
    if(verbose > 1) {
        G4cout << "### " << GetPhysicsName() << " Construct Processes " << G4endl;
    }
    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
    
    auto myParticleIterator=GetParticleIterator();
    myParticleIterator->reset();
    while( (*myParticleIterator)() )
    {
        G4ParticleDefinition* particle = myParticleIterator->value();
        G4String particleName = particle->GetParticleName();
        
        if (particleName == "e-") {
            
            G4double solvationHighLimit = -1.0;
            
            // Setup elastic scattering model
            G4String eScatteringModel = "champion";
            if ( fPm->ParameterExists("Ph/"+fName+"/Electron/SetElasticScatteringModel"))
                eScatteringModel = fPm->GetStringParameter("Ph/"+fName+"/Electron/SetElasticScatteringModel");
            eScatteringModel.toLower();
            
            if ( eScatteringModel == "champion" ) {
                G4DNAElastic* theDNAElasticProcess = new G4DNAElastic("e-_G4DNAElastic");
                theDNAElasticProcess->SetEmModel(new G4DNAChampionElasticModel());
                ph->RegisterProcess(theDNAElasticProcess, particle);
                solvationHighLimit = 7.4 * eV;
                
			} else if ( eScatteringModel == "elsepa" ) {
				G4DNAElastic* theDNAElasticProcess = new G4DNAElastic("e-_G4DNAElastic");
				theDNAElasticProcess->SetEmModel(new G4DNAELSEPAElasticModel());
				ph->RegisterProcess(theDNAElasticProcess, particle);
				solvationHighLimit = 10.0 * eV;
				
			} else if ( eScatteringModel == "screenedrutherford" ) {
                G4DNAElastic* theDNAElasticProcess = new G4DNAElastic("e-_G4DNAElastic");
                theDNAElasticProcess->SetEmModel(new G4DNAScreenedRutherfordElasticModel());
                ph->RegisterProcess(theDNAElasticProcess, particle);
                solvationHighLimit = 9.0 * eV;
                
            } else if ( eScatteringModel == "ueharascreenedrutherford" ) {
                G4DNAElastic* theDNAElasticProcess = new G4DNAElastic("e-_G4DNAElastic");
                theDNAElasticProcess->SetEmModel(new G4DNAUeharaScreenedRutherfordElasticModel());
                ph->RegisterProcess(theDNAElasticProcess, particle);
                solvationHighLimit = 9.0 * eV;
                
            } else if ( eScatteringModel == "cpa100" ) {
                G4DNAElastic* theDNAElasticProcess = new G4DNAElastic("e-_G4DNAElastic");
                G4VEmModel* model1 = new G4DNACPA100ElasticModel();
                G4VEmModel* model2 = new G4DNAChampionElasticModel();
                
                model1->SetActivationLowEnergyLimit(11*eV);
                model1->SetActivationHighEnergyLimit(255.955*keV);
                theDNAElasticProcess->SetEmModel(model1);
                
                model2->SetActivationLowEnergyLimit(255.955*keV);
                model2->SetActivationHighEnergyLimit(1*MeV);
                theDNAElasticProcess->AddEmModel(2,model2);
                ph->RegisterProcess(theDNAElasticProcess, particle);
                
                solvationHighLimit = 11 * eV;
                
            } else if ( eScatteringModel.contains("dna") ) {
                if ( eScatteringModel == "dnascreenedrutherford") {
                    G4DNAScreenedRutherfordElasticModel* e_modelDNARutherfordElastic = new G4DNAScreenedRutherfordElasticModel();
                    G4DNAPTBElasticModel* modelDNAPTBElastic = new G4DNAPTBElasticModel("THF/TMP/PY", particle);
                    
                    G4DNAModelInterface* e_elasticInteraction = new G4DNAModelInterface("e-_elastic_interaction");
                    e_elasticInteraction->RegisterModel(e_modelDNARutherfordElastic, particle);
                    e_elasticInteraction->RegisterModel(modelDNAPTBElastic);
                    e_elasticInteraction->RegisterModel(new G4DNAVacuumModel());
                    
                    G4DNAElastic* e_DNAElasticProcess = new G4DNAElastic("e-_G4DNAElastic");
                    e_DNAElasticProcess->SetEmModel(e_elasticInteraction);
                    
                    ph->RegisterProcess(e_DNAElasticProcess, particle);
                    
                } else if ( eScatteringModel == "dnacpa100" ) {
                    G4DNACPA100ElasticModel* e_modelDNACPA100Model = new G4DNACPA100ElasticModel();
                    G4DNAPTBElasticModel* modelDNAPTBElastic = new G4DNAPTBElasticModel("THF/TMP/PY", particle);
                    
                    G4DNAModelInterface* e_elasticInteraction = new G4DNAModelInterface("e-_elastic_interaction");
                    e_elasticInteraction->RegisterModel(e_modelDNACPA100Model, particle);
                    e_elasticInteraction->RegisterModel(modelDNAPTBElastic);
                    e_elasticInteraction->RegisterModel(new G4DNAVacuumModel());
                    
                    G4DNAElastic* e_DNAElasticProcess = new G4DNAElastic("e-_G4DNAElastic");
                    e_DNAElasticProcess->SetEmModel(e_elasticInteraction);
                    
                    ph->RegisterProcess(e_DNAElasticProcess, particle);
                    
                } else if ( eScatteringModel == "dnaueharascreenedrutherford" ) {
                    G4DNAUeharaScreenedRutherfordElasticModel* e_modelDNAUeharaScreenedRutherfordElastic = new G4DNAUeharaScreenedRutherfordElasticModel();
                    
                    G4DNAPTBElasticModel* modelDNAPTBElastic = new G4DNAPTBElasticModel("THF/TMP/PY", particle);
                    
                    G4DNAModelInterface* e_elasticInteraction = new G4DNAModelInterface("e-_elastic_interaction");
                    e_elasticInteraction->RegisterModel(e_modelDNAUeharaScreenedRutherfordElastic, particle);
                    e_elasticInteraction->RegisterModel(modelDNAPTBElastic);
                    e_elasticInteraction->RegisterModel(new G4DNAVacuumModel());
                    
                    G4DNAElastic* e_DNAElasticProcess = new G4DNAElastic("e-_G4DNAElastic");
                    e_DNAElasticProcess->SetEmModel(e_elasticInteraction);
                    
                    ph->RegisterProcess(e_DNAElasticProcess, particle);
                } else {
                    std::cout << "Error! Elastic model " << eScatteringModel << " not found " << std::endl;
                }
                
            } else { // WentzelVI
                G4eMultipleScattering* msc = new G4eMultipleScattering();
                msc->SetEmModel(new G4LowEWentzelVIModel());
                ph->RegisterProcess(msc, particle);
            }
            
            // Setup excitation model
            G4String eExcitationModel = "born";
            if ( fPm->ParameterExists("Ph/"+fName+"/Electron/SetExcitationModel") )
                eExcitationModel = fPm->GetStringParameter("Ph/"+fName+"/Electron/SetExcitationModel");
            eExcitationModel.toLower();
            
            G4String eIonisationModel = "born";
            if ( fPm->ParameterExists("Ph/"+fName+"/Electron/SetIonisationModel") )
                eIonisationModel = fPm->GetStringParameter("Ph/"+fName+"/Electron/SetIonisationModel");
            eIonisationModel.toLower();
            
            if ( eExcitationModel == "emfietzoglou" ) {
                // *** Excitation ***
                G4DNAExcitation* theDNAExcitationProcess = new G4DNAExcitation("e-_G4DNAExcitation");
                G4VEmModel* model1 = new G4DNAEmfietzoglouExcitationModel();
                G4VEmModel* model2 = new G4DNABornExcitationModel();
                
                model1->SetActivationLowEnergyLimit(8*eV);
                model1->SetActivationHighEnergyLimit(10*keV);
                model1->SetLowEnergyLimit(8*eV);
                model1->SetHighEnergyLimit(10*keV);
                
                model2->SetActivationLowEnergyLimit(10*keV);
                model2->SetActivationHighEnergyLimit(1*MeV);
                model2->SetLowEnergyLimit(10*keV);
                model2->SetHighEnergyLimit(1*MeV);
                
                theDNAExcitationProcess->AddEmModel(1,model1);
                theDNAExcitationProcess->AddEmModel(2,model2);
                
                ph->RegisterProcess(theDNAExcitationProcess, particle);
                
            } else if ( eExcitationModel == "cpa100" ) {
                // *** Excitation ***
                G4DNAExcitation* theDNAExcitationProcess = new G4DNAExcitation("e-_G4DNAExcitation");
                G4DNACPA100ExcitationModel* model1 = new G4DNACPA100ExcitationModel();
                G4DNABornExcitationModel* model2 = new G4DNABornExcitationModel();
                
                model1->SetActivationLowEnergyLimit(11*eV);
                model1->SetActivationHighEnergyLimit(255.955*keV);
                model1->SetLowEnergyLimit(11*eV);
                model1->SetHighEnergyLimit(255.955*keV);
                
                model2->SetLowEnergyLimit(255.955*keV);
                model2->SetHighEnergyLimit(1*MeV);
                model2->SetActivationLowEnergyLimit(255.955*keV);
                model2->SetActivationHighEnergyLimit(1*MeV);
                
                theDNAExcitationProcess->AddEmModel(1,model1);
                theDNAExcitationProcess->AddEmModel(2,model2);
                
                ph->RegisterProcess(theDNAExcitationProcess, particle);
                
            } else if ( eExcitationModel.contains("dna") ) {
                if ( eExcitationModel == "dnaemfietzoglou" ) {
                    G4DNAEmfietzoglouExcitationModel* e_modelDNAEmfietzoglouExcitation = new G4DNAEmfietzoglouExcitationModel();
                    G4DNAPTBExcitationModel* modelDNAPTBExcitation = new G4DNAPTBExcitationModel("THF/TMP/PY",particle);
                    G4DNAModelInterface* e_excitationInteraction = new G4DNAModelInterface("e-_excitation_interaction");
                    
                    e_excitationInteraction->RegisterModel(e_modelDNAEmfietzoglouExcitation,particle);
                    e_excitationInteraction->RegisterModel(modelDNAPTBExcitation);
                    e_excitationInteraction->RegisterModel(new G4DNAVacuumModel());
                    
                    G4DNAExcitation* e_DNAExcitationProcess = new G4DNAExcitation("e-_G4DNAExcitation");
                    e_DNAExcitationProcess->SetEmModel(e_excitationInteraction);
                    
                    ph->RegisterProcess(e_DNAExcitationProcess, particle);
                } else if ( eExcitationModel == "dnaborn" ) {
                    G4DNABornExcitationModel* e_modelDNABornExcitation = new G4DNABornExcitationModel();
                    G4DNAPTBExcitationModel* modelDNAPTBExcitation = new G4DNAPTBExcitationModel("THF/TMP/PY",particle);
                    G4DNAModelInterface* e_excitationInteraction = new G4DNAModelInterface("e-_excitation_interaction");
                    
                    e_excitationInteraction->RegisterModel(e_modelDNABornExcitation,particle);
                    e_excitationInteraction->RegisterModel(modelDNAPTBExcitation);
                    e_excitationInteraction->RegisterModel(new G4DNAVacuumModel());
                    
                    G4DNAExcitation* e_DNAExcitationProcess = new G4DNAExcitation("e-_G4DNAExcitation");
                    e_DNAExcitationProcess->SetEmModel(e_excitationInteraction);
                    
                    ph->RegisterProcess(e_DNAExcitationProcess, particle);
                } else {
                    std::cout << "Error! Excitation model " << eExcitationModel << " not found " << std::endl;
                }
            } else { // Born
                ph->RegisterProcess(new G4DNAExcitation("e-_G4DNAExcitation"), particle);
            }
            
            // Setup ionisation model
            G4bool useVarianceReduction = false;
            if ( fPm->ParameterExists("Vr/UseG4DNAVarianceReduction") ) {
                useVarianceReduction = fPm->GetBooleanParameter("Vr/UseG4DNAVarianceReduction");
            }
            
            if ( eIonisationModel == "emfietzoglou" ) {
                //G4DNAIonisation*
                G4VEmModel* modelA = new G4DNAEmfietzoglouIonisationModel();
                G4VEmModel* modelB = new G4DNABornIonisationModel();
                
                modelA->SetActivationLowEnergyLimit(10*eV);
                modelA->SetActivationHighEnergyLimit(10*keV);
                modelA->SetLowEnergyLimit(10*eV);
                modelA->SetHighEnergyLimit(10*keV);
                
                modelB->SetActivationLowEnergyLimit(10*keV);
                modelB->SetActivationHighEnergyLimit(1*MeV);
                modelB->SetLowEnergyLimit(10*keV);
                modelB->SetHighEnergyLimit(1*MeV);
                
                G4DNAIonisation* theDNAIonisationProcess = new G4DNAIonisation("e-_G4DNAIonisation");
                theDNAIonisationProcess->AddEmModel(1,modelA);
                theDNAIonisationProcess->AddEmModel(2,modelB);

                if ( !useVarianceReduction ) {
                    ph->RegisterProcess(theDNAIonisationProcess, particle);
                } else {
                    G4int numberOfSplit = fPm->GetIntegerParameter("Vr/ParticleSplit/NumberOfSplit");
                    G4String splitRegion = fPm->GetStringParameter("Vr/ParticleSplit/SplitElectronsInRegionNamed");
                
                    G4cout << "-- Secondary split for electrons created in ionisation process actived "
                    << "with split number " << numberOfSplit << "-- " << G4endl;
                
                    //G4ProcessManager* eman = particle->GetProcessManager();
                    TsSplitProcessG4DNA* splitProcess = new TsSplitProcessG4DNA(splitRegion, numberOfSplit);
                    splitProcess->RegisterProcess(theDNAIonisationProcess);
                    //eman->AddDiscreteProcess(splitProcess);
                    ph->RegisterProcess(splitProcess, particle);
                }
            } else if ( eIonisationModel == "cpa100") {
                G4VEmModel* modelA = new G4DNACPA100IonisationModel();
                G4VEmModel* modelB = new G4DNABornIonisationModel();
                
                modelA->SetActivationLowEnergyLimit(11*eV);
                modelA->SetActivationHighEnergyLimit(255.955*keV);
                modelA->SetLowEnergyLimit(11*eV);
                modelA->SetHighEnergyLimit(255.955*keV);
                
                modelB->SetActivationLowEnergyLimit(11*eV);
                modelB->SetActivationHighEnergyLimit(255.955*keV);
                modelB->SetLowEnergyLimit(255.955*keV);
                modelB->SetHighEnergyLimit(1*MeV);
                
                G4DNAIonisation* theDNAIonisationProcess = new G4DNAIonisation("e-_G4DNAIonisation");

                theDNAIonisationProcess->AddEmModel(1,modelA);
                theDNAIonisationProcess->AddEmModel(2,modelB);
                
                if ( !useVarianceReduction ) {
                    ph->RegisterProcess(theDNAIonisationProcess, particle);
                } else {
                    G4int numberOfSplit = fPm->GetIntegerParameter("Vr/ParticleSplit/NumberOfSplit");
                    G4String splitRegion = fPm->GetStringParameter("Vr/ParticleSplit/SplitElectronsInRegionNamed");
                
                    G4cout << "-- Secondary split for electrons created in ionisation process actived "
                    << "with split number " << numberOfSplit << "-- " << G4endl;
                
                    G4ProcessManager* eman = particle->GetProcessManager();
                    TsSplitProcessG4DNA* splitProcess = new TsSplitProcessG4DNA(splitRegion, numberOfSplit);
                    splitProcess->RegisterProcess(theDNAIonisationProcess);
                    eman->AddDiscreteProcess(splitProcess);
                }
            } else if ( eIonisationModel.contains("dna") ) {
                if ( eIonisationModel == "dnaemfietzoglou" ) {
                    G4DNAEmfietzoglouIonisationModel* e_modelDNAEmfietzoglouIonisation = new G4DNAEmfietzoglouIonisationModel();
                    G4DNAPTBIonisationModel* modelDNAPTBIonisation = new G4DNAPTBIonisationModel("THF/TMP/PY",particle);
                    
                    G4DNAModelInterface* e_ionisationInteraction = new G4DNAModelInterface("e-_ionisation_interaction");
                    e_ionisationInteraction->RegisterModel(e_modelDNAEmfietzoglouIonisation,particle);
                    e_ionisationInteraction->RegisterModel(modelDNAPTBIonisation);
                    e_ionisationInteraction->RegisterModel(new G4DNAVacuumModel());
                    
                    G4DNAIonisation* theDNAIonisationProcess = new G4DNAIonisation("e-_G4DNAIonisation");

                    theDNAIonisationProcess->SetEmModel(e_ionisationInteraction);
                    
                    if ( !useVarianceReduction ) {
                        ph->RegisterProcess(theDNAIonisationProcess, particle);
                    } else {
                        G4int numberOfSplit = fPm->GetIntegerParameter("Vr/ParticleSplit/NumberOfSplit");
                        G4String splitRegion = fPm->GetStringParameter("Vr/ParticleSplit/SplitElectronsInRegionNamed");
                
                        G4cout << "-- Secondary split for electrons created in ionisation process actived "
                        << "with split number " << numberOfSplit << "-- " << G4endl;
                
                        G4ProcessManager* eman = particle->GetProcessManager();
                        TsSplitProcessG4DNA* splitProcess = new TsSplitProcessG4DNA(splitRegion, numberOfSplit);
                        splitProcess->RegisterProcess(theDNAIonisationProcess);
                        eman->AddDiscreteProcess(splitProcess);
                    }
                    
                } else if ( eIonisationModel == "dnaborn" ) {
                    G4DNABornIonisationModel* e_modelDNABornIonisation = new G4DNABornIonisationModel();
                    G4DNAPTBIonisationModel* modelDNAPTBIonisation = new G4DNAPTBIonisationModel("THF/TMP/PY",particle);
                    
                    G4DNAModelInterface* e_ionisationInteraction = new G4DNAModelInterface("e-_ionisation_interaction");
                    e_ionisationInteraction->RegisterModel(e_modelDNABornIonisation,particle);
                    e_ionisationInteraction->RegisterModel(modelDNAPTBIonisation);
                    e_ionisationInteraction->RegisterModel(new G4DNAVacuumModel());
                    
                    G4DNAIonisation* theDNAIonisationProcess = new G4DNAIonisation("e-_G4DNAIonisation");

                    theDNAIonisationProcess->SetEmModel(e_ionisationInteraction);
                    
                    if ( !useVarianceReduction ) {
                        ph->RegisterProcess(theDNAIonisationProcess, particle);
                    } else {
                        G4int numberOfSplit = fPm->GetIntegerParameter("Vr/ParticleSplit/NumberOfSplit");
                        G4String splitRegion = fPm->GetStringParameter("Vr/ParticleSplit/SplitElectronsInRegionNamed");
                
                        G4cout << "-- Secondary split for electrons created in ionisation process actived "
                        << "with split number " << numberOfSplit << "-- " << G4endl;
                
                        G4ProcessManager* eman = particle->GetProcessManager();
                        TsSplitProcessG4DNA* splitProcess = new TsSplitProcessG4DNA(splitRegion, numberOfSplit);
                        splitProcess->RegisterProcess(theDNAIonisationProcess);
                        eman->AddDiscreteProcess(splitProcess);
                    }
                }
            } else {
                G4DNAIonisation* theDNAIonisationProcess = new G4DNAIonisation("e-_G4DNAIonisation");
		G4DNABornIonisationModel* mod = new G4DNABornIonisationModel();
                mod->SetHighEnergyLimit(1*MeV);
		mod->SelectFasterComputation(true);
		theDNAIonisationProcess->SetEmModel(mod);

                if ( !useVarianceReduction ) {
                    ph->RegisterProcess(theDNAIonisationProcess, particle);
                } else {
                    G4int numberOfSplit = fPm->GetIntegerParameter("Vr/ParticleSplit/NumberOfSplit");
                    G4String splitRegion = fPm->GetStringParameter("Vr/ParticleSplit/SplitElectronsInRegionNamed");
                
                    G4cout << "-- Secondary split for electrons created in ionisation process actived "
                    << "with split number " << numberOfSplit << "-- " << G4endl;
                
                    //G4ProcessManager* eman = particle->GetProcessManager();
                    TsSplitProcessG4DNA* splitProcess = new TsSplitProcessG4DNA(splitRegion, numberOfSplit);
                    splitProcess->RegisterProcess(theDNAIonisationProcess);
                    //eman->AddDiscreteProcess(splitProcess);
                    ph->RegisterProcess(splitProcess, particle);
                }
            }
            
            G4DNAElectronSolvation* solvation = new G4DNAElectronSolvation("e-_G4DNAElectronSolvation");
            
            auto solvationModel = G4DNASolvationModelFactory::GetMacroDefinedModel();
            if ( solvationHighLimit > 0 )
                solvationModel->SetHighEnergyLimit(solvationHighLimit);
            
            if ( fPm->ParameterExists("Ph/"+fName+"/Electron/SetHighEnergyLimitForSolvation") ) {
                solvationHighLimit = fPm->GetDoubleParameter("Ph/"+fName+"/Electron/SetHighEnergyLimitForSolvation","Energy");
                solvationModel->SetHighEnergyLimit(solvationHighLimit);
            }
            solvation->SetEmModel(solvationModel);
            
            ph->RegisterProcess(solvation, particle);
            
            // *** Vibrational excitation ***
            G4bool activeVibExcitation = true;
            if ( fPm->ParameterExists("Ph/"+fName+"/Electron/ActiveVibExcitation") )
                activeVibExcitation = fPm->GetBooleanParameter("Ph/"+fName+"/Electron/ActiveVibExcitation");
            
            if ( activeVibExcitation )
                ph->RegisterProcess(new G4DNAVibExcitation("e-_G4DNAVibExcitation"), particle);
            
            // *** Attachment ***
            G4bool activeAttachment = true;
            if ( fPm->ParameterExists("Ph/"+fName+"/Electron/ActiveAttachment") )
                activeAttachment = fPm->GetBooleanParameter("Ph/"+fName+"/Electron/ActiveAttachment");
            
            if ( activeAttachment )
                ph->RegisterProcess(new G4DNAAttachment("e-_G4DNAAttachment"), particle);
            
        } else if ( particleName == "proton" ) {
            G4String pScatteringModel = "default";
            if ( fPm->ParameterExists("Ph/"+fName+"/Proton/SetElasticScatteringModel") )
                pScatteringModel = fPm->GetStringParameter("Ph/"+fName+"/Proton/SetElasticScatteringModel");
            
            pScatteringModel.toLower();
            if ( pScatteringModel == "wentzelvi" ) {
                G4hMultipleScattering* msc = new G4hMultipleScattering();
                msc->SetEmModel(new G4LowEWentzelVIModel());
                ph->RegisterProcess(msc, particle);
            } else { 
                ph->RegisterProcess(new G4DNAElastic("proton_G4DNAElastic"), particle);
            }
            
            ph->RegisterProcess(new G4DNAExcitation("proton_G4DNAExcitation"), particle);

            G4String pIonisationModel = "default";
            if ( fPm->ParameterExists("Ph/"+fName+"/Proton/SetIonisationModel") )
                pIonisationModel = fPm->GetStringParameter("Ph/"+fName+"/Proton/SetIonisationModel");
            pIonisationModel.toLower();

            if ( pIonisationModel == "ritracks") {
                G4DNAIonisation* protonIonisationProcess = new G4DNAIonisation("proton_G4DNAIonisation");
                G4DNARuddIonisationExtendedModel* mod1 = new G4DNARuddIonisationExtendedModel();
                mod1->SetLowEnergyLimit(0*eV);
                mod1->SetHighEnergyLimit(500*keV);

                TsDNARuddIonisationExtendedRITRACKSModel* mod2 = new TsDNARuddIonisationExtendedRITRACKSModel();
                mod2->SetLowEnergyLimit(500*keV);
                mod2->SetHighEnergyLimit(500*MeV);

                protonIonisationProcess->AddEmModel(1, mod1);
                protonIonisationProcess->AddEmModel(2, mod2);
            }

            else {
                G4DNAIonisation* protonIonisationProcess = new G4DNAIonisation("proton_G4DNAIonisation");
                G4DNARuddIonisationExtendedModel* mod1 = new G4DNARuddIonisationExtendedModel();
                mod1->SetLowEnergyLimit(0*eV);
                mod1->SetHighEnergyLimit(500*keV);
    			
                G4DNABornIonisationModel* mod2;
                mod2= new G4DNABornIonisationModel();
                mod2->SetLowEnergyLimit(500*keV);
                mod2->SetHighEnergyLimit(100*MeV);
    	        mod2->SelectFasterComputation(true);
                
                TsDNARuddIonisationExtendedModel* mod3 = new TsDNARuddIonisationExtendedModel();
                mod3->SetLowEnergyLimit(100*MeV);
                mod3->SetHighEnergyLimit(500*MeV);
    			
                protonIonisationProcess->AddEmModel(1, mod1);
                protonIonisationProcess->AddEmModel(2, mod2);
                protonIonisationProcess->AddEmModel(3, mod3);
                ph->RegisterProcess(protonIonisationProcess, particle);
                
                ph->RegisterProcess(new G4DNAChargeDecrease("proton_G4DNAChargeDecrease"), particle);
            }
            
        } else if ( particleName == "hydrogen" ) {
            ph->RegisterProcess(new G4DNAElastic("hydrogen_G4DNAElastic"), particle);
            ph->RegisterProcess(new G4DNAExcitation("hydrogen_G4DNAExcitation"), particle);
	    G4DNAIonisation* theDNAIonisationProcess = new G4DNAIonisation("hydrogen_G4DNAIonisation");
            theDNAIonisationProcess->SetEmModel(new G4DNARuddIonisationExtendedModel());
            ph->RegisterProcess(theDNAIonisationProcess, particle);
            ph->RegisterProcess(new G4DNAChargeIncrease("hydrogen_G4DNAChargeIncrease"), particle);
            
        } else if ( particleName == "alpha" ) {
            G4String alphaEModel = "default";
            if ( fPm->ParameterExists("Ph/"+fName+"/Alpha/SetElasticScatteringModel") ) 
                alphaEModel = fPm->ParameterExists("Ph/"+fName+"/Alpha/SetElasticScatteringModel");
            
            alphaEModel.toLower();

            if (alphaEModel == "default") { 
                ph->RegisterProcess(new G4DNAElastic("alpha_G4DNAElastic"), particle); 
            } else {
                G4hMultipleScattering* msc = new G4hMultipleScattering();
                msc->SetEmModel(new G4LowEWentzelVIModel());
                ph->RegisterProcess(msc, particle);
            }

            ph->RegisterProcess(new G4DNAExcitation("alpha_G4DNAExcitation"), particle);
            G4DNAIonisation* theDNAIonisationProcess = new G4DNAIonisation("alpha_G4DNAIonisation"); 
            theDNAIonisationProcess->SetEmModel(new G4DNARuddIonisationExtendedModel());
            ph->RegisterProcess(theDNAIonisationProcess, particle);
            ph->RegisterProcess(new G4DNAChargeDecrease("alpha_G4DNAChargeDecrease"), particle);
            
        } else if ( particleName == "alpha+" ) {
            G4String alphaEModel = "default";
            if ( fPm->ParameterExists("Ph/"+fName+"/AlphaPlus/SetElasticScatteringModel") ) 
                alphaEModel = fPm->ParameterExists("Ph/"+fName+"/AlphaPlus/SetElasticScatteringModel");
            
            alphaEModel.toLower();
            if (alphaEModel == "default") { 
                ph->RegisterProcess(new G4DNAElastic("alpha+_G4DNAElastic"), particle); 
 
 	    } else { 
                G4hMultipleScattering* msc = new G4hMultipleScattering();
                msc->SetEmModel(new G4LowEWentzelVIModel());
                ph->RegisterProcess(msc, particle);
            } 

            ph->RegisterProcess(new G4DNAExcitation("alpha+_G4DNAExcitation"), particle);
            G4DNAIonisation* theDNAIonisationProcess = new G4DNAIonisation("alpha+_G4DNAIonisation");
            theDNAIonisationProcess->SetEmModel(new G4DNARuddIonisationExtendedModel());
            ph->RegisterProcess(theDNAIonisationProcess, particle);
            ph->RegisterProcess(new G4DNAChargeDecrease("alpha+_G4DNAChargeDecrease"), particle);
            ph->RegisterProcess(new G4DNAChargeIncrease("alpha+_G4DNAChargeIncrease"), particle);
            
        } else if ( particleName == "helium" ) {
            ph->RegisterProcess(new G4DNAElastic("helium_G4DNAElastic"), particle);
            ph->RegisterProcess(new G4DNAExcitation("helium_G4DNAExcitation"), particle);
            G4DNAIonisation* theDNAIonisationProcess = new G4DNAIonisation("helium_G4DNAIonisation");
            theDNAIonisationProcess->SetEmModel(new G4DNARuddIonisationExtendedModel());
            ph->RegisterProcess(theDNAIonisationProcess, particle);
            ph->RegisterProcess(new G4DNAChargeIncrease("helium_G4DNAChargeIncrease"), particle);
            
        } else if ( particleName == "GenericIon" ) {
            G4String ionEModel = "default";
            if ( fPm->ParameterExists("Ph/"+fName+"/GenericIon/SetElasticScatteringModel") ) 
                ionEModel = fPm->GetStringParameter("Ph/"+fName+"/GenericIon/SetElasticScatteringModel");
             
            ionEModel.toLower();

            if ( ionEModel == "wentzelvi" ) { 
                G4hMultipleScattering* msc = new G4hMultipleScattering();
                msc->SetEmModel(new G4LowEWentzelVIModel(), 1);
                ph->RegisterProcess(msc, particle);
            }
 
	    G4DNAIonisation* genericIonIonisationProcess = new G4DNAIonisation("GenericIon_G4DNAIonisation"); 
            if ( fPm->ParameterExists("Ph/"+fName+"/GenericIon/SetInelasticScatteringModel")) {
                ionEModel = fPm->GetStringParameter("Ph/"+fName+"/GenericIon/SetInelasticScatteringModel");
            }
           
            ionEModel.toLower();
            if ( ionEModel == "ruddextended" ) { 
	        TsDNARuddIonisationExtendedModel* mod1 = new TsDNARuddIonisationExtendedModel();
                genericIonIonisationProcess->AddEmModel(1, mod1);
                //mod3->SetLowEnergyLimit(100*eV);
                //mod3->SetHighEnergyLimit(1000*MeV);
	        //mod3->SetLowEnergyLimit(100*MeV);
            }

            ph->RegisterProcess(genericIonIonisationProcess, particle);
            
        } else if (particleName == "e+") {
            G4eMultipleScattering* msc = new G4eMultipleScattering();
            msc->SetStepLimitType(fUseDistanceToBoundary);
            G4eIonisation* eIoni = new G4eIonisation();
            eIoni->SetStepFunction(0.2, 100*um);
            
            ph->RegisterProcess(msc, particle);
            ph->RegisterProcess(eIoni, particle);
            ph->RegisterProcess(new G4eBremsstrahlung(), particle);
            ph->RegisterProcess(new G4eplusAnnihilation(), particle);
            
        } else if (particleName == "gamma") {
            G4PhotoElectricEffect* thePhotoElectricEffect = new G4PhotoElectricEffect();
            thePhotoElectricEffect->SetEmModel(new G4LivermorePhotoElectricModel());
            ph->RegisterProcess(thePhotoElectricEffect, particle);
            
            G4ComptonScattering* theComptonScattering = new G4ComptonScattering();
            theComptonScattering->SetEmModel(new G4LivermoreComptonModel());
            ph->RegisterProcess(theComptonScattering, particle);
            
            G4GammaConversion* theGammaConversion = new G4GammaConversion();
            theGammaConversion->SetEmModel(new G4LivermoreGammaConversionModel());
            ph->RegisterProcess(theGammaConversion, particle);
            
            G4RayleighScattering* theRayleigh = new G4RayleighScattering();
            ph->RegisterProcess(theRayleigh, particle);
        }
    }
    
    G4VAtomDeexcitation* de = new G4UAtomicDeexcitation();
    G4LossTableManager::Instance()->SetAtomDeexcitation(de);
}

