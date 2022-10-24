// Physics Module for TsEmDNAPhysics_opt1

#include "TsEmDNAPhysics_opt1.hh"
#include "TsParameterManager.hh"

#include "TsDNARuddIonisationExtendedModel.hh"
#include "TsDNARuddIonisationExtendedRITRACKSModel.hh"
#include "TsSplitProcessG4DNA.hh"

#include "G4SystemOfUnits.hh"
#include "G4DNAGenericIonsManager.hh"

#include "G4DNAElectronSolvation.hh"
#include "G4DNAElastic.hh"
#include "G4DNAELSEPAElasticModel.hh"
#include "G4DNAChampionElasticModel.hh"
#include "G4DNAScreenedRutherfordElasticModel.hh"

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
#include "G4ProcessManager.hh"

G4_DECLARE_PHYSCONSTR_FACTORY(TsEmDNAPhysics_opt1);

TsEmDNAPhysics_opt1::TsEmDNAPhysics_opt1(G4int ver, const G4String&)
: G4VPhysicsConstructor("TsEmDNAPhysics_opt1"), verbose(ver)
{
    G4EmParameters* param = G4EmParameters::Instance();
    param->SetDefaults();
    param->SetFluo(true);
    param->SetAuger(true);
    param->SetAugerCascade(true);
    param->SetDeexcitationIgnoreCut(true);
    SetPhysicsType(bElectromagnetic);
}


TsEmDNAPhysics_opt1::TsEmDNAPhysics_opt1(TsParameterManager* pM)
: G4VPhysicsConstructor("TsEmDNAPhysics_opt1"), fPm(pM), verbose(2)
{
    G4EmParameters* param = G4EmParameters::Instance();
    param->SetDefaults();
    param->SetFluo(true);
    param->SetAuger(true);
    param->SetAugerCascade(true);
    param->SetDeexcitationIgnoreCut(true);
    SetPhysicsType(bElectromagnetic);
}


TsEmDNAPhysics_opt1::~TsEmDNAPhysics_opt1()
{}


void TsEmDNAPhysics_opt1::ConstructParticle()
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


void TsEmDNAPhysics_opt1::ConstructProcess()
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
            if ( fPm->ParameterExists("Ph/Default/Electron/SetElasticScatteringModel"))
                eScatteringModel = fPm->GetStringParameter("Ph/Default/Electron/SetElasticScatteringModel");
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
		solvationHighLimit = 7.4 * eV;

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
            if ( fPm->ParameterExists("Ph/Default/Electron/SetExcitationModel") )
                eExcitationModel = fPm->GetStringParameter("Ph/Default/Electron/SetExcitationModel");
            eExcitationModel.toLower();
            
            G4String eIonisationModel = "born";
            if ( fPm->ParameterExists("Ph/Default/Electron/SetIonisationModel") )
                eIonisationModel = fPm->GetStringParameter("Ph/Default/Electron/SetIonisationModel");
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
                
                    G4ProcessManager* eman = particle->GetProcessManager();
                    TsSplitProcessG4DNA* splitProcess = new TsSplitProcessG4DNA(splitRegion, numberOfSplit);
                    splitProcess->RegisterProcess(theDNAIonisationProcess);
                    eman->AddDiscreteProcess(splitProcess);
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
                //G4DNAIonisation*
                G4DNAIonisation* theDNAIonisationProcess = new G4DNAIonisation("e-_G4DNAIonisation");
		G4DNABornIonisationModel* mod = new G4DNABornIonisationModel();
		mod->SelectFasterComputation(true);
		theDNAIonisationProcess->SetEmModel(mod);

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
            
            G4DNAElectronSolvation* solvation = new G4DNAElectronSolvation("e-_G4DNAElectronSolvation");
            
            G4VEmModel* solvationModel = new G4DNAOneStepThermalizationModel();
            if ( solvationHighLimit > 0 )
                solvationModel->SetHighEnergyLimit(solvationHighLimit);
            
            if ( fPm->ParameterExists("Ph/Default/Electron/SetHighEnergyLimitForSolvation") ) {
                solvationHighLimit = fPm->GetDoubleParameter("Ph/Default/Electron/SetHighEnergyLimitForSolvation","Energy");
                solvationModel->SetHighEnergyLimit(solvationHighLimit);
            }
            solvation->SetEmModel(solvationModel);
            
            ph->RegisterProcess(solvation, particle);
            
            // *** Vibrational excitation ***
            G4bool activeVibExcitation = true;
            if ( fPm->ParameterExists("Ph/Default/Electron/ActiveVibExcitation") )
                activeVibExcitation = fPm->GetBooleanParameter("Ph/Default/Electron/ActiveVibExcitation");
            
            if ( activeVibExcitation )
                ph->RegisterProcess(new G4DNAVibExcitation("e-_G4DNAVibExcitation"), particle);
            
            // *** Attachment ***
            G4bool activeAttachment = true;
            if ( fPm->ParameterExists("Ph/Default/Electron/ActiveAttachment") )
                activeAttachment = fPm->GetBooleanParameter("Ph/Default/Electron/ActiveAttachment");
            
            if ( activeAttachment )
                ph->RegisterProcess(new G4DNAAttachment("e-_G4DNAAttachment"), particle);
            
        } else if ( particleName == "proton" ) {
			G4bool protonStationaryPhysics = false;
			if ( fPm->ParameterExists("Ph/Default/Proton/ActivePhysicsStationary"))
				protonStationaryPhysics = fPm->GetBooleanParameter("Ph/Default/Proton/ActivePhysicsStationary");
			
            G4String pScatteringModel = "default"; //wentzelvi";
            if ( fPm->ParameterExists("Ph/Default/Proton/SetElasticScatteringModel") )
                pScatteringModel = fPm->GetStringParameter("Ph/Default/Proton/SetElasticScatteringModel");
            
            pScatteringModel.toLower();
            if ( pScatteringModel == "wentzelvi" ) {
                G4hMultipleScattering* msc = new G4hMultipleScattering();
                msc->SetEmModel(new G4LowEWentzelVIModel());
                ph->RegisterProcess(msc, particle);
            } else {
				G4DNAElastic* theDNAElasticProcess = new G4DNAElastic("proton_G4DNAElastic");
				theDNAElasticProcess->SetEmModel(new G4DNAIonElasticModel());
				((G4DNAIonElasticModel*)(theDNAElasticProcess->EmModel()))->SelectStationary(protonStationaryPhysics);
                ph->RegisterProcess(theDNAElasticProcess, particle);
            }
			
			G4DNAExcitation* theDNAExcitationProcess = new G4DNAExcitation("proton_G4DNAExcitation");
			theDNAExcitationProcess->SetEmModel(new G4DNAMillerGreenExcitationModel());
			theDNAExcitationProcess->SetEmModel(new G4DNABornExcitationModel());
			
			((G4DNAMillerGreenExcitationModel*)(theDNAExcitationProcess->EmModel()))->SetLowEnergyLimit(10*eV);
			((G4DNAMillerGreenExcitationModel*)(theDNAExcitationProcess->EmModel()))->SetHighEnergyLimit(500*keV);
			((G4DNAMillerGreenExcitationModel*)(theDNAExcitationProcess->EmModel()))->SelectStationary(protonStationaryPhysics);
			
			((G4DNABornExcitationModel*)(theDNAExcitationProcess->EmModel(1)))->SetLowEnergyLimit(500*keV);
			((G4DNABornExcitationModel*)(theDNAExcitationProcess->EmModel(1)))->SetHighEnergyLimit(500*MeV);
			((G4DNABornExcitationModel*)(theDNAExcitationProcess->EmModel(1)))->SelectStationary(protonStationaryPhysics);
			
            ph->RegisterProcess(theDNAExcitationProcess, particle);
            
            G4DNAIonisation* theDNAIonisationProcess = new G4DNAIonisation("proton_G4DNAIonisation");
			if ( !protonStationaryPhysics ) {
				G4String pioniModel = "geant4-dna";
				if ( fPm->ParameterExists("Ph/Default/Proton/SetIonisationCrossSection") )
					pioniModel = fPm->GetStringParameter("Ph/Default/Proton/SetIonisationCrossSection");

				pioniModel.toLower();
				if ( pioniModel == "geant4-dna" ) {
					G4VEmModel* mod1 = new G4DNARuddIonisationModel();
					mod1->SetLowEnergyLimit(0*eV);
					mod1->SetHighEnergyLimit(500*keV);
					
					G4VEmModel* mod2;
					mod2= new G4DNABornIonisationModel();
					mod2->SetLowEnergyLimit(500*keV);
					mod2->SetHighEnergyLimit(100*MeV);
					theDNAIonisationProcess->AddEmModel(1, mod1);
					theDNAIonisationProcess->AddEmModel(2, mod2);
				} else if (pioniModel == "ritracks" ) {	
					G4VEmModel* mod1 = new TsDNARuddIonisationExtendedRITRACKSModel();
					mod1->SetLowEnergyLimit(0*eV); 
					mod1->SetHighEnergyLimit(500*MeV);
					theDNAIonisationProcess->AddEmModel(1, mod1);
				} else {
					G4cerr << "TOPAS is exiting due to a serious in physics list setup" << G4endl;
					G4cerr << "Ph/Default/Proton/SetIonisationCrossSection has an invalid parameter value" << G4endl; 
					G4cerr << "Options: Geant4-DNA, RITRACKS " << G4endl;
					fPm->AbortSession(1);
				}

				ph->RegisterProcess(theDNAIonisationProcess, particle);

			} else {
				theDNAIonisationProcess->SetEmModel(new G4DNARuddIonisationExtendedModel);
				theDNAIonisationProcess->SetEmModel(new G4DNABornIonisationModel);
				
				((G4DNARuddIonisationExtendedModel*)(theDNAIonisationProcess->EmModel()))->SetLowEnergyLimit(0*eV);
				((G4DNARuddIonisationExtendedModel*)(theDNAIonisationProcess->EmModel()))->SetHighEnergyLimit(500*keV);
				((G4DNARuddIonisationExtendedModel*)(theDNAIonisationProcess->EmModel()))->SelectStationary(protonStationaryPhysics);
				
				((G4DNABornIonisationModel*)(theDNAIonisationProcess->EmModel(1)))->SetLowEnergyLimit(500*keV);
				((G4DNABornIonisationModel*)(theDNAIonisationProcess->EmModel(1)))->SetHighEnergyLimit(100*MeV);
				((G4DNABornIonisationModel*)(theDNAIonisationProcess->EmModel(1)))->SelectStationary(protonStationaryPhysics);
				
				ph->RegisterProcess(theDNAIonisationProcess, particle);
			}
			
			G4DNAChargeDecrease* theDNAChargeDecreaseProcess = new G4DNAChargeDecrease("proton_G4DNAChargeDecrease");
			theDNAChargeDecreaseProcess->SetEmModel(new G4DNADingfelderChargeDecreaseModel());
			((G4DNADingfelderChargeDecreaseModel*)(theDNAChargeDecreaseProcess->EmModel()))->SelectStationary(protonStationaryPhysics);
            ph->RegisterProcess(new G4DNAChargeDecrease("proton_G4DNAChargeDecrease"), particle);
            
        } else if ( particleName == "hydrogen" ) {
	    ph->RegisterProcess(new G4DNAElastic("hydrogen_G4DNAElastic"), particle);
            ph->RegisterProcess(new G4DNAExcitation("hydrogen_G4DNAExcitation"), particle);
            ph->RegisterProcess(new G4DNAIonisation("hydrogen_G4DNAIonisation"), particle);
            ph->RegisterProcess(new G4DNAChargeIncrease("hydrogen_G4DNAChargeIncrease"), particle);
            
        } else if ( particleName == "alpha" ) {
            
            G4hMultipleScattering* msc = new G4hMultipleScattering();
            msc->SetEmModel(new G4LowEWentzelVIModel());
            ph->RegisterProcess(msc, particle);

	    //ph->RegisterProcess(new G4DNAElastic("alpha_G4DNAElastic"), particle);            
            ph->RegisterProcess(new G4DNAExcitation("alpha_G4DNAExcitation"), particle);
            ph->RegisterProcess(new G4DNAIonisation("alpha_G4DNAIonisation"), particle);
            ph->RegisterProcess(new G4DNAChargeDecrease("alpha_G4DNAChargeDecrease"), particle);
            
        } else if ( particleName == "alpha+" ) {
            
            G4hMultipleScattering* msc = new G4hMultipleScattering();
            msc->SetEmModel(new G4LowEWentzelVIModel());
            ph->RegisterProcess(msc, particle);
            //ph->RegisterProcess(new G4DNAElastic("alpha+_G4DNAElastic"), particle); 
            ph->RegisterProcess(new G4DNAExcitation("alpha+_G4DNAExcitation"), particle);
            ph->RegisterProcess(new G4DNAIonisation("alpha+_G4DNAIonisation"), particle);
            ph->RegisterProcess(new G4DNAChargeDecrease("alpha+_G4DNAChargeDecrease"), particle);
            ph->RegisterProcess(new G4DNAChargeIncrease("alpha+_G4DNAChargeIncrease"), particle);
            
        } else if ( particleName == "helium" ) {
            ph->RegisterProcess(new G4DNAExcitation("helium_G4DNAExcitation"), particle);
            ph->RegisterProcess(new G4DNAIonisation("helium_G4DNAIonisation"), particle);
            ph->RegisterProcess(new G4DNAChargeIncrease("helium_G4DNAChargeIncrease"), particle);
            
        } else if ( particleName == "GenericIon" ) {
           
            if ( fPm->ParameterExists("Ph/Default/Ion/UseIonScatteringModelNamed" ) ) {
                G4String ionMSC = fPm->GetStringParameter("Ph/Default/Ion/UseIonScatteringModelNamed");
                ionMSC.toLower();
                if ( ionMSC == "wentzelvi" ) {
                    G4hMultipleScattering* msc = new G4hMultipleScattering();
                    msc->SetEmModel(new G4LowEWentzelVIModel(), 1);
                    ph->RegisterProcess(msc, particle);
                } else if ( ionMSC == "coulomb" ) {
                    G4hMultipleScattering* msc = new G4hMultipleScattering("ionmsc");
                    ph->RegisterProcess(msc, particle);
                } else {
                    G4cerr << "Error! TOPAS has stopped due an error in physics list" << G4endl;
                    G4cerr << "Parameter Ph/Default/Ion/UseIonScatteringModelNamed has wrong model name" << G4endl;
                    G4cerr << "Available models are: None, WentzelVI and Coulomb" << G4endl;
                    exit(1);
                }
            } 

	    G4DNAIonisation* genericIonIonisationProcess = new G4DNAIonisation("GenericIon_G4DNAIonisation"); 
	        G4VEmModel* mod = new TsDNARuddIonisationExtendedModel();
            genericIonIonisationProcess->AddEmModel(1, mod);

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

