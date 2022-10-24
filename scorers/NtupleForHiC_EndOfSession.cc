// EndSession for TOPAS
/*
*
*  Ingram & Henthorn et al
*  Hi-C TOPAS Geometry
*
*/


#include "NtupleForHiC_EndOfSession.hh"

NtupleForHiC_EndOfSession::NtupleForHiC_EndOfSession(TsParameterManager* pM)
{
    G4bool ActiveHiCScorer = false;

    if (pM->ParameterExists("Sc/HiCScorer/Active"))
        ActiveHiCScorer = pM->GetBooleanParameter("Sc/HiCScorer/Active");

    if (ActiveHiCScorer)
        Start(pM);
}

NtupleForHiC_EndOfSession::~NtupleForHiC_EndOfSession(){}

G4bool NtupleForHiC_EndOfSession::Start(TsParameterManager* pM)
{
    if (pM->ParameterExists("Sc/HiCScorer/WriteSSB")){
        WriteSSB=pM->GetBooleanParameter("Sc/HiCScorer/WriteSSB");
    }
    
    G4String SDDFilename="Output.sdd";
    if (pM->ParameterExists("Sc/HiCScorer/SDDOutputFile")){
        SDDFilename=pM->GetStringParameter("Sc/HiCScorer/SDDOutputFile");
    }
    
    G4String RawDamageFile="Output.phsp";
    if (pM->ParameterExists("Sc/HiCScorer/OutputFile")){
        RawDamageFile=pM->GetStringParameter("Sc/HiCScorer/OutputFile");
        RawDamageFile+=".phsp";
    }
    
    G4String VerticesFile="HiC_Vertices/4umIMR90_noLADs.vert.txt";
    if (pM->ParameterExists("Sc/HiCScorer/VerticesFile")){
        VerticesFile=pM->GetStringParameter("Ge/HiC/VerticesFile");
    }
    
    if (pM->ParameterExists("Sc/HiCScorer/OverwiteSDD")){
        OverwriteSDD=pM->GetBooleanParameter("Sc/HiCScorer/OverwiteSDD");
    }
    if (OverwriteSDD){remove(SDDFilename);}
    else {WriteSDDHeader=false;}
    
    
    std::cout<<"  ---  Session Ends, Begin Analysis"<<std::endl;
    
    //Read in raw hit data
    std::map<G4int,std::vector<NtupleForHiC_HitPoint*>> ExposureHits = ReadOutput(RawDamageFile);
    if (ExposureHits.size()==0){
        std::cout<<"   --   No Hits Scored"<<std::endl;
        return false;
    }
    
    //Cluster raw hits: sorted by 'exposure', then by 'clusterID'
    std::map<G4int,std::map<G4int,NtupleForHiC_Clusters*>> ExposureClusters = StartClustering(ExposureHits);
    
    //Add fraction along chromosome to breaks & calculate chromosome size
    std::map<G4int,G4int> ChromosomeSize;
    G4ThreeVector NucDim;
    AddChromosomeFractions(VerticesFile,
                           ChromosomeSize,
                           ExposureClusters,
                           NucDim);
    
    //Fill the phase space store with scored damages
    G4int DamageEntries=0;
    TsDamagePhaseSpaceStore*PSStore=FillPhaseSpaceStore(ChromosomeSize,
                                                        ExposureClusters,
                                                        DamageEntries);
    
    //Fill the phase space header
    FillPhaseSpaceHeader(ChromosomeSize,
                         PSStore,
                         pM,
                         NucDim,
                         DamageEntries);
    
    
    WriteSDD(SDDFilename,
             PSStore);
    
    return true;
}

std::map<int,std::vector<NtupleForHiC_HitPoint*>> NtupleForHiC_EndOfSession::ReadOutput(G4String filename)
{
    std::cout<<"   --   Reading Raw Data From: "<<filename<<std::endl;
    std::ifstream infile(filename);
    
    std::map<int,std::vector<NtupleForHiC_HitPoint*>> Hits;
    
    G4String line;
    while (getline(infile, line)){
        G4double Pos_x=0.0, Pos_y=0.0, Pos_z=0.0;
        G4int Cause=-1;
        G4int Strand=-1;
        G4int Chromosome=-1;
        G4int BeadID=-1;
        G4int EvtID=-1;
        G4int RunID=-1;
        std::stringstream stream(line.data());
        stream
        >>Pos_x
        >>Pos_y
        >>Pos_z
        >>Cause
        >>Strand
        >>Chromosome
        >>BeadID
        >>EvtID
        >>RunID;
        
        G4ThreeVector Pos(Pos_x*nm,Pos_y*nm,Pos_z*nm);
        
        NtupleForHiC_HitPoint*hit=new NtupleForHiC_HitPoint();
        hit->SetPosition(Pos);
        hit->SetStrand(Strand);
        hit->SetIsBack(true);
        hit->SetIsBase(false);
        hit->SetChromID(Chromosome);
        hit->SetVolumeCopy(BeadID);
        hit->SetEventID(EvtID);
        hit->SetRunID(RunID);
        
        if (Cause==0){hit->SetIsPhys(true);}
        if (Cause==1){hit->SetIsChem(true);}
        if (Cause==2){hit->SetIsPhys(true);hit->SetIsChem(true);}
        if (Cause==3){hit->SetIsChargeMig(true);}
        
        Hits[RunID].push_back(hit);
    }
    
    return Hits;
}


std::map<G4int,std::map<G4int,NtupleForHiC_Clusters*>> NtupleForHiC_EndOfSession::StartClustering(std::map<int,std::vector<NtupleForHiC_HitPoint*>> ExposureHits)
{
    
    std::map<G4int,std::map<G4int,NtupleForHiC_Clusters*>> ExposureClusters;
    
    for (auto data:ExposureHits){
        std::cout<<"   --   Running Clustering For Exposure: "<<data.first<<std::endl;
        NtupleForHiC_ClusterAlgo*clustering=new NtupleForHiC_ClusterAlgo();
        clustering->RunDBSCAN(data.second);
        ExposureClusters[data.first]=clustering->GetClusterMap();
        G4int nClust=ExposureClusters[data.first].size();
        G4int nSSB=ExposureClusters[data.first][0]->GetnBacks();
        G4int nDSB=nClust-1;
        
        //All isolated damages sorted into cluster 0
        //Create a new "cluster" for each of these
        //overwrite that "0" cluster that contained the isolated damages with a new isolated damage
        NtupleForHiC_Clusters*IsoDam=ExposureClusters[data.first][0];
        G4int nIso=IsoDam->GetHits().size();
        for (G4int i=0;i<nIso;i++){
            if (i==nIso-1){
                NtupleForHiC_Clusters*cluster=new NtupleForHiC_Clusters();
                cluster->AddPoint(IsoDam->GetHits()[i]);
                cluster->SetIsSSB(true);
                ExposureClusters[data.first][0]=cluster;
            }
            else {
                NtupleForHiC_Clusters*cluster=new NtupleForHiC_Clusters();
                cluster->AddPoint(IsoDam->GetHits()[i]);
                cluster->SetIsSSB(true);
                ExposureClusters[data.first][nClust]=cluster;
                nClust++;
            }
        }
        
        
        std::cout<<"    -    "<<data.second.size()<<" Hits, "<<nSSB<<" SSBs, "<<nDSB<<" DSBs"<<std::endl;
        delete clustering;
    }
    
    return ExposureClusters;
}

void NtupleForHiC_EndOfSession::AddChromosomeFractions(G4String HiCVerticesFile,
                                                       std::map<G4int,G4int> &ChromosomeSize,
                                                       std::map<G4int,std::map<G4int,NtupleForHiC_Clusters*>> &ExposureClusters,
                                                       G4ThreeVector &NucDim)
{
    //Read HiC vertices, calculate fraction along chromosome
    //change HiC vertices into map<map<
    //sorted by chrom, then beadID
    ChromosomeParser*chrom=new ChromosomeParser();
    chrom->ReadBeads(HiCVerticesFile);
    std::map<G4int,std::vector<ChromObj>> BeadDetails = chrom->GetBeads();
    ChromosomeSize = chrom->GetChromosomeSize();
    chrom->ResizeNucleus();
    NucDim=chrom->GetDimensions();
    delete chrom;
    
    std::map<G4int,std::map<G4int,ChromObj>> SortedBeadDetails;
    for (auto bead:BeadDetails){
        for (unsigned int i=0;i<bead.second.size();i++){
            SortedBeadDetails[bead.first][bead.second[i].VolumeCopyNum]=bead.second[i];
        }
    }

    
    std::cout<<"   --   Adding Fraction Along Chromosome for Damages"<<std::endl;
    for (auto exposure:ExposureClusters){
        for (auto dsb:exposure.second){
                //add fraction along chromosome here
                G4int HiCBeadCopy=dsb.second->GetHicBeadCopy();
                G4int Chromosome=dsb.second->GetChromID();
                G4double MaxFract=SortedBeadDetails[Chromosome][HiCBeadCopy].MaxFractAlong;
                G4double MinFract=SortedBeadDetails[Chromosome][HiCBeadCopy].MinFractAlong;
                
                G4double Rando = G4UniformRand();
                G4double FractAlong=MinFract+(Rando*(MaxFract-MinFract));
                G4bool Arm=SetChromosomeArm(Chromosome-1,
                                            FractAlong);

                ExposureClusters[exposure.first][dsb.first]->SetChromosomePosition(FractAlong);
                ExposureClusters[exposure.first][dsb.first]->SetChromArm(Arm);
                ExposureClusters[exposure.first][dsb.first]->SetChromatid(1);
        }
    }
    
    
    
    
}

G4bool NtupleForHiC_EndOfSession::SetChromosomeArm(G4int ChromID,
                                                   G4double FractionAlong)
{
    G4double FractionAtCentromere=-1.0;
    
    //Set the arm of damage
    //centromere position as fraction of chromosome (for chrom 1->22 + XY
    map<int,double> CentromereLocations;
    CentromereLocations[0]=123.4/249.0; //C1 copy 1
    CentromereLocations[1]=93.9/242.0; //C2 copy 1
    CentromereLocations[2]=90.9/198.0; //C3 copy 1
    CentromereLocations[3]=50.0/190.0; //C4 copy 1
    CentromereLocations[4]=48.8/182.0; //C5 copy 1
    CentromereLocations[5]=59.8/171.0; //C6 copy 1
    CentromereLocations[6]=60.1/159.0; //C7 copy 1
    CentromereLocations[7]=45.2/145.0; //C8 copy 1
    CentromereLocations[8]=43.0/138.0; //C9 copy 1
    CentromereLocations[9]=39.8/134.0; //C10 copy 1
    CentromereLocations[10]=53.4/135.0; //C11 copy 1
    CentromereLocations[11]=35.5/133.0; //C12 copy 1
    CentromereLocations[12]=17.7/114.0; //C13 copy 1
    CentromereLocations[13]=17.2/107.0; //C14 copy 1
    CentromereLocations[14]=19.0/102.0; //C15 copy 1
    CentromereLocations[15]=36.8/90.0; //C16 copy 1
    CentromereLocations[16]=25.1/83.0; //C17 copy 1
    CentromereLocations[17]=18.5/80.0; //C18 copy 1
    CentromereLocations[18]=26.2/59.0; //C19 copy 1
    CentromereLocations[19]=28.1/64.0; //C20 copy 1
    CentromereLocations[20]=12.0/47.0; //C21 copy 1
    CentromereLocations[21]=15.0/51.0; //C22 copy 1
    CentromereLocations[22]=61.0/156.0; //CX
    CentromereLocations[23]=123.4/249.0; //C1 copy 2
    CentromereLocations[24]=93.9/242.0; //C2 copy 2
    CentromereLocations[25]=90.9/198.0; //C3 copy 2
    CentromereLocations[26]=50.0/190.0; //C4 copy 2
    CentromereLocations[27]=48.8/182.0; //C5 copy 2
    CentromereLocations[28]=59.8/171.0; //C6 copy 2
    CentromereLocations[29]=60.1/159.0; //C7 copy 2
    CentromereLocations[30]=45.2/145.0; //C8 copy 2
    CentromereLocations[31]=43.0/138.0; //C9 copy 2
    CentromereLocations[32]=39.8/134.0; //C10 copy 2
    CentromereLocations[33]=53.4/135.0; //C11 copy 2
    CentromereLocations[34]=35.5/133.0; //C12 copy 2
    CentromereLocations[35]=17.7/114.0; //C13 copy 2
    CentromereLocations[36]=17.2/107.0; //C14 copy 2
    CentromereLocations[37]=19.0/102.0; //C15 copy 2
    CentromereLocations[38]=36.8/90.0; //C16 copy 2
    CentromereLocations[39]=25.1/83.0; //C17 copy 2
    CentromereLocations[40]=18.5/80.0; //C18 copy 2
    CentromereLocations[41]=26.2/59.0; //C19 copy 2
    CentromereLocations[42]=28.1/64.0; //C20 copy 2
    CentromereLocations[43]=12.0/47.0; //C21 copy 2
    CentromereLocations[44]=15.0/51.0; //C22 copy 2
    CentromereLocations[45]=10.4/57.0; //CY
    
    if (CentromereLocations.find(ChromID)!=CentromereLocations.end()){
        FractionAtCentromere = CentromereLocations[ChromID];
    }
    else{
        cout<<"Cant Find Chromosome: "<<ChromID<<", for arm"<<endl;
    }
    if (FractionAlong>FractionAtCentromere){
        return 1;
    } else {
        return 0;
    }
}


TsDamagePhaseSpaceStore* NtupleForHiC_EndOfSession::FillPhaseSpaceStore(std::map<G4int,G4int> ,
                                                                        std::map<G4int,std::map<G4int,NtupleForHiC_Clusters*>> &ExposureClusters,
                                                                        G4int &DamageEntries)
{
    TsDamagePhaseSpaceStore*PSStore=TsDamagePhaseSpaceStore::Instance();

    
    //Fill damage store
    
    //Get ptr to damage store from TsDamagePhaseSpaceStore
    TsDamagePhaseSpaceStore::dmgstore *DamageStore=&PSStore->GetDamageEventStore();
    for (auto exposure:ExposureClusters){
        G4bool NewExposure=true;
        //G4int ExposureID=exposure.first;
        //run over clusters in this exposure and sort by primary ID
        std::map<G4int,std::vector<NtupleForHiC_Clusters*>> PrimaryClusters;
        std::vector<std::vector<DrDamageEvent*>> PrimaryDamages;
        for (auto cluster:exposure.second){
            G4int PrimaryID=cluster.second->GetPrimaryID();
            PrimaryClusters[PrimaryID].push_back(cluster.second);
        }
        
        for (auto data:PrimaryClusters){
            std::vector<DrDamageEvent*> Damages;
            G4bool NewPrimary=true;
            for (unsigned int i=0;i<data.second.size();i++){
                
                //If not writing SSB, skip adding to PS store
                G4bool brea=false;
                while (data.second[i]->GetIsSSB() && !WriteSSB){
                    i++;
                    if (i==data.second.size()){
                        brea=true;
                        break;
                    }
                }
                if (brea){break;}
                
                DrDamageEvent*dmg=new DrDamageEvent();
                
                //Field 1: Exposure details
                //[0]: 2=New Exposure, 1=New Primary, 0=Same Primary
                //[1]: Primary ID
                G4int Field1_1=2;
                G4int Field1_2=data.second[i]->GetPrimaryID();
                if (!NewExposure){
                    if (!NewPrimary){Field1_1=0;}
                    if (NewPrimary){Field1_1=1;}
                }
                dmg->fNewEvent.push_back(Field1_1);
                dmg->fNewEvent.push_back(Field1_2);
                
                //Field 2: Position
                //[0]: Centre
                //[1]: Maximum extent
                //[2]: Minimum extent
                dmg->fPosition.push_back(data.second[i]->GetCentre());
                dmg->fPosition.push_back(data.second[i]->GetMaxPos());
                dmg->fPosition.push_back(data.second[i]->GetMinPos());
                
                
                //Field 3: Chromosome ID
                //[0]: 0=unspecified, 1=heterochromatin, 2=euchromatin, 3=free DNA, 4=viral
                //[1]: Chromosome number
                //[2]: Chromatid number
                //[3]: 0=short arm, 1=long arm
                dmg->fChromasomeID.push_back(0);
                dmg->fChromasomeID.push_back(data.second[i]->GetChromID());
                dmg->fChromasomeID.push_back(data.second[i]->GetChromatid());
                dmg->fChromasomeID.push_back(data.second[i]->GetChromArm());
                
                //Field 4: Position along chromosome
                dmg->fChromasomePosition=data.second[i]->GetChromosomePosition();
                
                //Field 5: Cause
                //[0] 0=direct, 1=indirect, 2=combo, 3=other
                //[1] n Direct
                //[2] n Indirect
                dmg->fCause=data.second[i]->GetDamageCause();
                
                //Field 6: Damage Types
                G4int nBases=data.second[i]->GetnBases();
                G4int nBacks=data.second[i]->GetnBacks();
                G4int nDSB=data.second[i]->GetIsDSB();
                dmg->fDamageTypes.push_back(nBases);
                dmg->fDamageTypes.push_back(nBacks);
                dmg->fDamageTypes.push_back(nDSB);
                
                
                
                Damages.push_back(dmg);
                NewPrimary=false;
                NewExposure=false;
                DamageEntries++;
            }
            
            PrimaryDamages.push_back(Damages);
            
        }
        
        DamageStore->push_back(PrimaryDamages);
    }
    
    return PSStore;
}

void NtupleForHiC_EndOfSession::FillPhaseSpaceHeader(std::map<G4int,G4int> ChromosomeSize,
                                                     TsDamagePhaseSpaceStore*PSStore,
                                                     TsParameterManager* pM,
                                                     G4ThreeVector NucDim,
                                                     G4int DamageEntries)
{
    
    
    //Chromosome Sizes in Mbp
    G4String ChromosomeDetails=to_string(ChromosomeSize.size())+",";
    for (auto chrom:ChromosomeSize){
        ChromosomeDetails+=to_string(chrom.second)+",";
    }
    
    //Sim details = HiC vertces input, Sensitive fraction
    G4String VerticesFile="HiC.txt";
    if (pM->ParameterExists("Ge/HiC/VerticesFile")){
        VerticesFile=pM->GetStringParameter("Ge/HiC/VerticesFile");
    }
    G4String SensFrac="0.141328949";
    if (pM->ParameterExists("Sc/HiCScorer/SensitiveFraction")){
        SensFrac=to_string(pM->GetUnitlessParameter("Sc/HiCScorer/SensitiveFraction"));
    }
    G4String SimDetails="Vertices: "+VerticesFile+", SensFrac: "+SensFrac;
    
    
    //Energy range
    G4double EMin = 5.0;
    if (pM->ParameterExists("Sc/HiCScorer/MinEneRange")){
        EMin=pM->GetDoubleParameter("Sc/HiCScorer/MinEneRange", "Energy") /eV;
    }
    G4double EMax=37.5;
    if (pM->ParameterExists("Sc/HiCScorer/MaxEneRange")){
        EMax=pM->GetDoubleParameter("Sc/HiCScorer/MaxEneRange", "Energy") /eV;
    }
    
    
    //N Primaries = histories * repeats
    G4int nPrim=
    pM->GetIntegerParameter("So/Beam/NumberOfHistoriesInRun")*
    pM->GetIntegerParameter("Tf/NumberOfSequentialTimes");
    
    //Data entries. Keep "a b c...", map sorts alphabetically
    std::map<G4String,G4int> DEnt;
    DEnt["a Classification"]=1;
    DEnt["b X,Y,Z"]=1;
    DEnt["c Chromosome IDs"]=1;
    DEnt["d Chromosome position"]=1;
    DEnt["e Cause"]=1;
    DEnt["f Damage types"]=1;
    DEnt["g Full break spec"]=0;
    DEnt["h DNA sequence"]=0;
    DEnt["i Lesion time"]=0;
    DEnt["j Particle types"]=0;
    DEnt["k Energies"]=0;
    DEnt["l Translation"]=0;
    DEnt["m Direction"]=0;
    DEnt["n Particle time"]=0;
    
    //Volumes - details fo the bounding geometry
    G4String Vol="1,";
    Vol+=to_string(NucDim[0]/um)+",";
    Vol+=to_string(NucDim[1]/um)+",";
    Vol+=to_string(NucDim[2]/um);
    
    
    
    
      
    //Fill header
    DrDamageHeader*Header=PSStore->GetDamageHeaderPtr();
    Header->fSDD_version="SDDv1.0";
    Header->fSoftware="TOPAS_Hi-C";
    Header->fAuthor="Nicholas T. Henthorn, nicholas.henthorn@manchester.ac.uk, 9/2/2021, Ingram & Henthorn et al, PLoS Comput. Biol. (2020), 10.1371/journal.pcbi.1008476";
    Header->fSimulation_details=SimDetails;
//    Header->fSource=;
//    Header->fSource_type=;
//    Header->fIncident_particles.push_back();
    Header->fMean_particle_energy.push_back(pM->GetDoubleParameter("So/Beam/BeamEnergy", "Energy"));
//    Header->fParticle_fraction.push_back();
    
    Header->fDamage_definition="0,0,3.4,-1,"+to_string(EMin)+","+to_string(EMax);
    Header->fIrradiation_target="Cell Nucleus";
    Header->fCell_cycle_phase="2,0";
    Header->fDNA_structure="0,1";
    Header->fProliferation_status="1";
    Header->fData_entries=DEnt;
    Header->fVolumes=Vol;
    
    
    
//    Header->fDNA_density=;
    Header->fIn_vitro_in_vivo=0;
//    Header->fTime=;
//    Header->fAdditional_Information=;
    
    
    
    
    Header->fChromosome_sizes=ChromosomeDetails;
    
    Header->fDamage_and_primary_count=to_string(DamageEntries)+","+to_string(nPrim);
    
    
}

void NtupleForHiC_EndOfSession::WriteSDD(G4String Filename,
                                         TsDamagePhaseSpaceStore*PSStore)
{
    //Print Header
    if (WriteSDDHeader){
        PSStore->GetDamageHeaderPtr()->PrintToFileHeader(Filename);
    }
    
    //Print Damage Entries
    for (unsigned int i=0;i<PSStore->GetDamageEventStore().size();i++){
        for (unsigned int j=0;j<PSStore->GetDamageEventStore()[i].size();j++){
            for (unsigned int k=0;k<PSStore->GetDamageEventStore()[i][j].size();k++){
                PSStore->GetDamageEventStore()[i][j][k]->PrintToFileEvent(Filename);
            }
        }
    }
}
