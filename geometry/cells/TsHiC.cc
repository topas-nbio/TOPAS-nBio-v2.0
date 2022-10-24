// Component for TsHiC
/*
*
*  Ingram & Henthorn et al
*  Hi-C TOPAS Geometry
*
*/

#include "TsHiC.hh"



TsHiC::TsHiC(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
                 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
}


TsHiC::~TsHiC()
{
}



G4VPhysicalVolume* TsHiC::Construct()
{
    BeginConstruction();
    
    //Read HiC vertices
    ChromosomeParser* chrom = new ChromosomeParser();
    G4String filename="HiC.txt";
    if (fPm->ParameterExists(GetFullParmName("VerticesFile"))){
        filename=fPm->GetStringParameter(GetFullParmName("VerticesFile"));
    
    }
    chrom->ReadBeads(filename);
    chrom->RecentreBeads();
    chrom->ResizeNucleus();
    
    std::map<G4int,std::vector<ChromObj>> Beads = chrom->GetBeads();
    G4ThreeVector Dimensions = chrom->GetDimensions();
    G4RotationMatrix* rot = chrom->GetRotation();
    delete chrom;
    
    
    //Set some global parameters
    ResolveParameters();

    
    //create bounding box
    G4Box* gCell = new G4Box(fName,
                             ExtraCytoplasmSize+transX+Dimensions.x()+3.0*um,
                             ExtraCytoplasmSize+transY+Dimensions.y()+3.0*um,
                             ExtraCytoplasmSize+transZ+Dimensions.z()+3.0*um);
    fEnvelopeLog = CreateLogicalVolume(gCell);
    fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);
    
    //Place the cytoplasm, rotate short dimension parrallel to z-axis
    G4VPhysicalVolume*pCytoplasm=BuildCytoplasm(Dimensions,
                                                ExtraCytoplasmSize,
                                                rot);
    
    //Place the nucleus, daughter of cytoplasm
    G4VPhysicalVolume*pNucleus=BuildNucleus(Dimensions,
                                            pCytoplasm);
    
    //Place the Hi-C polymer beads
    BuildChromosomes(pNucleus,
                     Beads);

    


    InstantiateChildren(fEnvelopePhys);


    return fEnvelopePhys;
}

//some global parameters
void TsHiC::ResolveParameters()
{
    //cell position
    G4String name1 = GetFullParmName("Cytoplasm/transX");
    if (fPm -> ParameterExists(name1)){
        transX = fPm->GetDoubleParameter(name1, "Length");
    }
    name1 = GetFullParmName("Cytoplasm/transY");
    if (fPm -> ParameterExists(name1)){
        transY = fPm->GetDoubleParameter(name1, "Length");
    }
    name1 = GetFullParmName("Cytoplasm/transZ");
    if (fPm -> ParameterExists(name1)){
        transZ = fPm->GetDoubleParameter(name1, "Length");
    }
    
    //exra size of cytoplasm around nucleus
    name1 = GetFullParmName("Cytoplasm/ExtraSpacing");
    if (fPm -> ParameterExists(name1)){
        ExtraCytoplasmSize = fPm->GetDoubleParameter(name1, "Length");
    }
    
    //chromosome colours
    VolumeColor[1]=G4Color(0.49247182453900407,0.0659001922311141,0.012988584024996808);
    VolumeColor[2]=G4Color(0.23140307968667773,0.8033827822443745,0.6393786641060017);
    VolumeColor[3]=G4Color(0.7585813217010804,0.9789711950363517,0.04221926337005577);
    VolumeColor[4]=G4Color(0.6673318619120212,0.9091785613101261,0.40749183413727263);
    VolumeColor[5]=G4Color(0.5438324400450564,0.8582016971322097,0.7298291678491887);
    VolumeColor[6]=G4Color(0.5141843924132756,0.4334353128660875,0.06194062319464788);
    VolumeColor[7]=G4Color(0.06864444688357363,0.14248557431722153,0.11511212481096811);
    VolumeColor[8]=G4Color(0.6933887874880108,0.01291826248331085,0.2549173411321378);
    VolumeColor[9]=G4Color(0.2251048540298951,0.28361025068038936,0.6255927899481798);
    VolumeColor[10]=G4Color(0.3977974846409681,0.6391693261056358,0.8309788808771028);
    VolumeColor[11]=G4Color(0.8059894651576868,0.5593126506944623,0.7954053403801448);
    VolumeColor[12]=G4Color(0.8408359866041606,0.8357053359530919,0.9483472082473238);
    VolumeColor[13]=G4Color(0.462512080629906,0.3085490759997732,0.6424884993194124);
    VolumeColor[14]=G4Color(0.7700282930625697,0.30796019975503697,0.05497354567683577);
    VolumeColor[15]=G4Color(0.6034410328438998,0.8370792540551809,0.3595951478515589);
    VolumeColor[16]=G4Color(0.7671376524276918,0.477774350838304,0.563476858200768);
    VolumeColor[17]=G4Color(0.9991196735013264,0.834617737773388,0.7727825764337234);
    VolumeColor[18]=G4Color(0.11998640724326837,0.531481038067903,0.03724809448701927);
    VolumeColor[19]=G4Color(0.514670247660917,0.879570429037138,0.5673084127632296);
    VolumeColor[20]=G4Color(0.530813288817109,0.7010723586948585,0.6685790711907813);
    VolumeColor[21]=G4Color(0.06211092281215558,0.7248169987691381,0.7049864602754835);
    VolumeColor[22]=G4Color(0.6126041271212809,0.7236663637285328,0.09573925444515718);
    VolumeColor[23]=G4Color(0.9572274585786695,0.9440878206384548,0.10094721430277032);
    VolumeColor[24]=G4Color(0.8576378536672312,0.6891210878060025,0.44688146558955866);
    VolumeColor[25]=G4Color(0.16186046181104918,0.2250166427782745,0.8523452328220089);
    VolumeColor[26]=G4Color(0.7703513648078352,0.3955846626017572,0.9480492827724442);
    VolumeColor[27]=G4Color(0.07660658776310814,0.4013660573615594,0.518513001579725);
    VolumeColor[28]=G4Color(0.35834326858230703,0.9536623040045678,0.745731737953009);
    VolumeColor[29]=G4Color(0.9055728481985799,0.03847836921902892,0.6707334874695682);
    VolumeColor[30]=G4Color(0.611247021751855,0.3810769724147606,0.03166671787052788);
    VolumeColor[31]=G4Color(0.5739452460644697,0.06256064849253895,0.4540457516000722);
    VolumeColor[32]=G4Color(0.6879484714161074,0.35224100875483144,0.1395820584440488);
    VolumeColor[33]=G4Color(0.23667177191333533,0.7202888386603635,0.32178340212327583);
    VolumeColor[34]=G4Color(0.6283905266814671,0.5193564822241657,0.2076971939335739);
    VolumeColor[35]=G4Color(0.5498084969736927,0.4075550715325341,0.3931864908374545);
    VolumeColor[36]=G4Color(0.5674445955039258,0.9665854807167373,0.7368834782260029);
    VolumeColor[37]=G4Color(0.7285140452838902,0.5464894067255806,0.7433571744920119);
    VolumeColor[38]=G4Color(0.1729566031630141,0.46342590975582987,0.3050829342703273);
    VolumeColor[39]=G4Color(0.11937788820013084,0.4343837228764069,0.5228209121396812);
    VolumeColor[40]=G4Color(0.21766414980751214,0.40219778022020103,0.3779886809373476);
    VolumeColor[41]=G4Color(0.5669540124932198,0.1977741541363588,0.008850956336362237);
    VolumeColor[42]=G4Color(0.6222613914635593,0.5943619587851444,0.7770025654255228);
    VolumeColor[43]=G4Color(0.3223062840770442,0.3035451870535274,0.06416675955458895);
    VolumeColor[44]=G4Color(0.7238608202842554,0.38051455892379116,0.36155794257635676);
    VolumeColor[45]=G4Color(0.8739029288734816,0.9517004867511928,0.13142495401866805);
    VolumeColor[46]=G4Color(0.5998712010463205,0.0173406650149196,0.856335111668236);
}



//place the cytoplasm volume
//has extra size compared to nucleus
//rotated so that small dimension is parallel to z-axis
G4VPhysicalVolume*TsHiC::BuildCytoplasm(G4ThreeVector Dimensions,
                                        G4double ExtraCytoplasm,
                                        G4RotationMatrix*rot)
{
    std::cout<<"  ---  Building Cytoplasm Volume"<<std::endl;
    
    G4String subComponentName1 = "Cytoplasm";

    G4ThreeVector*cytoplasmPos = new G4ThreeVector(transX,transY,transZ);
    
    G4Ellipsoid * gCytoplasm = new G4Ellipsoid("gCytoplasm",
                                                Dimensions.x()+ExtraCytoplasm,
                                                Dimensions.y()+ExtraCytoplasm,
                                                Dimensions.z()+ExtraCytoplasm);
    
    G4LogicalVolume* lCytoplasm = CreateLogicalVolume(subComponentName1, gCytoplasm);
    G4VPhysicalVolume* pCytoplasm = CreatePhysicalVolume(subComponentName1, lCytoplasm, rot, cytoplasmPos, fEnvelopePhys);

    std::cout<<"   --   Cytoplasm Dimensions (preRot): ("
    <<(Dimensions.x()+ExtraCytoplasm)/um<<","
    <<(Dimensions.y()+ExtraCytoplasm)/um<<","
    <<(Dimensions.z()+ExtraCytoplasm)/um<<") um"
    <<std::endl;
    
    return pCytoplasm;
}

G4VPhysicalVolume*TsHiC::BuildNucleus(G4ThreeVector Dimensions,
                                      G4VPhysicalVolume* pMother)
{
    std::cout<<"  ---  Building Nucleus Volume"<<std::endl;
    
    G4String subComponentName2 = "Nucleus";
    
    
    G4ThreeVector*nucleusPos = new G4ThreeVector(0,0,0);
    G4Ellipsoid * gNucleus = new G4Ellipsoid("gNucleus",
                                             Dimensions.x(),
                                             Dimensions.y(),
                                             Dimensions.z());
    
    
    G4LogicalVolume* lNucleus = CreateLogicalVolume(subComponentName2, gNucleus);
    G4VPhysicalVolume* pNucleus = CreatePhysicalVolume(subComponentName2, lNucleus, 0, nucleusPos, pMother);
    
    return pNucleus;
}


void TsHiC::BuildChromosomes(G4VPhysicalVolume* pMother,
                             std::map<G4int,std::vector<ChromObj>> Beads)
{
    std::cout<<"  ---  Building Chromosome Volumes"<<std::endl;
    
    G4String Material = "G4_WATER";
    G4int TotalBeads=0;
    for (auto chrom:Beads){
        //TOPAS appends copy number to end of name, use "-" to identify
        G4String SubCompName="Chromosome_"+std::to_string(chrom.first)+"-";
        G4String test=SubCompName.substr(0, SubCompName.size()-1);
        std::cout<<"   --   "<<test<<", "<<chrom.second.size()<<" beads"<<std::endl;
        TotalBeads+=chrom.second.size();
        G4VisAttributes*Vis=new G4VisAttributes(VolumeColor[chrom.first]);
        Vis->SetForceSolid(true);
        RegisterVisAtt(Vis);
        
        for (unsigned int i=0; i<chrom.second.size(); i++){
            G4ThreeVector* pos = new G4ThreeVector(chrom.second[i].Pos);
            G4Orb* gChromosome = new G4Orb("gChromosome", chrom.second[i].Radius);
            G4LogicalVolume* lChromosome = CreateLogicalVolume(SubCompName, Material, gChromosome);
            CreatePhysicalVolume(SubCompName, chrom.second[i].VolumeCopyNum, false, lChromosome, 0, pos, pMother);
            lChromosome->SetVisAttributes(Vis);
        }
    }

    
    std::cout<<"    -    Total beads = "<<TotalBeads<<std::endl;
    
}





