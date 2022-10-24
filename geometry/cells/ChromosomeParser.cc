// Extra Class for use by TsHiC
/*
*
*  Ingram & Henthorn et al
*  Hi-C TOPAS Geometry
*
*/

#include "ChromosomeParser.hh"

//read in HiC vertices file
void ChromosomeParser::ReadBeads(G4String filename){
    std::cout<<"  ---  Reading Bead Vertices"<<std::endl;
    

    std::cout<<"   --   File: "<<filename<<std::endl;
    std::ifstream infile(filename);
    G4int num=0;
    
    G4String line;
    while (getline(infile, line)){
        G4String part = line.substr(0,1);
        //ignore line starting with #
        if(part != "#"){
            G4String field;
            std::stringstream stream(line.data());
            G4String Chrom;
            G4int chromosome=-1, genomeLength=0;
            G4double x=0.0,y=0.0,z=0.0,radius=0.0;
            G4String copy;
            stream
            >>Chrom
            >>copy
            >>x
            >>y
            >>z
            >>radius
            >>genomeLength;

            num++;
                    
            //renumber x=23, y=46
            if (Chrom=="chrX"){Chrom="chr23";}
            if (Chrom=="chrY"){Chrom="chr46";}
            
            //convert string to right of "r" in "chr" to int
            chromosome=stoi(Chrom.substr(Chrom.find("r") + 1));
                    
            if (copy=="B"){chromosome+=23;}
                    
                    
            //sum length of chromosome
            if (ChromSize.find(chromosome)!=ChromSize.end()){
                    ChromSize[chromosome]+=genomeLength;
                }
                else{
                    ChromSize[chromosome]=genomeLength;
            }


            ChromObj obj;
            obj.ChromID=chromosome;
            obj.Pos=G4ThreeVector(x*um,y*um,z*um);
            obj.Radius=radius*um;
            obj.GenomeLength=genomeLength;
            Beads[chromosome].push_back(obj);
        }
        
    }
    
    //work out fraction of chromosome covered by each bead
    //and set the ID of the volume
    for (auto data:Beads){
        G4double MinFract=0.0;
        G4double MaxFract = 0.0;
        G4int MinBP=1;
        G4int MaxBP=0;
        
        for (unsigned int i=0;i<data.second.size();i++){
            G4double FractOfChrom = data.second[i].GenomeLength/(G4double)ChromSize[data.first];
            G4int CopyNum = i;
            MaxFract=MinFract+FractOfChrom;
            if (MaxFract>=1.0){MaxFract=0.999999999;} //max fract along cannot = 1
            MaxBP=MinBP+data.second[i].GenomeLength;
            Beads[data.first][i].FractionOfChrom=FractOfChrom;
            Beads[data.first][i].VolumeCopyNum=CopyNum;
            Beads[data.first][i].MinFractAlong=MinFract;
            Beads[data.first][i].MaxFractAlong=MaxFract;
            Beads[data.first][i].MinBPAlong=MinBP;
            Beads[data.first][i].MaxBPAlong=MaxBP;
            MinFract=MaxFract;
            MinBP=MaxBP+1;
        }
    }
    
    infile.close();
}

//This function recentres the HiC polymer beads about 0,0,0
void ChromosomeParser::RecentreBeads(){
    std::cout<<"  ---  ReCentring Beads"<<std::endl;
    G4double ComX=0.0,ComY=0.0, ComZ=0.0;
    G4double RadialCentre=0.0;
    G4int nBeads=0;
    
    //Find Max Beads
    G4double MaxX=std::numeric_limits<G4double>::min(), MaxY=std::numeric_limits<G4double>::min(), MaxZ=std::numeric_limits<G4double>::min();
    G4double MinX=std::numeric_limits<G4double>::max(), MinY=std::numeric_limits<G4double>::max(), MinZ=std::numeric_limits<G4double>::max();
    for (auto data:Beads){
        for (unsigned int i=0;i<data.second.size();i++){
            G4double radius=data.second[i].Radius;
            G4double CurrX=data.second[i].Pos.x();
            G4double CurrY=data.second[i].Pos.y();
            G4double CurrZ=data.second[i].Pos.z();
            if (CurrX+radius>MaxX){MaxX=CurrX+radius;}
            if (CurrY+radius>MaxY){MaxY=CurrY+radius;}
            if (CurrZ+radius>MaxZ){MaxZ=CurrZ+radius;}
            if (CurrX-radius<MinX){MinX=CurrX-radius;}
            if (CurrY-radius<MinY){MinY=CurrY-radius;}
            if (CurrZ-radius<MinZ){MinZ=CurrZ-radius;}
            
            RadialCentre+=data.second[i].Pos.mag()+radius;
            ComX+=CurrX;
            ComY+=CurrY;
            ComZ+=CurrZ;
            nBeads++;
        }
    }
    
    RadialCentre/=(G4double)nBeads;
    ComX/=(G4double)nBeads;
    ComY/=(G4double)nBeads;
    ComZ/=(G4double)nBeads;
    std::cout<<"   --   COM_Origi: "<<ComX/um<<","<<ComY/um<<","<<ComZ/um<<" um, Radial: "<<RadialCentre/um<<" um"<<std::endl;
    RadialCentre=ComX=ComY=ComZ=0.0;
    
    //central axis coordinate = max - diff/2
    G4double CentreX=MaxX-((MaxX-MinX)/2.0);
    G4double CentreY=MaxY-((MaxY-MinY)/2.0);
    G4double CentreZ=MaxZ-((MaxZ-MinZ)/2.0);
    G4ThreeVector Centre(CentreX,CentreY,CentreZ);
    
    //Translate beads in XYZ by pos-centre
    for (auto data:Beads){
        for (unsigned int i=0;i<data.second.size();i++){
            data.second[i].Pos=(data.second[i].Pos-Centre);
            G4double CurrX=data.second[i].Pos.x();
            G4double CurrY=data.second[i].Pos.y();
            G4double CurrZ=data.second[i].Pos.z();
            
            RadialCentre+=data.second[i].Pos.mag()+data.second[i].Radius;
            ComX+=CurrX;
            ComY+=CurrY;
            ComZ+=CurrZ;
        }
    }
    
    RadialCentre/=(G4double)nBeads;
    ComX/=(G4double)nBeads;
    ComY/=(G4double)nBeads;
    ComZ/=(G4double)nBeads;
    std::cout<<"   --   COM_Trans: "<<ComX/um<<","<<ComY/um<<","<<ComZ/um<<" um, Radial: "<<RadialCentre/um<<" um"<<std::endl;
}

//calculate the nucleus size to be most radial bead
void ChromosomeParser::ResizeNucleus(){
    std::cout<<"  ---  ReSizing Nucleus"<<std::endl;
    G4double MaxX=0.0;
    G4double MaxY=0.0;
    G4double MaxZ=0.0;
    
    G4double PosExtX=0.0;
    G4double PosExtY=0.0;
    G4double PosExtZ=0.0;
    G4double NegExtX=std::numeric_limits<G4double>::max();
    G4double NegExtY=std::numeric_limits<G4double>::max();
    double NegExtZ=std::numeric_limits<G4double>::max();
    for (auto data:Beads){
        for (unsigned int i=0;i<data.second.size();i++){
            G4double x=data.second[i].Pos.x();
            G4double y=data.second[i].Pos.y();
            G4double z=data.second[i].Pos.z();
            //calculate max radial
            if (MaxX<fabs(x)+data.second[i].Radius){MaxX=x+data.second[i].Radius;}
            if (MaxY<fabs(y)+data.second[i].Radius){MaxY=y+data.second[i].Radius;}
            if (MaxZ<fabs(z)+data.second[i].Radius){MaxZ=z+data.second[i].Radius;}
            
            //calculate max and min xyz (get length in dimensions)
            //later rotate so small axis perp to z
            if (PosExtX<x+data.second[i].Radius){PosExtX=x+data.second[i].Radius;}
            if (PosExtY<y+data.second[i].Radius){PosExtY=y+data.second[i].Radius;}
            if (PosExtZ<z+data.second[i].Radius){PosExtZ=z+data.second[i].Radius;}
            if (NegExtX>x-data.second[i].Radius){NegExtX=x-data.second[i].Radius;}
            if (NegExtY>y-data.second[i].Radius){NegExtY=y-data.second[i].Radius;}
            if (NegExtZ>z-data.second[i].Radius){NegExtZ=z-data.second[i].Radius;}
            
        }
    }
    
    Dimensions = G4ThreeVector(MaxX+2.0*um,MaxY+2.0*um,MaxZ+2.0*um);
    std::cout<<"   --   Nucleus Dimensions (um): "<<Dimensions/um<<std::endl;
    
    G4double xLength=PosExtX-NegExtX;
    G4double yLength=PosExtY-NegExtY;
    G4double zLength=PosExtZ-NegExtZ;
    
    //determine short axis and rotate it perp to z
    rot=new G4RotationMatrix();
    if (xLength<yLength &&
        xLength<zLength){
        rot->rotateY(90.0*deg);
        std::cout<<"   --   Rotate Cytoplasm About Y"<<std::endl;
    }
    if (yLength<xLength &&
        yLength<zLength){
        rot->rotateX(90.0*deg);
        std::cout<<"   --   Rotate Cytoplasm About X"<<std::endl;
    }
}
