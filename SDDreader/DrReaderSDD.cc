// Extra Class for use by TOPAS
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
// Created by John-William Warmenhoven on 15/10/18.
// Please cite: https://doi.org/10.1667/RR15209.1
//

#include "DrReaderSDD.hh"
#include <regex>
#include <fstream>

using namespace std;

DrReaderSDD::DrReaderSDD() {}

DrReaderSDD::~DrReaderSDD() {}

void DrReaderSDD::ReadSTDInput(G4String fileName, DrDamageHeader* ptrHeader,vEvents& eventList){
    ReadHeader(fileName, ptrHeader);
    ReadDamage(fileName, ptrHeader, eventList);
}

void DrReaderSDD::ReadHeader(G4String fileName, DrDamageHeader* header) {
    ifstream inFile(fileName);
    if (!inFile) {
        G4cerr << "ERROR: could not read damage data, file not found." << G4endl;
        exit(EXIT_FAILURE);
    }

    G4String com = "(?:\\s*\\,\\s*)";
    G4String col = "(?:\\s*\\;\\s*)";
    G4String any = "(?:(?:\\s*\\,\\s*)([^;]*))?";
    G4String post = "(?:\\s*;)";
    G4String rint = "(-*\\d+)";
    G4String rPInt = "(\\d+)";
    G4String rdbl = "(-*\\d+\\.*\\d*)";
    G4String rBool = "([0,1])";
    G4String c3vec = rdbl+"?"+com+rdbl+"?"+com+rdbl+"?";

    //@@@@ Grabbing free text
    regex version("(sdd\\sversion)"+any+post, regex_constants::icase);
    regex software("(software)"+any+post, regex_constants::icase);
    regex author("(author)"+any+post, regex_constants::icase);
    regex simDetails("(simulation\\sdetails)"+any+post, regex_constants::icase);
    regex source("(source)"+any+post, regex_constants::icase);
    regex sourceType("(source\\stype)"+any+post, regex_constants::icase);
    regex energyDist("(energy\\sdistribution)"+any+post, regex_constants::icase);

    regex doseFluence("(dose\\sor\\sfluence)"+any+post, regex_constants::icase);
    regex doseRate("(dose\\srate)"+any+post, regex_constants::icase);
    regex irrTarget("(irradiation\\starget)"+any+post, regex_constants::icase);
    regex chroSizes("(chromosome\\ssizes)"+any+post, regex_constants::icase);
    regex additional("(additional\\sinformation)"+any+post, regex_constants::icase);

    //@@@@ Grabbing specific structures
    regex incidentParticles("(incident\\sparticles)(((?:"+com+rPInt+")+))?"+post, regex_constants::icase);
    regex meanPartEnergy("(mean\\sparticle\\senergy)(((?:"+com+rdbl+")+))?"+post, regex_constants::icase);
    regex partFract("(particle\\sfraction)(((?:"+com+rdbl+")+))?"+post, regex_constants::icase);
    regex volumes("(volumes)(?:"+com+"("+rint+com+c3vec+com+c3vec+com+rint+com+c3vec+com+c3vec+"))?"+post, regex_constants::icase);
    regex dnaDensity("(dna\\sdensity)(?:"+com+rdbl+")?"+post, regex_constants::icase);
    regex cellCycle("(cell\\scycle\\sphase)(?:"+com+"("+rint+com+"?"+rdbl+"?))?"+post, regex_constants::icase);
    regex dnaStruct("(dna\\sstructure)(?:"+com+"("+rint+com+"?"+rint+"?))?"+post, regex_constants::icase);
    regex vitroVivo("(in\\svitro\\s\\/\\sin\\svivo)(?:"+com+rPInt+")?"+post, regex_constants::icase);
    regex proliferation("(proliferation\\sstatus)(?:"+com+"("+rint+any+"))?"+post, regex_constants::icase);
    regex microEnv("(microenvironment)(?:"+com+"("+rdbl+"?"+com+"?"+rdbl+"?))?"+post, regex_constants::icase);
    regex damageDef("(damage\\sdefinition)(?:"+com+"("+rint+"?"+com+rBool+"?"+com+rdbl+"?"+com+rdbl+"?"+com+rdbl+"?"+com+rdbl+"?))?"+post, regex_constants::icase);
    regex simTime("(time)(?:"+com+rdbl+")?"+post, regex_constants::icase);
    regex damageCount("(damage\\sand\\sprimary\\scount)(?:"+com+"("+rint+"?"+com+rint+"?))?"+post, regex_constants::icase);
    regex dataEntries("(data\\sentries)"+com+rBool+com+rBool+com+rBool+com+rBool+com+rBool+com+rBool+com+rBool+com+rBool+com+rBool+com+rBool+com+rBool+com+rBool+com+rBool+com+rBool+post, regex_constants::icase);

    //@@@@ Grabbing end of header
    regex endHead("(\\*\\*\\*endofheader\\*\\*\\*)"+post, regex_constants::icase);
    G4bool endHeadB{false};

    G4String line;
    smatch matchName;

    while (getline(inFile, line)) {
        //@@@@ If end of header has yet to be found we are still processing header fields
        //@@@@ If end of header has been found no longer need to regex_search for header
        //@@@@ fields so skip to save time.
        if(!endHeadB){
            //@@@@ Simple string fields [0] is whole match [1] is name of field and [2]
            //@@@@ is the string found for that field.
            if(regex_search(line,matchName,version)){
                header->fSDD_version = matchName[2];
            }
            if(regex_search(line,matchName,software)){
                header->fSoftware = matchName[2];
            }
            if(regex_search(line,matchName,author)){
                header->fAuthor = matchName[2];
            }
            if(regex_search(line,matchName,simDetails)){
                header->fSimulation_details = matchName[2];
            }
            if(regex_search(line,matchName,source)){
                header->fSource = matchName[2];
            }

            if(regex_search(line,matchName,sourceType)){
                header->fSource_type = std::stoi(matchName[2].str());
            }
            if(regex_search(line,matchName,incidentParticles)){
                G4String tempS = matchName[2].str();
                regex csInt("(?:\\s*\\,\\s*)(\\d+)");
                sregex_iterator next(tempS.begin(),tempS.end(),csInt);
                sregex_iterator end;
                while(next != end){
                    smatch match = *next;
                    header->fIncident_particles.push_back(std::stoi(match[1].str()));
                    next++;
                }
            }
            if(regex_search(line,matchName,meanPartEnergy)){
                G4String tempS = matchName[2].str();
                regex csInt("(?:\\s*\\,\\s*)(\\d+\\.*\\d*)");
                sregex_iterator next(tempS.begin(),tempS.end(),csInt);
                sregex_iterator end;
                while(next != end){
                    smatch match = *next;
                    header->fMean_particle_energy.push_back(std::stod(match[1].str()));
                    next++;
                }
            }
            if(regex_search(line,matchName,energyDist)){
                header->fEnergy_distribution = matchName[2];
            }
            if(regex_search(line,matchName,partFract)){
                G4String tempS = matchName[2].str();
                regex csInt("(?:\\s*\\,\\s*)(\\d+\\.*\\d*)");
                sregex_iterator next(tempS.begin(),tempS.end(),csInt);
                sregex_iterator end;
                while(next != end){
                    smatch match = *next;
                    header->fParticle_fraction.push_back(std::stod(match[1].str()));
                    next++;
                }
            }
            if(regex_search(line,matchName,doseFluence)){
                header->fDose_or_fluence = matchName[2];
            }
            if(regex_search(line,matchName,doseRate)){
                header->fDose_rate = matchName[2];
            }
            if(regex_search(line,matchName,irrTarget)){
                header->fIrradiation_target = matchName[2];
            }
            if(regex_search(line,matchName,chroSizes)){
                header->fChromosome_sizes = matchName[2];
            }
            if(regex_search(line,matchName,additional)){
                header->fAdditional_Information = matchName[2];
            }
            if(regex_search(line,matchName,volumes)){
                header->fVolumes = matchName[2];
            }
            if(regex_search(line,matchName,dnaDensity)){
                header->fDNA_density = stod(matchName[2].str());
            }
            if(regex_search(line,matchName,cellCycle)){
                header->fCell_cycle_phase = matchName[2];
            }
            if(regex_search(line,matchName,dnaStruct)){
                header->fDNA_structure = matchName[2];
            }
            if(regex_search(line,matchName,vitroVivo)){
                header->fIn_vitro_in_vivo = stoi(matchName[2].str());
            }
            if(regex_search(line,matchName,proliferation)){
                header->fProliferation_status = matchName[2];
            }
            if(regex_search(line,matchName,microEnv)){
                header->fMicroenvironment = matchName[2];
            }
            if(regex_search(line,matchName,damageDef)){
                header->fDamage_definition = matchName[2];
            }
            if(regex_search(line,matchName,simTime)){
                header->fTime = stod(matchName[2].str());
            }
            if(regex_search(line,matchName,damageCount)){
                header->fDamage_and_primary_count = matchName[2];
            }
            if(regex_search(line,matchName,dataEntries)){
                header->fData_entries["Classification"] = stoi(matchName[2].str());
                header->fData_entries["Position"] = stoi(matchName[3].str());
                header->fData_entries["ChromosomeID"] = stoi(matchName[4].str());
                header->fData_entries["ChromosomePos"] = stoi(matchName[5].str());
                header->fData_entries["Cause"] = stoi(matchName[6].str());
                header->fData_entries["DamageTypes"] = stoi(matchName[7].str());
                header->fData_entries["FullBreakSpec"] = stoi(matchName[8].str());
                header->fData_entries["DNASequence"] = stoi(matchName[9].str());
                header->fData_entries["LesionTime"] = stoi(matchName[10].str());
                header->fData_entries["ParticleTypes"] = stoi(matchName[11].str());
                header->fData_entries["Energies"] = stoi(matchName[12].str());
                header->fData_entries["Translation"] = stoi(matchName[13].str());
                header->fData_entries["Direction"] = stoi(matchName[14].str());
                header->fData_entries["ParticleTime"] = stoi(matchName[15].str());
            }
            if(regex_search(line,matchName,endHead)){
                G4cout<<"Found end of header: "<<matchName[0]<<G4endl;
                endHeadB=true;
            }
        }
    }
}

void DrReaderSDD::ReadDamage(G4String fileName, DrDamageHeader* header, vEvents& eventList){

    ifstream inFile(fileName);
    if (!inFile) {
        G4cerr << "ERROR: could not read damage data, file not found." << G4endl;
        exit(EXIT_FAILURE);
    }

    G4String col = "(?:\\s*\\;\\s*)";
    G4String com = "\\s*\\,\\s*";
    G4String sla = "\\s*\\/\\s*";
    G4String Rint = "-?\\d+";
    G4String pint = "\\d+";
    G4String dbl = "-?\\d\\.?\\d*(?:e\\+-?\\d+)?";
    G4String R3v_c = dbl+com+dbl+com+dbl;
    G4String R3v_s = dbl+sla+dbl+sla+dbl;

    //@@@@ Damage fields
    map<G4String,G4String> lookFor;
    lookFor["Classification"] = "("+pint+"(?:"+com+Rint+")?)"+col;
    lookFor["Position"] = "("+R3v_c+"(?:"+sla+R3v_c+"){2})"+col;
    lookFor["ChromosomeID"] = "("+pint+com+pint+com+pint+com+pint+")"+col;
    lookFor["ChromosomePos"] = "("+dbl+")"+col;
    lookFor["Cause"] = "("+pint+"(?:"+com+pint+com+pint+")?)"+col;
    lookFor["DamageTypes"] = "("+pint+com+pint+com+pint+")"+col;
    lookFor["FullBreakSpec"] = "((?:"+pint+com+pint+com+pint+sla+")+)"+col;
    lookFor["DNASequence"] = "((?:"+pint+com+pint+com+pint+sla+")+)"+col;
    lookFor["LesionTime"] = "("+dbl+"(?:"+sla+dbl+")*)"+col;
    lookFor["ParticleTypes"] = "("+pint+"(?:"+com+pint+")*)"+col;
    lookFor["Energies"] = "("+dbl+"(?:"+com+dbl+")*)"+col;
    lookFor["Translation"] = "("+R3v_s+"(?:"+com+R3v_s+")*)"+col;
    lookFor["Direction"] = "("+R3v_s+"(?:"+com+R3v_s+")*)"+col;
    lookFor["ParticleTime"] = "("+dbl+"(?:"+com+dbl+")*)"+col;

    G4String newDamageS = "";

    if(header->fData_entries["Classification"] == 1){ newDamageS += lookFor["Classification"];}
    if(header->fData_entries["Position"] == 1){ newDamageS += lookFor["Position"];}
    if(header->fData_entries["ChromosomeID"] == 1){ newDamageS += lookFor["ChromosomeID"];}
    if(header->fData_entries["ChromosomePos"] == 1){ newDamageS += lookFor["ChromosomePos"];}
    if(header->fData_entries["Cause"] == 1){ newDamageS += lookFor["Cause"];}
    if(header->fData_entries["DamageTypes"] == 1){ newDamageS += lookFor["DamageTypes"];}
    if(header->fData_entries["FullBreakSpec"] == 1){ newDamageS += lookFor["FullBreakSpec"];}
    if(header->fData_entries["DNASequence"] == 1){ newDamageS += lookFor["DNASequence"];}
    if(header->fData_entries["LesionTime"] == 1){ newDamageS += lookFor["LesionTime"];}
    if(header->fData_entries["ParticleTypes"] == 1){ newDamageS += lookFor["ParticleTypes"];}
    if(header->fData_entries["Energies"] == 1){ newDamageS += lookFor["Energies"];}
    if(header->fData_entries["Translation"] == 1){ newDamageS += lookFor["Translation"];}
    if(header->fData_entries["Direction"] == 1){ newDamageS += lookFor["Direction"];}
    if(header->fData_entries["ParticleTime"] == 1){ newDamageS += lookFor["ParticleTime"];}

    regex newDamage(newDamageS);
    G4String line;
    smatch matchName;

    vector<vector<DrDamageEvent*> > vvStore;
    vector<DrDamageEvent*> vStore;

    while (getline(inFile, line)) {
        sregex_iterator nextNewDamage(line.begin(),line.end(),newDamage);
        sregex_iterator endNewDamage;
        while(nextNewDamage != endNewDamage){
            matchName = *nextNewDamage;

            //@@@@ matchName[0] will always be the full match so
            //@@@@ we start "i" at 1.
            G4int i{1};
            auto* tempEvent = new DrDamageEvent();
            if(header->fData_entries["Classification"] == 1){
                G4String tempS = matchName[i].str();
                regex rgxItt(Rint);
                sregex_iterator next(tempS.begin(),tempS.end(),rgxItt);
                sregex_iterator end;
                while (next != end){
                    smatch match = *next;
                    tempEvent->fNewEvent.push_back(stoi(match[0].str()));
                    next++;
                }
                i++;
            }
            if(header->fData_entries["Position"] == 1){
                G4String tempS = matchName[i].str();
                regex rgxItt("("+dbl+")(?:"+com+")("+dbl+")(?:"+com+")("+dbl+")");
                sregex_iterator next(tempS.begin(),tempS.end(),rgxItt);
                sregex_iterator end;
                while(next != end){
                    smatch match = *next;
                    tempEvent->fPosition.emplace_back(stod(match[1].str()),stod(match[2].str()),stod(match[3].str()));
                    next++;
                }
                i++;
            }
            if(header->fData_entries["ChromosomeID"] == 1){
                G4String tempS = matchName[i].str();
                regex rgx("("+pint+")(?:"+com+")("+pint+")(?:"+com+")("+pint+")(?:"+com+")("+pint+")");
                smatch match;
                regex_match(tempS,match,rgx);
                tempEvent->fChromasomeID = {stoi(match[1].str()),stoi(match[2].str()),stoi(match[3].str()),stoi(match[4].str())};
                i++;
            }
            if(header->fData_entries["ChromosomePos"] == 1){
                tempEvent->fChromasomePosition = stod(matchName[i].str());
                i++;
            }
            if(header->fData_entries["Cause"] == 1){
                G4String tempS = matchName[i].str();
                regex rgxItt(pint);
                sregex_iterator next(tempS.begin(),tempS.end(),rgxItt);
                sregex_iterator end;
                vector<G4double> tempVec;
                while(next != end){
                    smatch match = *next;
                    tempEvent->fCause.push_back(stoi(match[0].str()));
                    next++;
                }
                i++;
            }
            if(header->fData_entries["DamageTypes"] == 1){
                G4String tempS = matchName[i].str();
                regex rgx("("+pint+")(?:"+com+")("+pint+")(?:"+com+")("+pint+")");
                smatch match;
                regex_match(tempS,match,rgx);
                tempEvent->fDamageTypes = {stoi(match[1].str()),stoi(match[2].str()),stoi(match[3].str())};
                i++;
            }
            if(header->fData_entries["FullBreakSpec"] == 1){
                G4String tempS = matchName[i].str();
                regex rgxItt("("+pint+")(?:"+com+")("+pint+")(?:"+com+")("+pint+")(?:"+sla+")");
                sregex_iterator next(tempS.begin(),tempS.end(),rgxItt);
                sregex_iterator end;
                while(next != end){
                    smatch match = *next;
                    tempEvent->fFullBreakStructure.push_back({stoi(match[1].str()),stoi(match[2].str()),stoi(match[3].str())});
                    next++;
                }
                i++;
            }
            if(header->fData_entries["DNASequence"] == 1){
                G4String tempS = matchName[i].str();
                regex rgxItt("("+pint+")(?:"+com+")("+pint+")(?:"+com+")("+pint+")(?:"+sla+")");
                sregex_iterator next(tempS.begin(),tempS.end(),rgxItt);
                sregex_iterator end;
                while(next != end){
                    smatch match = *next;
                    tempEvent->fDNASequence.push_back({stoi(match[1].str()),stoi(match[2].str()),stoi(match[3].str())});
                    next++;
                }
                i++;
            }
            if(header->fData_entries["LesionTime"] == 1){
                G4String tempS = matchName[i].str();
                regex rgxItt(dbl);
                sregex_iterator next(tempS.begin(),tempS.end(),rgxItt);
                sregex_iterator end;
                while(next != end){
                    smatch match = *next;
                    tempEvent->fLesionTime.push_back(stod(match[0].str()));
                    next++;
                }
                i++;
            }
            if(header->fData_entries["ParticleTypes"] == 1){
                G4String tempS = matchName[i].str();
                regex rgxItt(pint);
                sregex_iterator next(tempS.begin(),tempS.end(),rgxItt);
                sregex_iterator end;
                while(next != end){
                    smatch match = *next;
                    tempEvent->fParticleTypes.push_back(stoi(match[0].str()));
                    next++;
                }
                i++;
            }
            if(header->fData_entries["Energies"] == 1){
                G4String tempS = matchName[i].str();
                regex rgxItt(dbl);
                sregex_iterator next(tempS.begin(),tempS.end(),rgxItt);
                sregex_iterator end;
                while(next != end){
                    smatch match = *next;
                    tempEvent->fEnergies.push_back(stod(match[0].str()));
                    next++;
                }
                i++;
            }
            if(header->fData_entries["Translation"] == 1){
                G4String tempS = matchName[i].str();
                regex rgxItt("("+dbl+")(?:"+sla+")("+dbl+")(?:"+sla+")("+dbl+")");
                sregex_iterator next(tempS.begin(),tempS.end(),rgxItt);
                sregex_iterator end;
                while(next != end){
                    smatch match = *next;
                    tempEvent->fTranslation.emplace_back(stod(match[1].str()),stod(match[2].str()),stod(match[3].str()));
                    next++;
                }
                i++;
            }

            if(header->fData_entries["Direction"] == 1){
                G4String tempS = matchName[i].str();
                regex rgxItt("("+dbl+")(?:"+sla+")("+dbl+")(?:"+sla+")("+dbl+")");
                sregex_iterator next(tempS.begin(),tempS.end(),rgxItt);
                sregex_iterator end;
                while(next != end){
                    smatch match = *next;
                    tempEvent->fDirection.emplace_back(stod(match[1].str()),stod(match[2].str()),stod(match[3].str()));
                    next++;
                }
                i++;
            }
            if(header->fData_entries["ParticleTime"] == 1){
                G4String tempS = matchName[i].str();
                regex rgxItt(dbl);
                sregex_iterator next(tempS.begin(),tempS.end(),rgxItt);
                sregex_iterator end;
                while(next != end){
                    smatch match = *next;
                    tempEvent->fParticleTime.push_back(stod(match[0].str()));
                    next++;
                }
                i++;
            }
            //@@@@ sort into storage structures
            G4int sort = tempEvent->fNewEvent[0];
            if(sort == 2 || sort == 5){
                //@@@@ new exposure
                // take care of previous
                if(!vStore.empty()){
                    vvStore.push_back(vStore);
                    eventList.push_back(vvStore);
                    vStore.clear();
                    vvStore.clear();
                }

                // start new
                vStore.push_back(tempEvent);
            }
            else if(sort == 1 || sort == 4){
                //@@@@ new primary
                if(!vStore.empty()){
                    vvStore.push_back(vStore);
                    vStore.clear();
                }
                vStore.push_back(tempEvent);
            }
            else if(sort == 0 || sort == 3){
                //@@@@ same primary
                vStore.push_back(tempEvent);
            }
            else{
                G4cout<<"ERROR: Damage classification not as expected! was "
                      <<tempEvent->fNewEvent[0]<<" but expected 5-0."
                      <<G4endl;
                exit(EXIT_FAILURE);
            }
            nextNewDamage++;
        }
    }
    if(!vStore.empty()){
        vvStore.push_back(vStore);
        eventList.push_back(vvStore);
    }
}
