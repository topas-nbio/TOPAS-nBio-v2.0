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

#pragma once

#include "DrDamageHeader.hh"


typedef std::vector<std::vector<std::vector<DrDamageEvent*> > > vEvents;

class DrReaderSDD {
public:
    DrReaderSDD();
    ~DrReaderSDD();

    //@@@@ Use this function to read in a a SDDv1.0 file.
    /*
     * To use this function the user must provide a seperate class which
     * will store the information in the form of a DrDamageHeader object
     * and a list of DrDamageEvent objects structured as vEvents typedef
     * above. The user can use the TOPAS-nBio class of
     * TsDamagePhaseSpaceStore for this purpose but it has been kept
     * seperate to preserve choice in how to structure the storage class.
     *
     * Once the reader has populated the data structures they should be
     * accessed through whatever storage class the user has provided as
     * the reader itself does not return anything. This is done so all
     * memmory management can be done in the storage class.
     *
     * An example of how this will look:
     *  -----------------------------------------------------------------
     *  auto reader = DrReaderSDD();
     *  auto TsPSS = TsDamagePhaseSpaceStore::Instance();
     *  reader->ReadSTDInput(fileName,
     *                      TsPSS->GetDamageHeaderPtr(),
     *                      TsPSS->GetDamageEventStore());
     *  auto header = TsPSS->GetDamageHeaderPtr();
     *      //DO STUFF WITH HEADER
     *  auto eventList = TsPSS->GetDamageEventStore();
     *      //DO STUFF WITH LIST OF EVENTS
     *  delete TsPSS;
     *  -----------------------------------------------------------------
     *
     *  On the vEvent structure:
     *  The innermost vector stores all damages associated with the same
     *  primary. The next vector out stores the data for all the primaries
     *  associated with the same exposure. The final vector stores each
     *  exposure. So an event list is:
     *      eventList[Exposure][Primary][DamageEvent]
     */
    void ReadSTDInput(G4String fileName,DrDamageHeader* header, vEvents& eventList);

    //@@@@ Reads and parses the header fields. Still some work to be done
    //@@@@ in updating the regex as some fields are still just stored as
    //@@@@ free text rather than being parsed into actual parameters.
    void ReadHeader(G4String fileName, DrDamageHeader* header);

    //@@@@ Reads and parses the damage events according to the header
    //@@@@ field defining the data entries.
    void ReadDamage(G4String fileName, DrDamageHeader* header, vEvents& eventList);
};
