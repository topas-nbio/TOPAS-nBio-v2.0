# DaMaRiS
```
********************************************************************************
Software: DaMaRiS(nBio) v2020.07.23
Contact: Dr. John-William Warmenhoven
Contact email: john.warmenhoven@manchester.ac.uk
Initial Date: 26/01/2016
Last Updated: 23/07/2020

PRECISE (Proton Research at the Christie & Division of Cancer Sciences)
Division of Cancer Sciences
School of Medical Sciences
Faculty of Biology, Medicine and Health
University of Manchester
UK
Website: www.bmh.manchester.ac.uk/research/domains/cancer/proton/

This  code  implementation is the result of  the  scientific and technical work
of the PRECISE Group in the Division of Cancer Sciences at the University of
Manchester. By using,  copying,  modifying or  distributing the software (or any
work based  on the software)  you  agree  to  acknowledge its use  in resulting  
scientific publications. To acknowledge this work please cite the following
publications:

DOI: 10.1016/j.dnarep.2019.102743
DOI: 10.1038/s41598-018-21111-8
DOI: 10.1038/s41598-019-42901-8

This work is funded by the Manchester Biomedical Research Centre, EPSRC Proton
Therapy Network EP/N027167, BioProton EPS0243444, and EP/J500094. STFC Advanced
radiotherapy Network and EU Integrating Activity INSPIRE (Grant No. 730983).
********************************************************************************
```

## Quick Start

To run a DaMaRiS example navigate to examples/damaris and use the following
command:

`/path/to/topas/executable ./DaMaRiS.run`

runDr.run is a topas parameter file used to set up and run DaMaRiS. It is
extensively commented to give you an idea of what parameters you can play around
with to alter various aspects of the simulation.

## DaMaRiS User Settings

### Diffusion Parameters
These parameters control how the DSB molecules diffuse through the simulation
and are included from the file *Parameters_DaMaRiS_Motion.txt*:

`i:Ch/DaMaRiS/DiffusionMode = 1`

**0** sets the motion to normal Geant4-DNA Brownian motion with default
diffusion coefficients of 1.4 nm^2/s for all DSBs. **1** sets the motion to
sub-diffusion implemented in the DREP code as a CTRW model.

`u:Ch/DaMaRiS/DiffusionCoefficientForJump = 2.808e11`

This is the diffusion coefficient used by the molecule to determine its
motion during the one 'jump' it is allowed to do between waiting times.
This is probably the only parameter which should be changed in order to
modify the speed of sub-diffusion for the DSBs.

`u:Ch/DaMaRiS/DiffusionCoefficientForTrapped = 0`

This is the diffusion coefficient used by the molecule to determine its
motion during whilst it is waiting. This could be set to something small
to introduce some amount of 'wiggle' whilst the particle is waiting BUT
then be sure to properly characterise the motion and see if it still is
sub-diffusive.

`d:Ch/DaMaRiS/MinWaitingTime = 1e9 ps`

Minimum time the particle is allowed to be trapped for. Larger values will
speed up the simulation time but also reduce the overall mobility of DSBs.

### Geometry Parameters

`d:Ch/DaMaRiS/BoundingCellOrNucleusRadius = Ge/Target/RMax um`

This is important to set correctly! This is the radius of the sphere which the
DSBs are not allowed to leave (ie. 'the nucleus'). I have not been able to
confine them in any other manner so for now no other shape than a sphere is
allowed for the nucleus and it has to be centered at the world origin.

### Simulation Parameters
`i:Ch/DaMaRiS/AlternativeRunID = 0`

This is a number which is appended onto the end of all output files from DaMaRiS
runs.

`i:Ch/DaMaRiS/BiologyRepeatNumber = 10`

How many times do you want to repeat the repair simulation. 30+ seems to give OK
statistics.

`d:Ch/DaMaRiS/DaMaRiSStageTimeEnd = 86400 s`

How long should the repair simulation go on for. With current time constants
Pre-synaptic recruitment has reached a maximum at around ~30s and most DSB
repair happens within ~30 min. I would consider 24 hours a standard run. The
duration is extended in the simulation by 0.1% to ensure last time point is
captured.

`s:Ch/DaMaRiS/PathwayFileName = "pathway_NHEJ.in"`

Sets the biological pathway to be used in the simulation. pathway_NHEJ.in will
simulate repair by NHEJ only pathway_HR.in will simulate repair by HR and NHEJ

### Measurement Parameters

`d:Ch/DaMaRiS/ObserveDurationForMSD = 300 s`

`d:Ch/DaMaRiS/ObserveStepSizeForMSD = 1 s`

These parameters are used to set the time frame and frequency over which the
motion of DSBs are sampled in order to investigate their MSD. For sub-diffusive
motion this should scale as MSD = t^a, where a < 1. From literature I would
expect a to be around 0.5 in order to mimic DSB motion in real cells.

`s:Ch/DaMaRiS/ExplicitBinningFileName = "Null"`

If the user wishes to output the state of the system at specific times then this
parameter can be used to point DaMaRiS at a file containing those times.
Otherwise the default bins will be set as:
- From 0 second to 1 minute in 1 second steps.
- From 1 minute to 10 minutes in 10 second steps.
- From 10 minutes to 1 hour in 1 minute steps.
- From 1 hour till the end in 10 minute steps.

### Damage Placement Parameters
These settings allow you to place damage into the simulation through various
means with the standard, expected, method being to place from a SDD file.
DaMaRiS will only use *one* of these methods. If no method is supplied it is
assumed that there is a preceding damage simulation which will use a damage
phase space store to communicate DSB placements to DaMaRiS.

**1) From File**

`s:Ch/DaMaRiS/STDFormatDamageFileName = "damage.in"`

If you are reading in from a file the name of the file must be specified here.
DaMaRiS will store the parsed damage data extracted during the first run in
memory and use that in subsequent repeats. The damage file must be in SDD v1.0
format.

`b:Ch/DaMaRiS/turnOffTime = "True"`

Turns on/off placing DSBs using the timing information in the damage file

`i:Ch/DaMaRiS/SelectFromExposureNumber = -1`

If the damage file contains multiple exposures the default behaviour is to
select one of them at random to populate the repair simulation. If instead you
wish to investigate a specific exposure the it can be selected here. Numbers
start at 1 for the first exposure identified in the file. Setting the number to
-1 will cause random selection.

**2) DSB *Ends* at Origin**

`i:Ch/DaMaRiS/DSBOriginNumber = -1`

If DSBOriginNumber is >= 0 then this number of "DSBEnd" objects will be built at
 the origin.

**3) DSB *End* with Offset**

`d:Ch/DaMaRiS/DSBOffset = 0.0 nm`

If DSBOffset > 0.0 then this wil build a single "DSBEnd" object the specified
number of nm away from the origin on the x axis.

**4) *DSBs* in a Column**

`i:Ch/DaMaRiS/DSBColumnNumber = -1`

`d:Ch/DaMaRiS/DSBColumnRadius = -1.0 nm`

If DSBColumnNumber >= 0 then that number of DSBs will be built randomly in a
column of specified radius along the z axis through the nucleus.

**5) *2 DSBs* with a Specified Separation and Delay**

`d:Ch/DaMaRiS/DSBSeparation = -1.0 nm`

`d:Ch/DaMaRiS/DSBTimeDelay = 0.0 s`

If DSBSeparation >= 0.0 this will place two DSBs that far apart with a specified
delay in placing the second DSB

## Code Repository

DaMaRiS is a developing framework and will be continually developed. For the
latest version of the code please check at:

https://gitlab.com/PRECISE-RT/releases/damaris_topas-nbio

## References
* DOI: 10.1039/C8RA10168J						       
* DOI: 10.1038/s41598-018-21111-8					       
* DOI: 10.1038/s41598-019-42901-8
* DOI: 10.1016/j.dnarep.2019.102743
