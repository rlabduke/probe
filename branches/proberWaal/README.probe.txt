README.probe                                                 JMW - 7/25/01

The program "probe" generates "contact dots" at points on the van der Waals
surface of atoms which are in close proximity to other atoms [1]; reading
atomic coordinates in protein databank (PDB) format files and writing
color-coded dot lists (spikes where atoms clash) for inclusion in a
kinemage.

Directly based on the "sp" program by Zalis and Richardson (following the
work of Connolly), the approach is to place a small probe (typically of
radius 0.25A) at points along the van der Waals surface of a selected set
of atoms and determine if this probe also contacts atoms within a second
"target" set. A flexible method for selecting the source and target atoms
is available along with command line flags for altering the probe radius
and dot density. Although probe can generate "surface dots" were there are
no nearby atoms, its primary use is to analyze atomic packing. For packing
analysis and structure validation, probe can generate contact surfaces
within a set of atoms ("SELF dots").

For meaningful use of probe in the study of molecular structures,
coordinates for all hydrogen atoms must be included in the model. Modeling
with "implicit hydrogens" is inadequate since the vast majority of steric
interactions which constrain conformational choices take place among
hydrogens. A program called reduce, also available from the Richardson lab,
uses simple geometric considerations to add hydrogens to a PDB file and
optimize their orientations [2].

Probe has many options which modify the way output is formatted. Instead of
kinemage format, it can write graphical information in O or XtalView format.
It can calculate a table of dot information with contact score values and
percent dot coverage. Finally, it can produce a detailed "unformatted"
description of each dot, including source and target atom names, distances,
atom types, and partial scores. Because probe is very flexible, it is
helpful to develop a working knowledge of options and especially selection
criteria.

NEW FEATURES:
v2.4 (7/11/01) - Buttons are no longer generated for each element type
                 by default. To generate these buttons use the -element flag.

v2.5 (7/25/01) - Selections can now refer to negative residue numbers.
                 (Sometimes you need to include extra parentheses or a
		 space to prevent the selection from being treated as
		 a command line option.)

v2.6 (10/28/2011) - Introduced the -condense flag, which when used with the -u option
                    will give one line per source atom - target atom pair. Also it will give 
                    the number of dots in that interaction, as an additional column.
		    (Changes made by Swati Jain)

USAGE:
Probe was designed for UNIX and the commands described below follow the
UNIX conventions. For a brief description of probe features, run "probe"
without any options. The command "probe -h" will give a more complete
description of program options.

In its most basic form, the syntax is

 probe input.pdb >> outputDots.kin

which will generate SELF dots for all atoms in the input file except
alternate (e.g., B or C) conformations and append them to the end of the
kinemage file. Note the ">>" redirection symbol which stands for APPEND;
in the normal case, prekin would be used to make a kinemage of the
molecular structure and probe would be used to append the dot information.

A more extensive set of command line options is available in the format

 probe [-flags] "pattern1" ["pattern2"] input.pdb [more.pdbs...] [>>outfile]

(the parts in square brackets may be optional; by default the results go to
standard output).

There are four modes set by command line flags (-SELF is assumed if not given):

 probe -SELF "pattern1"            inputfiles >>kinfile  # Intersect 1->1 
 probe -BOTH "pattern1" "pattern2" inputfiles >>kinfile  # Intersect 1->2 and 2->1
 probe -ONCE "pattern1" "pattern2" inputfiles >>kinfile  # Intersect 1->2
 probe -OUT  "pattern1"            inputfiles >>kinfile  # External surface 
 
(How the selected atoms interact is listed above as a comment after the
hash mark.)

By default, HET groups and waters are included in the dot calculations
but *NOT* mainchain to mainchain interactions. These settings may be
changed with the -NOHET, -NOWATER and -MC flags.

The flag -U is used to dump 'unformatted' dot information which can
be sent to other programs or scripts for analysis. The flag -STDBONDS
will make probe consult an internal table when deciding the bonding
pattern. This is used in modeling where impossible conformations may
be analyzed without the problem of improper bonding patterns being
inferred from atomic distances.

PATTERNS:
The use of patterns to specify the interaction being examined is
illustrated with the following examples:

  probe "altA blt40" 1filH.pdb >>lowBdot.kin

calculates self packing in all atoms from the file 1filH.pdb with a
temperature factor less-than 40 and an alternate conformation code of
blank or "A" (-self is the default and the pattern is in quotes because it
contains a space). This is a useful pattern for validating a structure
because it ignores atoms which may have poorly determined coordinates. In
other situations, the pattern could be replaced with "all" to select all
the atoms.

To identify the interface between chain E and chain I in the file enzH.pdb

  probe -both "chainE" "chainI" enzH.pdb >> interface.kin

To create a table of contact statistics use -count.

  probe -count -self "all" mypdbH > dotinfo.table

The example also shows the use of a single '>' mark; the UNIX signal to
overwrite (!) rather than append to the output file.

Even more dot information for each dot can be tabulated with -unformated

  probe -unformated -self "all" mypdbH > rawinfo.table

If you need just one line/sot per source atom - target atom pair, use -condense
 
 probe -unformated -condense -self -mc MC mypdbH > reducedinfo.table

You can create surface dots

  probe -out all 1filH.pdb >> surfacedots.kin

These dots are equivalent to the non-reentrant part of a Connolly surface.
When using surface dots, it is sometimes useful to expand or contract
the probe radius using the -rad#.# flag (e.g. -rad1.4 for a water size probe,
or -rad0.0 to see a sphere-like representation of residues).

Finally, here is a sequence of prekin and probe commands which can create
a kinemage where each category of contact is broken down separately. The
patterns used give some sense of the level of control probe permits.

   prekin -lots input.pdb outputdot.kin
   probe -3 -lens -q -name scsc -self "sc alta blt40 ogt33" input.pdb >> outputdot.kin
   probe -3 -lens -q -name scmc -both "sc alta blt40 ogt33" "mc alta blt40 ogt33" input.pdb >> outputdot.kin
   probe -3 -lens -q -name mcmc -mc -self "mc alta blt40 ogt33" input.pdb >> outputdot.kin
   probe -3 -lens -q -name wathet -both "het,water alta blt40 ogt65,(not water ogt33)" \
              "not(het,water) alta blt40 ogt33" input.pdb >> outputdot.kin

AUTOBONDROT

Starting with version 2.0, probe has been extended to read specially marked
up fragments of a PDB file which describe dihedral rotations as well as
other transformations. The command line flag -autobondrot preceeds the
filename and causes the file to be interpreted as a script for scanning
a range of conformations. A description of the format of these rotation
scripts or .rotscr files is in the file README.autobondrot.txt.

TROUBLESHOOTING

If you don't have probe or if you have an old copy, you can get the latest
release from ftp://kinemage.biochem.duke.edu; binary executable files are
available for several operating systems along with source code. Make sure
you download .exe or .tar or .gz files as BINARY. If you download an .exe
file, you will probably wish to rename it to just "probe" and put it into
a directory which is listed in your PATH environmental variable. For UNIX
or LINUX you will also have to make it executable with the command:
	chmod +x probe

The source code should compile easily on almost all UNIX like systems.
Copy the "Makefile" for your system (cp Makefile.xxx Makefile), check
for any system specific issues (altering as required) and then type
	make probe

The most common problem using probe is specifying selection patterns.
Remember self dots (-self, the default) takes one pattern and then
the filename, while interface dots (-both) take two patterns before
the filename.

Output from probe is generally designed to be appended to the end of a
kinemage file. If you just want to see the dots without creating a
model first, add "@kinemage 1" as the first line to the dotfile (either
by hand or using probe flag: -kin) and mage can now display it.


REFERENCES

1) Word, et. al. (1999) Visualizing and Quantifying Molecular
   Goodness-of-Fit: Small-probe Contact Dots with Explicit Hydrogens,
   J. Mol. Biol. 285, 1711-1733.

2) Word, et. al. (1999) Asparagine and Glutamine: Using Hydrogen Atom
   Contacts in the Choice of Side-chain Amide Orientation, J. Mol. Biol.
   285, 1735-1747.

CONTACTS

We hope this helps you get started looking at molecular contact surfaces.
To find the latest version of probe, see http://kinemage.biochem.duke.edu.
A comprehensive description of the small-probe method is found in [1] and
[2].

Mike Word and David Richardson
(mike.word@duke.edu, dcr@kinemage.biochem.duke.edu)

Richardson Lab
Biochemistry Department
Duke University
Durham, NC USA 27710
