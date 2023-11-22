/*

					SOCKET
					 v3.02

					31-10-01
				     John Walshaw

		School of Biological Sciences, University of Sussex,
		Falmer, Brighton, East Sussex BN1 9QG, United Kingdom

			Funded by The Medical Research Council

			   Now at the John Innes Centre
			     john.walshaw@bbsrc.ac.uk


Description: this C program identifies coiled-coil motifs in a Protein Data Bank
(Berman et al, Nucleic Acids Res v28 pp235-42, 2000) file. Also required as
input is a DSSP (Kabsch & Sander, Biopolymers v22 pp 2577-637) file derived from
the PDB file.


*/

/*					preproc.h
					---------
*/


/* last modified: 11-3-1
	#define AMINO_ACIDS 21 added here, so it need be declared only once
	(preproc.h is #include-d by statchar.c,	and must also be #include-d
	before aminoa1.h is parsed)
*/

/*				2 sections:
		PRE-PROCESSOR STANDARD LIBRARY 'INCLUDE'S
		PRE-PROCESSOR CONSTANTS
*/


/* PRE-PROCESSOR STANDARD LIBRARY 'INCLUDE'S: ------------------------------- */


#include <stdio.h>
#include <ctype.h>
#include <math.h>


/* PRE-PROCESSOR CONSTANTS:  ------------------------------------------------ */

#define MAX_FILE_NAME_LENGTH 100	/* maximum no. of characters in the
					file-name strings */

#define MAX_LINE_WIDTH 140		/* this should be at least as great as
					the number of characters in the longest
					line that will be read in from any input
					file - used by fgets() statements */

/* the next 14 constants specify the sizes of arrays which hold various
	structural features; they therefore represent the maximum number
	of each feature which can be handled by the program */

#define MAX_ATOMS 100000			/* the maximum number of atoms which are
					in alpha-helices or within E residues of
					the ends of alpha-helices, where E is
					the helix-extension specified by -e */

#define MAX_RESIDUES_ALL 100000		/* the maximum number of residues of
					*any* type of secondary structure; used
					only to size local arrays in function
					pre_parse_dssp(), which is called only
					once */

#define MAX_RESIDUES 20000		/* the maximum number of residues which
					are in alpha-helices or within E
					residues of the ends of alpha-helices,
					where E is the helix-extension specified
					by -e */

#define MAX_HELICES 1000			/* the maximum number of alpha-helices*/

#define MAX_HELIX_LENGTH 600		/* the maximum number of residues in a
					single alpha-helix */

#define MAX_KNOBS 6000			/* the maximum number of knobs */

#define MAX_KNOBS_PER_HELIX 600		/* the maximum number of knobs in a
					single helix; effectively, this is the
					maximum length of a single uninterrupted
					helix which can be dealt with */

#define MAX_COMPKNOBS 4			/* the maximum number of 'complementary'
					knobs which can be stored, for a given
					knob; both pairwise- and cyclic-
					complementary knobs are stored for each
					knob */

#define MAX_PAIRS 150			/* the maximum number of pairs of
					adjacent helices which pack with knobs-
					into-holes; all such pairs are assessed
					and then boiled down into assemblies
					of N helices, where N is the order of
					the coiled coil */

#define MAX_DAISY_CHAINS 300		/* a 'daisy-chain' is an arrangement of
					N cyclic-complementary knobs, ie one
					knob per helix, where N > 2, forming a
					layer perpendicular to the coiled-coil
					axis */

#define MAX_DAISIES 40			/* the maximum number of daisies (knobs)
					in a daisy-chain (see above); the most
					observed in a solved coiled-coil
					structure is 5 to date; I reckon 6 is
					highly unlikely, because the angle
					subtended by the two knobs-into-holes
					seams, offset by a single residue
					position in the sequence, on a helix
					surface, is close to 90 (tetramer) and
					108 (pentamer) but not 120 (hexamer) */

#define MAX_COILED_COILS 250		/* the maximum number of coiled coils;
					this however is not simply the number of
					'true' coiled coils, but also the
					number of 'sub-assemblies' within them
					(which are identified first); sub-
					assemblies which are wholly part of
					larger assemblies are subsequently
					ignored N.B. this number really ought
					to be at least as large as MAX_PAIRS,
					but currently is not */

#define MAX_HELICES_PER_COIL 40		/* is effectively the same as
					MAX_DAISIES */

#define MAX_DUPLICATES 50		/* it is possible for a large sidechain
					to fit into more than one hole of 4
					sidechains, if the packing-cutoff is on
					the high side (the default of 7.0Å is
					quite liberal); almost certainly the
					2 holes would each be on a different
					helix, and such arrangements are not
					part of classic coiled coil architecture
					(and are also very rare);
					this constant is the maximum number of
					residues which are 'duplicate' knobs
					*/

#define MAX_SEQCHAINS 63		/* the maximum number of PDB chains
					which can be handled; it is 63 because
					its the alphabet in both cases, plus
					ten digits plus ' '; strictly, upper
					and lower case chain identifiers should
					represent the same thing, but in some
					non-native PDB files generated by
					applying symmetry operations to a PDB,
					there are > 26 chains and one convention
					is to use lower case chain identifiers
					as well as upper */


/* the next 4 constants represent miscellaneous features */

#define RASMOL_COLOURS 8		/* if a RasMol script is written, each
					helix is rendered in a different colour;
					when all the colours are used up, the
					next helix colour wraps to the first one
					in the list; the colour which are
					recognized by RasMol could well be very
					platform/graphics-card-dependent; also,
					the number of available colours can
					change dynamically depending on how
					many colours have already been mapped
					by other windows on the desktop, e.g. in
					an X or OpenGL environment; moral is, to
					keep RASMOL_COLOURS small, and use
					mainly primary and secondary colours */

#define RASMOL_WRAP 8			/* if a RasMol script is written,
					various 'sets' are defined, e.g. which
					residues belong to which helix, which
					knobs belong to which helix, which
					helices	belong to which coiled coils,
					etc; some of these definitions require
					lists of many comma-delimited elements,
					but there is a limit to the length of a
					command-line (255 characters?); so
					RASMOL_WRAP defines the maximum number
					of elements written per line; each of
					these is defined as a subset, and then
					the correct set defined as the union of
					these subsets; this causes problems if
					the number of elements is >
					RASMOL_WRAP^2, so it should not be set
					too small; but if its set too big, then
					the usual problem of overrunning the
					permitted command line length occurs */

#define DEFAULT_CUTOFF2 7.0		/* a crucial constant: this is the
					value of the packing-cutoff used if the
					program is run without the '-c' option
					which overrides it; units are Ångstroms
					*/

#define DEFAULT_KNOB_THRESHOLD 3	/* this is probably the most crucual
					constant: it defines the difference
					between a 'true' knob, i.e. one which
					is involved in a complementary knobs-
					into-holes interaction, and a 'one-way'
					knob; initially, all knobs are
					identified and assigned as either type 1
					or type 2 (across- or in-hole) and then
					those which are involved in
					complementary interactions with each
					other are promoted to type 3 or 4 (or
					even higher, in some rare non-canonical
					arrangements); from then on, types 1 and
					2 are completely ignored; this constant
					defines the lowest knobtype which is
					'true', i.e. complementary; it should
					never be changed (but can be with '-k',
					which is a legacy of earlier versions)*/


#define SPACER_PATTERNS 2		/* this constant is redundant */

/* the next 3 constants relate to the system for interpreting command-line
arguments (which may not be the most straightforward method in C, but was ported
from a very flexible Perl implementation, which allows flags and/or parameters
of any name-length to be specified in any order, with or without spaces
between parameter-names and their values; anyway, it works) */


/* the number of amino acids (used by functions relating 1-letter
and 3-letter codes) */

#define AMINO_ACIDS 21

/* the number of heterogen-type non-standard amino acids; e.g. CSE,
MET, which are specified in PDB files as HETERO atoms but correspond
to a standar amino acid (CYS, MET respectively) */

#define HETEROGENS 74

/* the maximum number of heterogen-type non-standard amino acids (see
above); up to HETEROGENS_MAX - HETEROGENS extra types can be read
from MODRES records in the PDB file (even more than 2 in one file
will be pretty rare, however) */

#define HETEROGENS_MAX 84

/* 27-6-1 the number of known solvent heterogens (i.e., don't appear
as part of the polypeptide chain). This is not crucial, but is so
that such heterogens are recognized, so that unknown amino-acid
heterogens can be flagged */

#define NON_AA_HETEROGENS 1012


/* the null seqchain[] value (PDB serial no , i.e. resSeq, of first residue of
	a chain) */

#define NULL_SEQCHAIN -99999

#define MAXQUALNAMELENGTH 20		/* used by argument.h - the
					maximum number of characters in the
					name of a command-line argument;
					ought to be specified in that header
					file, instead of here */

#define FLAGS 9				/* the number of available 'flags' i.e.
					binary command-line arguments, specified
					by '-<flagName>' */

#define PARS 11				/* the number of 'parameters', i.e.
					command-line arguments coupled with
					values, specified by
					'-<parameterName> <parameterValue>' */


/* this is the banner stamped on several of the output files */

#define SOCKET_TITLE "SOCKET v3.02 02-11-01 John Walshaw, University of Sussex\n"

