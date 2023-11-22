/*

					SOCKET
					 v3.03

					global.h

					13-7-06
				  John Walshaw

		School of Biological Sciences, University of Sussex,
		Falmer, Brighton, East Sussex BN1 9QG, United Kingdom

			Funded by The Medical Research Council


PLEASE REFERENCE:	Walshaw, J. & Woolfson, D.N. (2001), "SOCKET: ",
		Journal of Molecular Biology, ... (), pp ...-...


Description: this C program identifies coiled-coil motifs in a Protein Data Bank
(Berman et al, Nucleic Acids Res v28 pp235-42, 2000) file. Also required as
input is a DSSP (Kabsch & Sander, Biopolymers v22 pp 2577-637) file derived from
the PDB file.



*/

/*					global.h
					--------
*/



/*	 			6 sections:
		GLOBAL VARIABLES: enumerated types
		GLOBAL VARIABLES: command-line interpreter
		GLOBAL VARIABLES: side chain-position descriptors
		GLOBAL VARIABLES: 'relational database' of structural entities
		GLOBAL VARIABLES: miscellaneous
		GLOBAL VARIABLES: file-handles
*/

#include "aminoa1.h"

/* GLOBAL VARIABLES: enumerated types --------------------------------------- */

/* 3 enumerated types: one representing the available flags, the second the
available parameters, and the third a boolean type; N.B. the actual names of
the flags/parameters the user specifies on the command-line are defined in
flagname[] and parname[] */

enum flags {flag_debug,flag_a,flag_b,flag_i,flag_l,flag_p,flag_q,flag_u,flag_v};
enum pars {par_c,par_d,par_e,par_f,par_j,par_k,par_o,par_r,par_s,par_t,par_w};
enum boolean {false,true};

/* flags:
	debug	- writes out loads of very verbose info
	a	- use all atoms, not just heavy (ie include hydrogens if any),
			for determining side chain centres-of-volume
	b	- not used - previously, was used to select an alternative
			algorithm for 'bundling' the coiled coil regions into
			coiled coils
	i	- include C-alphas as part of a side chain
	l	- list all the knob-hole interactions and knob complementarity
			- and numerous other details
	p	- 'private' output - names of input/output files are not
			mentioned
	q	- 'quiet' output
	u	- check unmatched heterogen names against a list of known
		  heterogens
	v	- verbose output - lists the helices, the helical residues
		  and their centres of volume and ends; and many other details

parameters:
	c	- packing-cutoff (Ångstroms)
	d	- the name of the 'summary' file; this feature is of very
			limited use, has not been properly supported in recent
			versions and is UNSTABLE - legacy feature, don't use it!
	e	- helix extension (in residues); each helix is extended by e
			residues at each end (if these residues exist)
	f	- name of input PDB file (mandatory)
	j	- NOT USED
	k	- the lowest knob-type defined as complementary - should always
			be 3 - see DEFAULT_KNOB_THRESHOLD - legacy feature,
			don't use it!
	o	- name of the optional long output file, one residue per line
	r	- name of the optional output RasMol script file
	s	- name of the input DSSP file (mandatory)
	t	- name (4-characters) of the atom type used as reference atom
			- (0) see below; default is REFATOM0, and changing it to
			anything else is unlikely to be of much use (legacy
			feature)
	w	- the name of the 'helix register' file; now redundant, as
			nothing is actually written to it- legacy feature, don't
			use it!
	
*/


/* GLOBAL VARIABLES: command-line interpreter ------------------------------- */

/* the names of the available flags, as specified by the user on the command-
line when preceded directly by '-' */

char *flagname[FLAGS]; /* dummy declaration, so that flagname[] is known to
				all source files which #include global.h ; the values
				cannot be initialized here, because that must be done
				only once, or else the compiler gets upset; it is
				done in statchar.c */


/* corresponding to each of the elements of flagname[] is this boolean array
(each flag is either 'on' or 'off') */

enum boolean setflag[FLAGS];


/* the names of the available parameters, as specified by the user on the
command-line when preceded directly by '-' */

char *parname[PARS]; /* dummy declaration, so that parname[] is known to
				all source files which #include global.h ; the values
				cannot be initialized here, because that must be done
				only once, or else the compiler gets upset; it is
				done in statchar.c */



/* this array holds the actual values corresponding to each element of
parname[]; all values are strings, but non-string variables are later read from
them with sscanf() if appropriate */

char *par[PARS];


/* GLOBAL VARIABLES: side chain-position descriptors ------------------------ */

/* describing the position of a side chain:
	the initial concept is to use 3 descriptors:
	(0) is the coordinate of the end of the side chain nearest the backbone
	(1) is the coordinate of the furthest end of the side chain
	(2) is the coordinate of the middle of the side chain

	(0) is the only one of the three which is always an explicit atom, which
	is the alpha-carbon (CA); see REFATOM0 below; N.B. alpha-carbons are
	used by default, but a different atom type can be specified with '-t'

	- however, it turned out that (2) alone was capable of distinguishing
	between knob-into-holes packing and non-knobs-into-holes packing;
	so (0) and (1) were not used

	later features introduced were:
		distinguishing between knobs-into-holes and
		'knobs-across-holes', which requires the end of the side chain,
		so (1) was used for this;
		calculating core-packing angles, which requires alpha-alpha
		carbon vectors for the holes, and alpha-beta carbon vectors
		for the knobs; so (0) was used for this, and a fourth descriptor
		was introduced, the beta-carbon (CB) - see REFATOM3 below; like
		(0) this is an explicit atom

	(1) vs (1b): numerous amino acid R-groups are forked at the end, and so
		they have 2 'end-atoms', so each residue has a 2D array, though
		the second is not always used; if both are used, then 1b is the
		mean coordinate of the 2; otherwise 1b is the same as the
		coordinate of the sole end-atom (for R-groups unbranched at the
		end)

	(2)	the 'middle' of the side chain is the mean coordinate of all the
		heavy (i.e. non-hydrogen) side chain atoms; the alpha-carbon is
		not included as a side chain atom unless the residue is Glycine,
		(in which case (0), (1) and (2) are all the same), or the '-i'
		option is used; if hydrogen atoms are present in the PDB
		structure file, they can be included in the calculations by
		using the '-a' flag

	(3)	is the beta-carbon atom

	(0),(1) and (3) are represented by integer arrays
		- their values are atom-keys (serial no's of atoms)
	(1b) and (2) are represented by real arrays
		- their values are x, y and z coordinates
	
*/

/* a different type of reference atom (0) can be specified by using -t; it is
unlikely that this would ever be of much use, however (its a legacy feature);
REFATOM0 defines the *default* reference atom type (0); the type which is used
ends up in the string refatom0type; while the later-introduced REFATOM3 is
unchangeable by the user */

#define REFATOM0 " CA "
#define REFATOM3 " CB "

char refatom0type[5];	/* stores the actual 4-character string which represents
			the atom-type used as reference atom (0), e.g. " CA ",
			which is the default (see REFATOM0) */

int	refatom0[MAX_RESIDUES],		/* refatom0 is the list of C-alpha (or
					REFATOM0) atoms, one for each residue,
					referenced by index of atom array;
					values stored are atom-keys */

	refatom1[2][MAX_RESIDUES],	/* refatom1 is the list of (1 or 2)
					'end-atoms' of each residue, referenced
					by index of atom array (end-atoms for
					each residue type are in array
					refatomtype, below); values stored are
					 atom-keys */

	refatom3[MAX_RESIDUES];		/* refatom3 is the list of C-beta atoms
					(needed for calculating knob-into-hole
					packing geometry); values stored are
					atom-keys */


float	refatom1B[MAX_RESIDUES][3],	/* refatom1B is the list of coordinates
					of the 'pseudo atom' at the end of each
					residue's side chain. If there is only
					one actual end atom (see above), this
					pseudo-atom has the same coordinates as
					this end atom. If there are two actual
					end atoms, then the pseudo atom is the
					mid-point between them. This pseudo atom
					is needed to define a single coordinate
					of a residue's end */


	refatom2[MAX_RESIDUES][3];	/* refatom2 is the list of coordinates
					of the 'centre of mass' of each
					residue's side chain- actually mean
					coords of each side chain atom, without
					weighting by mass */



/* order of amino acids in type "enum amino_acid" is:
 {X,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y} */

static char *refatomtype[2][AMINO_ACIDS] = {
				{"    ", " CB ", " SG ", " OD1", " OE1",
					" CZ ", " CA ", " CE1", " CD1",
					" NZ ", " CD1", " CE ", " OD1",
					" CG ", " OE1", " NH1", " OG ",
					" OG1", " CG1", " CH2", " OH " },
				{"    ", "    ", "    ", " OD2", " OE2",
					"    ", "    ", " NE2", "    ",
					"    ", " CD2", "    ", " ND2",
					"    ", " NE2", " NH2", "    ",
					" CG2", " CG2", "    ", "    "}
			};

/* (the refatomtype array is to identify which atoms' coordinates should be
read into the 2 elements of refatom1) */


/* GLOBAL VARIABLES: 'relational database' of structural entities: -------------
	The structural elements read from the PDB are stored in a large number of
	arrays, which effectively constitute a relational database; the simplest
	element is the atom, and only atoms which are in alpha-helices are
	considered (or within E residues of the ends of helices, where E, the
	'helix-extension', is {0,1,2} residues).

	The structural elements are:
		atom, residue, helix, knob, helix-pair, coiled coil

	-most are described by several arrays, indexed by a key which starts
	at zero; ideally the same integer variable would always be used to
	index the arrays of a particular element, but this is far from true;
	sorry.

	In the following descriptions, the fields in the ATOM records of PDB
	files are often referred to (see PDB Contents Guide)
*/


/* arrays describing atoms ................................................. */

/*	the key of each array is the serial no. of each atom (atom-key),
	starting at zero; only atoms from residues which are in alpha-helices
	are stored, and the rest are ignored completely */

int	atom_no[MAX_ATOMS],		/* atom serial number (PDB 'serial'
						field) */

	atom_res[MAX_ATOMS];		/* the residue-key of the residue to
						which this atom belongs */

char	atom_name[MAX_ATOMS][5];	/* atom name (PDB 'name' field) */

float	coord[MAX_ATOMS][3];		/* orthogonal coordinates of atom
						(PDB 'x', 'y' and 'z' fields) */

/* N.B. there is no GLOBAL record of the total number of atoms which have been
read in; the variable atom_index is local to main(), which is sloppy as all
other 'tables' of the RDB have their size stored as a global variable (except
for the residue table); the original design was to use nearly all local
variables and pass them as function arguments, but this would have lead to some
very long and complicated function calls; a local atom_index is a legacy thing,
and should be put right; sorry */


/* arrays describing residues ............................................... */
	/* these arrays are indexed by the residue-key, which starts at zero */
	/* helix_residue_no and helix_residue_aacode are so named to labour the
	point that they store only residues which are in alpha-helices; there
	are no other types of residues dealt with by the program */

int	helix_residue_no[MAX_RESIDUES], /* residue sequence number (PDB
						'resSeq' field) */

	helix_no[MAX_RESIDUES],		/* the helix-key of the helix to which
						this residue belongs */

	helix_residue_aacode[MAX_RESIDUES],	/* the value of the amino_acid
						enumerated type (see
						aminoa1.h) corresponding to
						this residue */
						
	n_contacts[MAX_RESIDUES],	/* the number of contacts this residue
					makes with other residue(s); a contact
					is when the two refatom2's of a pair of
					side chains in *different* helices are
					separated by no more than the packing-
					cutoff */

	contact[MAX_RESIDUES][10];	/* the n_contacts contacts are stored
					here; up to 10 contacts can be stored,
					but if there are more than 5 then the
					packing-cutoff is almost certainly too
					high; the values are the residue-keys
					of the contacting residues */

char	helix_residue_name[MAX_RESIDUES][5],	/* residue name (PDB 'resName'
						field; this is not entirely
					redundant information with respect to
					helix_residue_aacode, because the latter
					describes one of 21 standard amino
					acids; non-standard residues have a
					helix_residue_aacode corresponding to
					'X' whereas helix_residue_name stores
					the 4-character string (usually a
					3-letter code) actually in the PDB file,
					which is useful for when outputting
					residue	details for the user */

	tad_register[MAX_RESIDUES][MAX_COILED_COILS],	/* the heptad register
							(0..6, corresponding to
						a..g) of the residue, if any; in
						complex assemblies, a residue
						can contribute to more than one
						coiled coil, so its a 2D array,
						indexed by the coiled-coil-key*/

	helix_residue_iCode[MAX_RESIDUES]; /* XXX one for the future... sort this out */

/* N.B. there is no GLOBAL record of the total number of helical residues which
have been read in; the variable residue_index is local to main(), which is
sloppy - see notes on atom_index above */


/* arrays describing helices ................................................ */
	/* these are indexed by the helix-key, which starts at zero */

int	helix_start[MAX_HELICES],	/* the residue sequence number (PDB
					'resSeq' field) of the most N-terminal
					residue of the helix; *NOT* the
					residue-key of that residue */

	helix_end[MAX_HELICES],		/* the residue sequence number (PDB
					'resSeq' field) of the most C-terminal
					residue of the helix; *NOT* the
					residue-key of that residue */

	n_knobs[MAX_HELICES],		/* the number of residues in the helix
					which are knobs of any type */

	n_hole_res[MAX_HELICES],	/* the number of holes in the helix;
					(not the number of residues which
					constitute those holes */

	helix_order[MAX_HELICES],	/* the highest order of any coiled coil
					to which the helix belongs (in complex
					assemblies, a helix can belong to more
					than one coiled coil) */

	n_knobtype[MAX_HELICES][7],	/* the numbers of knobs of each type
					which this helix has; knobs with type=6
					are possible (these are nearly always
					long side chains); the tally is
					cumulative, so that the number of knobs
					of type n1 also contribute to the
					tallies	of type n2, where n2 < n1 */

	n_holetype[MAX_HELICES][7];	/* the numbers of holes in the helix
					containing each type of knob; see
					n_knobtype above */

char	helix_chain[MAX_HELICES],	/* the chain identifier (PDB 'chainID'
					field) */

	helix_start_iCode[MAX_HELICES], /* the insertion code (PDB 'iCode'
					field) of the most N-terminal residue of the
					helix */

	helix_end_iCode[MAX_HELICES]; /* the insertion code (PDB 'iCode'
					field) of the most C-terminal residue of the
					helix */



/* scalars describing helices ............................................... */
int	helix_index;			/* the total number of helices read in;
					this starts at zero and is then
					incremented as each helix is read by
					read_helical_dssp() (if helix extension
					is 0) or pre_parse_dssp (if helix
					extension > 0); after that it is
					constant */


/* arrays describing polypeptide chains ..................................... */

int seqchain[MAX_SEQCHAINS];		/* the (PDB) serial number of the first
					residue of each chain */


/* arrays describing knobs .................................................. */
	/* these are indexed by the knob-key, which starts at zero;
	usually integer variable k is used as the index */

int	knob[MAX_KNOBS], 		/* the residue-key of the residue
					whose side chain is the knob */

	knobtype[MAX_KNOBS],		/* the type of the knob : a null
					value 0, or 1,2,3,4, or in rare
					cases 5 or 6 ('double knobs') */

	hole[MAX_KNOBS][4], 		/* the residue-keys of the 4 residues
					which form the hole into which the
					knob fits; there are always 4 hole
					residues (if a residue has more than
					4 contacts - see above - then the best
					4 are selected to make the hole) */

	n_compknob[MAX_KNOBS],		/* the number of complementary knobs*/

	compknob[MAX_KNOBS][MAX_COMPKNOBS],	/* the knob-key of all the
						knobs complementary to this
						knob */

	knob_order[MAX_KNOBS];		/* the number of knobs in the
					complementary arrangement to which the
					knob belongs; i.e. 2 for pairwise
					complementary, 3,4 or 5 for 3-, 4-,
					or 5- stranded coiled coils; in complex
					assemblies, a knob can be part of more
					than one cyclic arrangement (i.e. more
					than one coiled coil); in this event
					the knob_order is the highest order of
					any of the coiled coils (this frequently
					happens in 4- and 5-stranded, and less
					often in 3-stranded, where a core knob
					is cyclically complementary with other
					core knobs (order > 2), but pairwise-
					complementary with a peripheral knob
					(order = 2) ) */

float  angle[MAX_KNOBS],		/* the core-packing angle of the knob,
					in degrees; see function packing_angle()
					*/

      contact_distance[MAX_RESIDUES][10], /* added by JW 13-7-6.
	                                This records the distances between
					residues (centres-of-volume) which are
				        'in contact', i.e. CoVs separated by
					the cutoff distance or less */

       hole_distance[MAX_KNOBS][4];    /* the distance between the centre of
	                                volume of the knob sidechain and the
					hole sidechain (0..3) JW 13-7-6 */


/* scalars describing knobs ................................................. */
int	knob_index;			/* the total number of knobs (of any
					type) identified, i.e. side chains
					fitting into a hole of 4 side chains on
					another helix */


/* arrays describing sets of duplicate knobs ................................ */
int	duplicate_knobs[MAX_DUPLICATES][3];	/* it is possible for a side
						chain to fit into more than one
					hole (see MAX_DUPLICATES); it is
					inconceivable that a side chain could
					fit into more than 2 holes, unless the
					packing-cutoff were stupidly high, but
					the 2nd dimension has 3 elements anyway;
					two different knob-keys could in fact
					describe one and the same residue, and
					such relationships are store here; each
					value of duplicate_knobs is a knob-key*/

/* scalars describing sets of duplicate knobs ............................... */
int	n_duplicate_knobs;		/* the number of sets of 'different'
					knobs (i.e. different knob-keys) which
					are in fact the same residue (residue-
					key) */


/* arrays describing daisy-chains ........................................... */
	/* these are indexed by the daisy-chain-key, which starts at zero */

int	daisy_chain[MAX_DAISY_CHAINS][MAX_DAISIES],	/* daisies in daisy
							chains; each daisy-chain
						is a list of daisies, i.e.
						knob-keys */

	daisy_chain_cc[MAX_DAISY_CHAINS];	/* the coiled-coil-key of the
						coiled coil to which the daisy
						chain belongs (there is a many
						to one mapping of daisy chains
						to coiled coils) */

/* scalars describing daisy-chains .......................................... */
int	daisy_chains;			/* the number of arrangments of
					cyclically-complementary knobs (i.e.
					layers in coiled coils of order > 2) */


/* arrays describing coiled coils ........................................... */
	/* these are indexed by the coiled-coil-key, which starts at zero */

int	coiled_coil[MAX_COILED_COILS][MAX_HELICES_PER_COIL],	/* the list of
								helix-keys of
						the helices which constitute
						the coiled coil */

	coiled_coil_tally[MAX_COILED_COILS],	/* the number of times that the
						coiled coil (i.e., combination
						of helices) occurs in the lists
						of pairwise-complementaty knobs
						and of daisy chains; a tally of
						only 1 daisy chain is required
						to identify a coiled coil of
						order > 2, whereas the arbitrary
						lower threshold of 2 pairwise
						layers is imposed for
						identifying 2-stranded coiled
						coils (a pair of helices with
						only a single pairwise layer of
						knobs-into-holes is marked
						'IGNORED' in the list of 'coiled
						coils' */

	coiled_coil_helices[MAX_COILED_COILS],	/* the number of helices in each
						coiled coil */

	coiled_coil_subset[MAX_COILED_COILS],	/* some arrangements of helices
						are subsets of others; e.g. in
						many 4-stranded coiled coils,
						adjacent pairs of helices form
						pairwise interactions as well
						as contributing to cyclic
						interactions with other pairs of
						helices; initially, each of
						these pairs is identified as a
						'coiled coil' in its own right,
						but is then later identified as
						being wholly a subset of the
						larger, 4-stranded assembly; the
						coiled_coil_subset of the '2-
						stranded coiled coil' is
						therefore set to the
						coiled-coil-index of the
						4-stranded coiled coil; the
						latter's coiled_coil_subset is
						set to -1 to denote it as being
						a complete, true coiled coil */

	coiled_coil_orientation[MAX_COILED_COILS],	/* the overall
							orientation (0 =
						parallel, 1 = antiparallel) of
						the coiled coil; this is the
						result of logically ORing all
						the orientations of adjacent
						pairs of helices which belong
						to the coiled coil */

	coiled_coil_max_length[MAX_COILED_COILS],	/* each helix in a
							coiled coil has a span
						of knobs-into-holes packing,
						from the most N-terminal knob
						and/or hole residue to the most
						C-terminal; there can be missing
						layers (no knobs-into-holes
						interactions) between them; for
						various	reasons, the span on
						each helix of a coiled coil is
						not necessarily exactly the
						same; coiled_coil_max_length is
						simply the length of the longest
						span (which contributes to this
						coiled coil; the helix might
						have independent heptad(s) which
						are involved in another coiled
						coil, if its a complex assembly)
						of any of the helices */

	coiled_coil_begin[MAX_COILED_COILS][MAX_HELICES_PER_COIL][2],	/* this
									stores
						the most N-terminal limit of the
						span of each helix which belongs
						to the coiled coil; the helix
						is indexed by the order in which
						they appear in the coiled coil,
						NOT the helix-key; there are two
						definitions of the limit,
						corresponding to the two elements
						of the third dimension: 0 holds
						a *knob-key* of the most
						N-terminal knob; 1 holds the
						*residue-key* of the most N-
						terminal residue of any hole
						into which a knob fits */

	coiled_coil_end[MAX_COILED_COILS][MAX_HELICES_PER_COIL][2];	/* same
									as array
						coiled_coil_begin, but defines
						the most C-terminal limits */

float   coiled_coil_mean_length[MAX_COILED_COILS];	/* the mean of all the
							span-lengths of all the
						helices in this coiled coil
						(only counting spans which
						contribute to this coiled coil,
						in the case of complex
						assemblies where a helix belongs
						to more than one coiled coil;
						see coiled_coil_max_length
						(above) */

/* scalars describing coiled coils .......................................... */

int	n_total_ccs,		/* the total number of 'coiled coils' in the
				structure, i.e. the number of arrangments of
				2 or more helices packing with knobs-into-holes;
				some of these arrangements are subsets of larger
				coiled coils, and so are ignored, while others
				are ignored because they have only one layer of
				pairwise interactions and no cyclic interactions
				(see coiled_coil_tally) */

	n_true_ccs,		/* the number of the assemblies (see n_total_ccs
				above) which: have at least 1 cyclic layer, or
				at least 2 pairwise layers, and which are not
				subsets of any other assemblies */

	coiled_coils;		/* a bit of a mess; this is redundant, as it is
				only used by function determine_order(), and
				should be local to it; its value is returned to
					n_total_ccs (above) */

/* GLOBAL VARIABLES: miscellaneous: ----------------------------------------- */

/* scalars describing miscellaneous structure-related parameters ............ */

int	knob_threshold,		/* the lowest knob-type defined as
				complementary; the value used is either defined
				by -k or else DEFAULT_KNOB_THRESHOLD is used;
				really, the value 3 should only ever be used,
				and a constant used instead of this variable
				with -k disabled */

	extend;			/* the helix-extension in residues; the value
				used is either defined by -e or else 0 is used;
				values of 0, 1 or 2 are acceptable */

float	/*cutoff0,*/		/* refer to 'GLOBAL VARIABLES: side
				chain-position descriptors; each side chain is
				defined by the position of its proximal end
				(near the backbone, ie alpha-carbon), its distal
				end and its middle; originally it was
				anticipated that all 3 descriptors would be used to define
				knobs-into-holes packing, but using only the
				middle (descriptor 2) was found to be
				sufficient; so only cutoff2 (in Ångstroms) is
				used; cutoff0 and cutoff1 have never been
				employed; see notes for functions
				measure_CA_distance() and
				measure_end_distance_B() */
	/* cutoff1,*/
	cutoff2;

char cutoff_warning;		/* this is set to a space character; its name
				is for historical reasons */


/* Multi purpose string variables used in several functions.................. */

char	record_type[8],			/* used in read_helical_dssp and
					read_helical_pdb) */

	textstring[MAX_LINE_WIDTH],	/*used mainly in read_helical_dssp and
					read_helical_pdb) */

	code[MAX_FILE_NAME_LENGTH];	/* stores the 'code' identifying the
					PDB structure; if the name of the input
					PDB file (specified by -f) appears to be
					standard, i.e.
					(pdb)?\d\w\w\w(.(pdb|ent))?, then the
					4-letter code \d\w\w\w is used;
					otherwise the name of the file (minus
					any leading path) is used; value is
					determined by function get_code */

/* boolean .................................................................. */
enum boolean is_helical;	/* this is used by functions read_helical_pdb,
				read_helical_dssp and pre_parse_dssp; it really
				ought to be local to them */


/* number of types of amino-acid residues which are known to be recorded as
heterogens in PDB files (i.e. HETATM not ATOM records); the list of these
types is defined in function statchar(), but can be extended if MODRES
records are encountered which specify residue types not in this list; so,
the number is not static */

int	n_heterogens;

/* GLOBAL VARIABLES: file-handles ------------------------------------------- */
	/* These are opened by function check_files */

/* input-files .............................................................. */
FILE	*pdb_file,		/* used by read_helical_pdb function, to read
				atom and residue data */

	*dssp_file;		/* used by read_helical_dssp, and also by
				pre_parse_dssp if the helix-extension > 0;
				the positions of helices in the sequence is
				read by these functions */

/* output-files ............................................................. */

FILE	*data_outfile,		/* redundant; specified by -d */
	*long_outfile,		/* 'long', one-residue-per line outfile (-o) */
	*rasmol_file,		/* RasMol script file (-r) */
	*helix_register_file;	/* redundant; specified by -w */

