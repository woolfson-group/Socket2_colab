/*

					SOCKET
					 v3.03
					 
					 prototyp.h

					13-07-06
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

/*					prototyp.h
					----------
*/


/*	 			1 section:
			FUNCTION PROTOTYPES
*/

/* FUNCTION PROTOTYPES: ----------------------------------------------------- */

/* set_flagsNpars: this turns the command line arguments (held in argc and
	argv) into the appropriate settings in GLOBAL arrays par and setflag;
	the indices of these arrays are respectively members of the GLOBAL
	enumerated-type arrays pars and flags (the sizes of which are defined by
	constants PARS and FLAGS); the 'names' of these arguments (specified on
	the command line) which are interpreted are held in GLOBAL arrays
	parname and flagname; the values of par are set to strings, which might
	be subsequently read to provide integer or real variables if
	appropriate; the values of setflag are the GLOBAL enumerated-type
	boolean. This might be a long-winded way of interpreting command-line
	arguments, but see notes re pre-processor constants. */

int set_flagsNpars(int argc, char *argv[]); 


/* determine_centre_of_mass: processes the list of all atoms in alpha-helices
	(atom_index atoms have been read by read_helical_pdb() ),and the
	coordinates of the appropriate atoms for each residue are used to
	calculate the mean coordinates of the side chain, i.e. the values of
	GLOBAL refatoms2; which atoms are used is determined by -a and -i */

void determine_centre_of_mass(int atom_index);


/* determine_end: processes the list of alpha-helix residues (there are
	residue_index residues) and their atoms (there are atom_index atoms),
	to determine the values of GLOBAL refatom1, which is either the
	coordinate of the single terminal atom of the side chain, or the mean of
	the two	equivalent terminal atoms of the side chain */

void determine_end(int residue_index, int atom_index);


/* measure_helix_pair: applies measure_residue_pair() to all combinations of
	residues R1, R2 where R1 is on helix1 and R2 is on helix2; helix1 and
	helix2 are helix-keys */

void measure_helix_pair(int helix1, int helix2, int residue_index);


/* measure_residue_pair: measures three distances between pairs of atoms in
	residues res1 and res2 (residue-keys); these distances are between the
	reference atoms of the same type, i.e. 0, 1 and 2, where 0 is alpha-
	carbon by default, 1 is side chain 'end' and 2 is side chain 'middle';
	only (2) is used to determine contacts, and (0) and (1) are only
	calculated if -v or -l are set, because they are printed to the screen
	and not actually used for anything else; if -v or -l, and there is a
	contact, then the distance between each possible pair of the end atoms
	is also	calculated, again only for show. */

void measure_residue_pair(int res1, int res2);


/* measure_CA_distance: retruns the distance between the two type (0) reference
	atoms (which are alpha-carbons by default) of residues indentified by
	residue-keys res1 and res2; at the outset it was anticipated that
	knob-into-holes packing would be distinguished by a combination of
	CA-CA, centre-centre, and end-end distances, but CA-CA distances were
	found to be unnecessary (as were end-end); they are retained for
	purposes of output in case anyone wants to do something useful with
	them, such as characterizing in more detail the characteristics of
	true knobs-into-holes packing */

float measure_CA_distance(int res1, int res2);


/* measure_end_distance: a cosmetic function which measures and prints the
	distances between each pair of atoms on the 'ends' of residues res1 and
	res2 (residue-keys); e.g. if each side chain has 2 end-atoms defined in
	GLOBAL refatom1 array, then 4 distances are calculated; pretty much
	redundant */

void measure_end_distance(int res1, int res2);


/* measure_end_distance_B: also of little use, this returns the difference
	between the single 'end' coordinates of the side chains res1 and res2
	(residue keys); see notes for measure_CA_distance() */

float measure_end_distance_B(int res1, int res2);


/* measure_centre_distance: this function returns the distance between the
	centres (mean atomic coordinates, held in GLOBAL refatom2) of residues
	identified by residue-keys res1 and res2; this is the basis of defining
	side chain-side chain contacts, which in turn define knobs-into-holes
	packing */

float measure_centre_distance(int res1, int res2);


/* measure_k_end_h_CA: returns the mean distance between the end (GLOBAL
	refatom1B) of a knob (identified by knob-key knobid) side chain and each
	of the 4 alpha-carbons (GLOBAL refatom0) of the residues forming the
	hole into which it fits; this is how types 1 and 3 ('knobs across
	holes') are distinguished from types 2 and 4 ('knobs-into-holes'); the
	latter have the result less than or equal to the insertion-cutoff (set
	at 7.0Å) */

float measure_k_end_h_CA(int knobid);


/* distance: returns the distance between two coordinates, each represented as
	a 3-element float array; used by the measure_*_diatance() functions
	described above */

float distance(float coords1[], float coords2[]);


/* check_complementarity: after all the knobs fitting into holes have been
	identified (the knobs will have been assigned either type 1 or 2),
	all the knobs (the total number of knobs is held in GLOBAL knob_index)
	are checked against one another for complementarity, using functions
	complementary() to identify pairwise-, and then check_daisy_chain() to
	identify cyclic-complementarity; complementary knobs have their (GLOBAL)
	knobtype increased to 3 or 4 (or even 5 or 6; see n_knobtype[] array),
	and the complementary knob(s) (if any) of each knob are stored in
	compknob[] array. ** N.B. refer to check_daisy_chain() **
	The call to check_daisy_chain() specifies an order of 0, and a
	1-dimensional array, which is a single element of the 2-dimensional
	GLOBAL array daisy_chain[][], is passed to check_daisy_chain(); the
	1-dimensional array local to check_daisy_chain() is also called
	daisy_chain[] (but its 1-DIMENSIONAL!). The GLOBAL daisy_chain[][]
	array is simply a list (indexed by the first []) of daisy chains. The
	element passed is simply the next empty daisy chain in the list. The
	number of daisy chains is held by the GLOBAL daisy_chains, and the
	current value of daisy_chains represents the next (empty) slot. So the
	call from check_complementarity() always passes
	daisy_chain[daisy_chains] to check_daisy_chain(). If the result is
	negative (-1 or -2 returned), then the components (if any) of the failed
	chain are effectively wiped from daisy_chain[daisy_chains] (because the
	next call to check_daisy_chain() will re-use this slot); otherwise, the
	newly found daisy-chain is checked to see if its the same as any
	previously found daisy chains (earlier in the GLOBAL daisy_chain[][]
	list), but with the knobs simply listed in a different order; if it is
	indeed the same, then daisy_chain[daisy_chains] is wiped. Otherwise, it
	is a genuinely new daisy chain, and so daisy_chains is incremented. This
	continues until all knobs have been used as the argument for a call to
	check_daisy_chain() from check_complementarity(). */

void check_complementarity();


/* check_duplication: it is possible, but rare, for a long side chain to make
	contacts with 2	different sets of 4 side chain on different helices -
	that is , it fits into 2 different holes, one on each helix (e.g. in the
	1aik MMOL structure, Leu 565:N fits into 
	GLN 563:A, LEU 566:A, GLN 567:A, VAL 570:A
	and
	ILE 635:C, TYR 638:C, THR 639:C, ILE 642:C);
	565:N will have been recorded as two different knobs, and this function
	detects such instances (not much is done with them, but it is useful to
	tell the user what they are) */

void check_duplication();


/* complementary: returns the knob-key of the knob residue which is also
	the holeresno'th hole residue of the knob identified by knob-key
	knobno, or returns -1 if there is no such knob; holesresno can only
	be 0,1,2,3 but only values of 1 or 2 make sense (the sides of the
	hole) */

int complementary(int knobno, int holeresno);


/* packing_angle: returns the angle between these two vectors:
	a) C-alpha to C-beta of knobres
	b) C-alpha of holeres1 to C-alpha of holeres2
	- knobres, holeres1 and holres2 are the residue-keys;
                                      _ _     _  _            _     _
	the angle is given by arccos( a.b / (|a||b|) ), where a and b are the
	two vectors; the result is converted to degrees */

float packing_angle(int knobres, int holeres1, int holeres2);


/* add_contact: the (GLOBAL) contact[] array stores, for each residue
	(identified by residue-key), the residue-keys of all the residues in
	contact with it; this function updates the contact array for residue
	res1, by adding res2 to it; if res2 is the 5th or later contact, a
	warning is issued */
/* JW 13-7-6 ; added third argument, which is the distance separating the centres
        of volume of the two residues */

void add_contact(int res1, int res2, float distance);


/* report_kih: after all the contacts have been assigned, the cases of residues
	with 4 or more contacts are assigned as knobs; each is assigned the
	knob-type 1 or 2, as determined by measure_k_end_h_CA(); the knob
	residues are put into the array knob[], the residues of their
	corresponding holes are put in hole[] and the number of knobs is
	stored in GLOBAL knob_index */

void report_kih(int residue_index);


/* reset_contacts: initializes number of contacts of each residue
	(n_contacts[], indexed by residue-key) to zero; is called only once;
	should probably be done in function initialize() */

void reset_contacts();


/* best_kih: If there are > 4 residues in contact with a side chain, then find
	the best (most 'hole-like') group of four; the first group of 4 putative
	hole residues with a 3,1,3 spacing (ie X--XX--X) will be used as the
	hole; if no such group exists then any group where the second and third
	residues are consecutive will be used. In either case, the contact[j][]
	array will be rearranged so that the 4 residues constituting the best
	group will occupy the first 4 elements. This is because later functions
	which test the hole spacing pattern assume that there is an x,1,y
	spacing pattern (anything but the first 4 elements are ignored); groups
	which do not match this will be discarded */

void best_kih(int residue_index);


/* write_files: creates the 'long' and 'data' output files, specified by -o and
	-d respectively, using file handles long_outfile and data_outfile */

void write_files(int residue_index);


/* find_register: the longest and most complicated function, which assesses each
	true coiled coil (see determine_order() ) and assigns to it not only
	the register of each residue within knobs-into-holes packing regions,
	but also the orientation. The source of the function is heavily
	commented, but here's a summary. It begins by assessing, within each
	coiled coil, the orientation of each pair of helices which interact
	with each other by means of knobs-into-holes (i.e. parallel or
	antiparallel). If the pair constitutes a 2-stranded coiled coil, then
	this gives the orientation of the coiled coil; if there are other
	helices in the coiled coil, then all the pairwise orientations are
	OR'd (1 = antiparallel) to give the orientation of the coiled coil.
	For each pair of interacting helices (even in the context of a coiled
	coil of order > 2), the N- and C-terminal limits of the helix-helix
	interaction are determined, in terms of both the terminal knobs, and
	the terminal hole residues with which terminal knobs interact. The
	relative orientation of these termini are used to determine the
	orientation of the two helices (parallel or antiparallel); however a
	second method based on a crude approximation of the helix axes (angle
	between them is given by arccos( a.b / (|a||b|) ), where a and b are
	the vectors representing the axes) is used as a check (this is done
	first; the other method overrides it if they are different). The point
	is, that with very short helices which are quite steeply inclined, it is
	not always clear what the orientation is; such helices are unlikely to
	have knobs-into-holes packing, but htey can do with the aid of long
	side chains, and it is possible for the orientation, with respect to the
	N- or C- terminality of a knob and its respective hole, to be opposite
	to the orientation based totally on the geometry of the helix axes.

	Once the orientations of neighbouring helices in the coiled coils have
	been established, they are used, in conjunction with the relative
	sequence-position of the side of a hole which is the complementary knob,
	to explicitly assign a heptad register to each knob.
	E.g. in a 2-stranded antiparallel coiled coil, knob K fits into a hole
	whose sides are residues i and i + 1. Knob K' is complementary to knob
	K (i.e. K' fits into a hole, on of the sides of which is K). Suppose
	that K' is i (rather than i + 1). This means that K has register a,and
	K' in register d. The rules are quite straightforward but not described
	here. a and d knobs can be assigned explicitly; in coiled coils with > 2
	helices, e and g residues can be knobs-in-holes, and such cases can be
	explicitly assigned. After explicit assignments have been done, they are
	used to infer the register of non-knob residues.

	If a RasMol script has been requested, set-definitions of the true
	coiled coils are written to it, along with the lists of residues which
	have each of the 7 registers.
	*/

int find_register(int residue_index);


/* terminal_orientation: this returns the relative orientation (0 for parallel,
	1 for antiparallel) of the two helices whose helix-keys are helix1 and
	helix2. Each helix is represented by a fairly crude axis, from the
	reference atom type (0) (i.e. alpha-carbon) of the N-terminal residue to
	the equivalent of the C-terminal residue. The interhelical angle is
	therefore       _ _     _  _                   _     _
		arccos( a.b / (|a||b|) ) where vectors a and b represent the
	two helix axes. Called by find_register() */

int terminal_orientation(int helix1, int helix2, int residue_index);


/* orientation_of_helices: LOCAL to find_register(), a table of relative helix
	orientations is stored in the arrays orientation_first_helix[],
	orientation_second_helix[] and orientation[]. The number of rows in the
	table is stored in variable pairs, also LOCAL. The values in the first
	two are helix-keys, while the values in the second are 0 and 1. The
	relative orientation of two helices might be referred to several times,
	for example if they belong to complex assemblies in which they
	effectively belong to more than one coiled coil simultaneously.
	Therefore they are determined once and then stored in the table. This
	function simply retrieves the relative orientation from this table (the
	3 arrays are passed to it, along with the helix-keys of the two helices
	queried); called by find_register() */ 

int orientation_of_helices(int helix1, int helix2,
	int orientation[], int orientation_first_helix[],
	int orientation_second_helix[], int pairs);

/* set_orientation_of_helices: this function resets the relative orientation of
	helices identified by helix-keys helix1 and helix2 in the table (see
	orientation_of_helices() above). This happens if the secondary
	determination (relying on terminality of knobs and their holes)
	disagrees with the result from terminal_orientation() (above); refer to
	find_register(), which calls it */

void set_orientation_of_helices(int helix1, int helix2, int new_orientation,
	int orientation[], int orientation_first_helix[],
	int orientation_second_helix[], int pairs);


/* check_extremes_of_hole: given the knob-key of a knob, each of the four
	residues of the hole into which	it fits are checked; if any are the most
	extreme (N- or C-terminal) residues on that helix (the one with the
	hole, not the knob) -evaluated by comparing with the current values of
	the GLOBAL coiled_coil_begin[][][1] and	coiled_coil_end[][][1] arrays -
	then one (or both) of these arrays is updated appropriately;
	called by find_register() */

void check_extremes_of_hole(int knob, int c);


/* relative_register: returns the appropriate register (as a char), given a
	register and a sequence offset;
	e.g. relative_register('e',3) = 'a' */

char relative_register(char reg, int offset);


/* check_daisy_chain: this is called by check_complementarity(), and by itself.
	Its purpose is to check a knob to test whether it is in a 'daisy-chain',
	i.e. a cyclic arrangment of knobs-into-holes (the hallmark of coiled
	coils with > 2 strands).
	It works by checking EITHER of the two residues (which, is determined by
	the variable direction, which is either 1 or -1) which form the sides of
	the hole into which the query knob (identified by knob-key thisknob)
	fits, to see if the hole residue is itself a knob; if it is, then the
	function calls itself with this latter knob as the argument. Each
	recursive call lengthens a chain of complementarity by one residue; the
	length of this chain is stored in the variable order, so-called because,
	if indeed the chain is found to form a closed loop, then the length of
	the loop is the order of the coiled coil.
	The components of the lengthening chain are stored in the 1-dimensional
	array daisy_chain, whose elements are the knob-keys of the knobs in the
	chain. A chain is found to be closed if one of the two residues, which
	form the sides of the hole into which knob identified by knob-key
	thisknob fits, is the original knob, i.e. the FIRST one in the
	daisy_chain[] array (it does not matter which of the two sides is the
	knob, just as long as one of them is; i.e. direction is not relevant
	here).
	If all the knobs in the list (there are GLOBAL knob_index knobs) are
	checked, and the chain does not close, then the function returns a value
	of -1.
	** N.B. in complex arrangements, it is possible for a chain to close
	itself, but not head-to-tail; i.e. the first knob is at the end of a
	linear chain which feeds into a cyclic arrangement. If this happens,
	then the function returns -2.
	** N.B. if a chains which close themselves immediately are ignored; that
	is if the closed chain has a length of 2, it is a pairwise-complementary
	interaction. Pairwise-complementarity is checked separately, before
	cyclic-complementarity. It would be neater to check them all at once,
	but checks for daisy chains were implemented much later than the
	pairwise-check. But it could probably do with being tidied up to this
	end, though its unlikely it would make the algorithm significantly
	faster.
	** N.B. the first call of check_daisy_chain() is from
	check_complementarity(); this call specifies an order of 0, and the
	1-dimensional array, used as the daisy_chain[] array locally by
	check_daisy_chain(), IS A SINGLE ELEMENT OF A GLOBAL 2-DIMENSIONAL ARRAY
	- WHICH IS ALSO CALLED 'daisy_chain' - SORRY. See
	check_complementarity()	*/

int check_daisy_chain(int thisknob, int order, int daisy_chain[], int direction);


/* determine_order: this function uses the lists of complementary knobs
	(pairwise-complementary and daisy-chains) to produce an initial list of
	'coiled-coils'. For these purposes, a 'coiled-coil' is either a pair of
	helices which interact by means of pairwise-complementary holes, or
	groups of N helices whose residues contribute to N-order daisy chains;
	some of the pairs of helices might also belong to daisy chains, and so
	are designated as subsets of the latter; this produces a list of true
	coiled coils (which are not subsets of any other coiled coils, and which
	have more than one layer of complementary knobs, if they are 2-stranded
	coiled coils); the order of each coiled coil is determined */

int determine_order(/*int residue_index*/);

/* define_ras_coils: if a RasMol file is being created, this function defines
	the set of coiled coils (the union of all coiled coils which arent
	subsets of other coiled coils) */

void define_ras_coils();


/* function get_code: filename[] is input (the name of a user-specified file,
	which may or may not include a path specification) and code[] is
	effectively the result of this function, ie a 4-letter PDB-style code.
	Such a 4-letter code is returned only if the file name (ignoring the
	path) is of the form either: \d\w{3} of pdb\d\w{3} - otherwise, the full
	file-name (minus path) is returned, which is treated as a pseudo-PDB
	'code' */

void get_code(char filename[], char code[]);


/* substring: takes a substring of string parent, and puts it into string sub
	(surely there is a C function which does this??); called only by
	get_code() */

void substring(char parent[], int start, int length, char sub[]);


/* dumpknobs: prints details of the knob whose knob-key is aknob, or all knobs
	if aknob == -1; used only for debugging (if -debug is specified) */

void dumpknobs(int aknob);


/* initialize: this function is called once by main(), which passes argc and
	argv to it; initialize() calls set_flagsNpars(argc,argv) to convert the
	command-line arguments into the correct values of par and setflag;
	initialize() reads int and float variables from some of the par strings,
	otherwise these variables are set to their default values. The various
	GLOBAL arrays are set to null values */
	
void initialize(int argc, char *argv[]);


/* check_files: the files specified by the command line (-f and -s are
	mandatory) are opened, and the program dies if any of them cannot be;
	the 4-letter code (if there is one) in the PDB input file name (-f) is
	derived, by calling get_code() */

void check_files();


/* pre_parse_dssp: this function is called only if the helix-extension
	parameter (held in the GLOBAL extension) is > 0 (if unspecified by -e,
	it is set to 0); it reads in the positions, in the sequence(s) of the
	PDB chain(s), of the alpha-helices, from the DSSP file (helix-residues
	are specified by 'H'). The N-terminal position of each helix is then
	decreased by extension residues likewise, the helix grows by extension
	residues at the C-terminal (but see notes for prune_extended_helices()).
	This is so that the function read_helical_dssp() will know which are the
	residues to be treated as alpha-helical (even if they don't have an 'H').
	The function sets the values of the GLOBAL helix_index, and the GLOBAL
	arrays helix_start, helix_end and helix_chain; N.B. main() passes these
	GLOBAL arrays to the function, which uses the same names for the local
	reference; so they could just be omitted as function parameters */

void pre_parse_dssp(int *helix_index, int helix_start[], char helix_start_iCode[],
	int helix_end[], char helix_end_iCode[], char helix_chain[], int extend);

/* read_helical_dssp: reads in the alpha-helix residues from the DSSP file, and
	returns the total number of helical residues read in; if extension > 0,
	then the table consisting of GLOBAL helix_start, helix_end, helix_chain,
	which has been previously determined by pre_parse_dssp(), is used to
	indicate which residues should be read in; otherwise, this table will be
	undefined, and read_helical_dssp() simply reads only those residues
	which have an 'H' secondary-structure type, filling in the table as it
	proceeds. RasMol commands defining sets corresponding to complete
	helices are then issued if appropriate. */

int read_helical_dssp(int extend);


/* within_helix: this is called only by read_helical_dssp, and only if
	extend > 0; it uses the helix table (helix_start, helix_end,
	helix_chain, which correspond to GLOBAL arrays of the same name, and
	number of rows n_helices, to return the helix-key of the helix to which
	a residue belongs, or -1 if the residue belongs to no helix. The residue
	is specified by the resno, which is the value of the PDB 'resSeq' field
	(NOT a residue-key) and the chain, corresponding to PDB 'chainID' field.
	*/

int within_helix(int resno, char iCode, char chain, int helix_start[],
	char helix_start_iCode[], int helix_end[], char helix_end_iCode[],
	char helix_chain[], int n_helices);

/* prune_extended_helices: if extend > 0 (set by -e; see pre_parse_dssp() ) then
	after pre_parse_dssp() has been called to define the groups of residues
	which are alpha-helical, and read_helical_dssp() has been called to
	read in the amino acids which fall within these limits, then
	prune_extended_helices() is called to check that the limits specified by
	GLOBAL helix_start[], helix_end[], helix_chain[], correspond to residues
	which actually exist. E.g. if a helix starts at residue 2, and
	extension == 2, then the beginning of the helix would be residue 0 (but
	there isn't one). So, this function modifies the helix_start and
	helix_end values appropriately. */

void prune_extended_helices(int residue_index, int helix_index,
	int helix_start[], int helix_end[], char helix_chain[]);


/* read_helical_pdb: reads in the names (PDB 'name'), atom no's (PDB 'serial'),
	PDB 'chainID', and coordinates PDB 'x', 'y' and 'z', of each atom
	belonging to an alpha-helix residue; these are respectively read into
	GLOBALS atom_name, atom_no, coord[0..2], and the corresponding
	residue-key (GLOBAL atom_res) is determined; the function returns the
	total number of alpha-helix atoms read in; if there are any hydrogen
	atoms, they are not read unless -a was specified */

int read_helical_pdb();


/* find_knobs_and_holes: does what it says, firstly by calling
	measure_helix_pair() for each unique combination of helices to find all
	inter-helical residue-residue contacts; then adjusting any groups of > 4
	contacts per residue with best_kih(); then assigning a type of 1 or 2 to
	each of these knobs with report_kih(); then checking for pairwise- and
	cyclically- complementary knobs with check_complementarity(), and
	reporting duplicate knobs with check_duplication(); the list of knobs in
	holes is then displayed, and if a RasMol script file has been requested
	with -r, RasMol commands defining the set of knobs for each helix, and
	holes for each helix, are written to file. */


void find_knobs_and_holes(int residue_index);


/* null_refatom2: if atom with serial no atomno has all 3 coordinates == 9999.99
	then a true result is returned; otherwise false. The initialize()
	function sets all atomic coordinates in GLOBAL coord array to 9999.99 at
	the outset (an undefined null value, which might well be 0.00, would be
	no good of course). */

int null_refatom2(int atomno);
