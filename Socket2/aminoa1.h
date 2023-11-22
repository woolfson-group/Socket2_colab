/*

					SOCKET
					 v3.02

					aminoa1.h

					24-10-01
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

/*					aminoa1.h
					---------
*/

/* JW 3-4-98
	defines enumerated type 'amino_acid', so that variables of this type can
	be used directly, e.g. 
		enum amino_acid aa1;
		aa1 = K;
	etc.
*/

/* last modified: 11-3-1
	initializations of amino_acid[1] and amino_acid[3] replaced by
	dummy definitions; see below;
	also, #define AMINO_ACIDS 21 is now done in preproc.h, so it need
	be declared only once (preproc.h is #include-d by statchar.c,
	and must also be #include-d before aminoa1.h is parsed)
*/

/* The amino acid codes in the enum type amino_acid can be specified in
	any order (you might want to change it so that they are grouped
	by property eg hydrophobicity). Type enum amino_acid is mapped
	from the corresponding ascii characters by the int array
	map_alpha_to_amino_acid[]. If you do change the order, make sure
	that amino_acid1[] and amino_acid3[] are in the same order...

	MAKE SURE THAT 'X' THE FIRST CODE THOUGH.

	The other thing is, make sure you have a line calling the function
	'aa_map()' in your code.

*/
enum amino_acid {X,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y};
enum amino_acid aa;
char *amino_acid1[AMINO_ACIDS]; /* dummy declaration, so that amino_acid1[]
				is known to all source files which #include global.h ;
				the values cannot be initialized here, because that
				must be done only once, or else the compiler gets upset;
				it is done in statchar.c */



char *amino_acid3[AMINO_ACIDS]; /* dummy declaration; see above */

int map_alpha_to_amino_acid[26];

char heterogen3[HETEROGENS_MAX][4]; /* dummy declaration; see above */

int map_heterogen_no_to_amino_acid[HETEROGENS_MAX]; /* dummy declaration; see above */

char *solvent_heterogen3[NON_AA_HETEROGENS]; /* dummy declaration; see above */


/* end */
