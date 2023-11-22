/*

					SOCKET
					 v3.02

					register.c

					02-11-01
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

/*					register.c
					----------
*/


#include "socket.h"
#include <stdlib.h>
/* 6-11-00 */
/* Minor alteration to output to preempt details of assessing
coiled-coil orientation
Corrected a mistake in the RasMol echo statements (the coiled-coils
banner was printed once for each coiled coil, but is now printed
only once)
- and eliminated some more messages for each coiled coil which should
not appear if the -q option is used
*/
/* 2-11-00 */
/*
Eliminated unwanted redundant diagnostic message, which reported details
of knob #5
User-friendly echo statements added to RasMol script if requested,
stating the helices which comprise each coiled coil
*/
/* 8-2-00 */
/* Corrected bug in loop starting line 1068, which puts together the strings
	representing the hole helix, and the type, of each knob residue;
	previously, where a residue effectively exists as more than one
	knob, not necessarily the right knob-details were shown in these
	strings (this only affected the output, not any calculations) */

/* 21-9-99 */

/* A complete re-writing of the find_register function, using the geometry of
individual knob-into-hole interactions, instead of the sequential spacing
between knobs */

/* Global variables used:

	coiled_coil[MAX_COILED_COILS][6]
		- stores the helix id of all helices (up to 6)
		contributing to this coiled coil
	coiled_coil_tally[MAX_COILED_COILS]
		- the number of times this combination of
		helices is thrown up in the daisy chain list
	coiled_coil_helices[MAX_COILED_COILS]
		- the number of helices in each coil
	coiled_coil_subset[MAX_COILED_COILS]
		- the coiled coil (if any) of which the helices in
		each coiled coil are a subset
	coiled_coil_begin[MAX_COILED_COILS][6][2]
		- stores the id of the most N-terminal knob
		contributing to this coiled coil, for each helix
	coiled_coil_end[MAX_COILED_COILS][6][2]
		- stores the id of the most C-terminal knob
		contributing to this coiled coil, for each helix
	coiled_coil_orientation[MAX_COILED_COILS]
		- stores the overall orientation of each coiled coil; if all
		the helices are parallel, the coiled coil is parallel;
		otherwise its antiparallel
	coiled_coils - the number of different coiled coils
	helix_start[MAX_HELICES] - *PDB number* of a helix' N-terminal residue
	helix_end[MAX_HELICES] - *PDB number* of a helix' C-terminal residue
	helix_order[MAX_HELICES] - the number of helices in the coiled coil
	code,cutoff2,extend,helix_start[MAX_HELICES], helix_end[MAX_HELICES]
*/

int find_register(int residue_index)
	{
	int c,d,daisy,helix,helix1,helix2,helix_orientation,i,j,k,l,r,pairs,side,true_ccs;

	int orientation_first_helix[MAX_PAIRS], /*(ID of first helix) */
	    orientation_second_helix[MAX_PAIRS], /*(ID of second helix)*/
	    orientation[MAX_PAIRS] /*(orientation: 0 for parallel, 1 for antiparallel)*/
	;
	/* - all these indexed by serial number of the helix-helix interaction. */


	/* these variable store, for a given helix 'helix' in a given coiled coil 'c', *knob IDs*
	of the N-most terminal and C-most terminal knobs of helix 'helix'. */

	int c_h_begin[MAX_HELICES_PER_COIL], c_h_end[MAX_HELICES_PER_COIL];

	/* these variable store, for a given helix 'helix' in a given coiled coil 'c', *knob IDs*
	of the knobs (on all the other helices of the coiled coil) which are complementary to the
	the N-most terminal and C-most terminal knobs of helix 'helix'. This is a means of
	determining whether the current helix is parallel or antiparallel to each of the others */

	int c_h_begin_partner[MAX_HELICES_PER_COIL], c_h_end_partner[MAX_HELICES_PER_COIL];

	/* the above means of orientation-determination may contradict the first stab which is
	done by means of pseudo-helix-axes; the following variable holds the number of helix-helix
	pairs which have their orientation changed in this way; if > 0, the orientations of all
	the coiled coils (dependent on all the pair orientations of their constituent helices)
	are re-assigned */

	int reorientate;

	static char* orientation_name[2] = {"", "anti"};
	static char* core = {"da"};
	static char* flank = {"ge"};
	static char* roman[51] = {"","i","ii","iii","iv","v","vi","vii","viii","ix","x",
				  "xi","xii","xiii","xiv","xv","xvi","xvii","xviii","xix","xx",
				  "xxi","xxii","xxiii","xxiv","xxv","xxvi","xxvii","xxviii","xxix","xxx",
				  "xxxi","xxxii","xxxiii","xxxiv","xxxv","xxxvi","xxxvii","xxxviii","xxxix","xl",
				  "xli","xlii","xliii","xliv","xlv","xlvi","xlvii","xlviii","xlix","l"};

	char result[4][MAX_HELIX_LENGTH], previous_was_knob, current_is_knob;
	static char* result_name[4] = {"sequence","register","partner ","knobtype"};

	/* when the patterns of knobs in holes are displayed, the helices are labelled with
	capital letters, eg 'X','Y','Z'. If a coiled coil has two helices, they will be labelled
	X and Y; if 3, X, Y and Z; if 4 W, X, Y and Z; if 5, V, W, X, Y and Z. */
	static char* alphabase = {"  XXWVU"};

	/* stores length of each tad of current coiled-coil region of current helix of current
		coiled coil*/
	int tadlength[200];

	/* number of tads in the current coiled-coil region of current helix of current coiled coil */
	int tad_index;

	/* number of non-canonical interrupts in the current coiled-coil region of current helix of
		current coiled coil */
	int nonc_breaks;

	/* total length of the current coiled-coil region of current helix of
		current coiled coil */
	int region_length;

	/* total number of non-canonical interrupts in all the coiled coils */
	int total_nonc_breaks = 0;

	/* length of longest non-canonical interrupts in any of the coiled coils */
	int longest_nonc = 0;

	/* workspace strings used for appending non-string variables to strings */
	char tmpstr[15],tmpstr2[15];

	/* initialize the orientation arrays to nul */
	for (i = 0; i < MAX_PAIRS; i++)
		{
		orientation_first_helix[i] = -1;
		orientation_second_helix[i] = -1;
		orientation[i] = -1;
		}

	pairs = 0;

	true_ccs = 0;

	reorientate = 0;

	if (par[par_r] != NULL) fprintf(rasmol_file,"echo\necho *** COILED COILS ***\necho\n");

	for (c = 0; c < coiled_coils; c++)
		{/* REDUNDANT BRACKET */
		if ((coiled_coil_subset[c] == -1) && ((coiled_coil_helices[c] > 2) || (coiled_coil_tally[c] > 1)))
			{
			coiled_coil_orientation[c] = 0;
			if (!setflag[flag_q]) printf("\n\ncoiled coil %2d:\n",c);

			/* this represents a new putative 'true coiled coil' in that
			   its set of helices is not a subset of another coiled coils' */

			/* the following loop makes a first determination of the
			   orientations of all the pairs of helices within this coiled
			   coil */

			for (helix1 = 0; helix1 < coiled_coil_helices[c] - 1; helix1++)
				{/* REDUNDANT BRACKET */
				for (helix2 = 1; helix2 < coiled_coil_helices[c]; helix2++)
					if (helix1 != helix2)
						{
						/* has this pair already had its orientation done? */
						i = 0;
						while ((i < pairs) &&
							((orientation_first_helix[i] != coiled_coil[c][helix1]) ||
							 (orientation_second_helix[i] != coiled_coil[c][helix2])))
							i++;
						if ((orientation_first_helix[i] != coiled_coil[c][helix1]) ||
						    (orientation_second_helix[i] != coiled_coil[c][helix2]))
							/* it hasn't been done yet */
							{
							/* does it need to be done - ie are there any contacts between these
							two helices ? */
							k = 0;
							while (
								(k < knob_index) &&
								( ( (helix_no[knob[k]] != coiled_coil[c][helix1]) ||
										(helix_no[hole[k][0]] != coiled_coil[c][helix2]) )
								&&( (helix_no[knob[k]] != coiled_coil[c][helix2]) ||
										(helix_no[hole[k][0]] != coiled_coil[c][helix1]) ) )
							      )
								k++;

							if ( ( (helix_no[knob[k]] == coiled_coil[c][helix1]) &&
									(helix_no[hole[k][0]] == coiled_coil[c][helix2]) )
								||( (helix_no[knob[k]] == coiled_coil[c][helix2]) &&
									(helix_no[hole[k][0]] == coiled_coil[c][helix1]) ) )

								{
								/* yes, the two helices have at least one knob-in-hole interaction, so
								determine their relative orientation */
								if (setflag[flag_v] || setflag[flag_l])
									printf("helices %2d and %2d are in contact\n",
										coiled_coil[c][helix1],coiled_coil[c][helix2]);
								if (pairs == MAX_PAIRS)
									{
									printf("Maximum number of helix pairs (%d) exceeded\n",MAX_PAIRS);
									exit(1);
									}
								orientation_first_helix[pairs] = coiled_coil[c][helix1];
								orientation_second_helix[pairs] = coiled_coil[c][helix2];
								orientation[pairs] =
									terminal_orientation(coiled_coil[c][helix1],
										coiled_coil[c][helix2],residue_index);
								if (orientation[pairs]) coiled_coil_orientation[c] = orientation[pairs];
								if (!setflag[flag_q]) printf("%sparallel\n",orientation_name[orientation[pairs]]);
								pairs++;
								}
							}
						}
				} /* end of orientation-determining loop REDUNDANT BRACKET*/

			if (!setflag[flag_q]) printf("this coiled coil is %sparallel\n",orientation_name[coiled_coil_orientation[c]]);

			/* next stage is to go through all the helices which form this
			   coiled coil, and attempt to assign a register to all the
			   residues within it. Note that a separate register is
			   assigned to each residue *for each coiled coil to which it
			   belongs*; in the majority of structures, a helix will
			   belong to only one coiled coil, but in structures like
			   Arac (eg 2arc) and colicin 1A (1cii) there are helices
			   which belong to 2. In Stat3B, a *pair* of helices belongs to
			   2 different 3-stranded antiparallel coiled coils (1bg1); in
			   HIV gp41 (1aik), each pair of central core helices belongs to
			   2 different 3-stranded arrangements, and each helix belongs to
			   3 different 3-stranded coiled coils.

			   !! REMOVE The register is stored in array tad_register[MAX_RESIDUES][5],
			   indexed respectively by the residue serial number, and an
			   identifier of the coiled coil to which this register relates.
			   NB this is **NOT** the coiled coil serial number, because
			   there might be loads of them (most are subsets of other coiled
			   coils) . Which is the coiled coil in question is held in the
			   int array tad_register_cc[5]. So for example, if the first
			   'true' coiled coil to be processed is coiled coil 6, then
			   tad_register_cc[0] == 6 , and tad_register[r][0] means the
			   register of residue r in the context of coiled coil 6. !!REMOVE!!


			   The register is stored in array tad_register[MAX_RESIDUES][MAX_COILED_COILS],
			   indexed respectively by the residue serial number, and the serial
			   number of the coiled coil to which this register relates. */


			/* initialize the register of this coiled coil before processing
			   any of the helices */

			for (r = 0; r < residue_index; r++)
				tad_register[r][c] = ' ';


			/* initialize the coiled_coil_begin and coiled_coil_end arrays */
				
			for (helix = 0; helix < coiled_coil_helices[c]; helix++)
				for (i = 0; i < 2; i++)
					{
					coiled_coil_begin[c][helix][i] = -1;
					coiled_coil_end[c][helix][i] = -1;
					}


			/* helix indexes which helix of the coiled coil is being
			   processed (NOT the serial ID of the helix, which is
			   coiled_coil[c][helix] */

			for (helix = 0; helix < coiled_coil_helices[c]; helix++)
				{ /* NOT REDUNDANT BRACKET */
				/* cycle thru all the residues in this helix.
				if a residue is a knob, then confirm that it is
				in this coiled coil; ie, if the coiled coil is of
				order > 2, is this knob in a daisy chain of the
				same order, which belongs to this coiled coil?
				If so, it is a core knob (register = a or d) and its
				register can be assigned by the geometry of the hole
				into which it fits: one of the side residues must also
				be in the same daisy chain in this coiled coil, and
				whether this residue is sequentially before or after
				the other side residue, combined with the orientation
				of the two helices, specifies whether its a or d.

				If the coiled coil is of order > 2 but the order of
				this knob = 2 and the helix of the hole into which this
				knob fits is also in this coiled coil, then, with respect
				to this coiled coil, the knob is an e or g. NB the
				same pair of pairwise-complementary knobs can simultaneously
				belong to 2 different coiled coils, with 2 different
				registers, such as in Stat3B or HIV gp41. Which register
				is correct wrt this coiled coil is determined by which
				of the 2 side residues of the hole into which this knob
				fits is a member of a daisy-chain which belongs to this
				coiled coil: again, whether its before or after the other
				side residue, in conjunction with the orientation of the
				2 helices, determines the register.
				If neither of the two side residues is in a daisy chain of
				the same order as the coiled coil - this is possible if a
				3- or 4- stranded coiled coil 'frays' at the end, so that
				say, 2 of the helices continue to form KiH interactions with
				each other in the absence of the 3rd and/or 4th helix/es.
				In such a case, it is best not to try and assign a register
				to the knob. This is because, such a region can be considered
				in two different ways. The interactions can be treated as
				either part of the higher order coiled coil, but the cyclic
				interaction is not complete in this layer; or, as a region
				of 2-stranded coiled coil, in which case the register would
				be independent of the higher-order coiled coil. However, the
				latter option would only make sense if this region was before
				or after the complete 3- (or 4-) stranded section, not in the
				middle of it. Because it is expected that in some layers,
				not all expected knobs and holes are present, the former
				is the best course. An assignment could possibly be done by
				determining which of the two sides of the hole is not the
				complementary knob, but this won't really be necessary.

				If the coiled coil is of order 2, then the register is
				determined by which of the two side residues (of the hole
				into which the knob fits) is the complementary knob to
				this knob. Once again, whether its before or after the
				other side residue, combined with helix orientation, gives
				the answer. */

				if (setflag[flag_v] || setflag[flag_l]) printf("\thelix #%d (%d)\n",helix,coiled_coil[c][helix]);

				for (i = 0; i < coiled_coil_helices[c]; i++)
					{
					c_h_begin[i] = -1;
					c_h_end[i] = -1;
					}

				/* r cycles thru all the residues */

				for (r = 0; r < residue_index; r++)
					if (helix_no[r] == coiled_coil[c][helix])
						/* residue r is part of the helix'th helix
						   of coiled coil c */
						{
						if (setflag[flag_v]) printf("\t\tresidue #%4d (%4d:%c, iCode='%c')\n",r,helix_residue_no[r],
													helix_chain[helix_no[r]],helix_residue_iCode[r]);
						for (k = 0; k < knob_index; k++)
							if ((knobtype[k] > 2) && (knob[k] == r))
							    {
							    if (setflag[flag_v]) printf("\t\t\t= knob %d (type %d, order %d)\n",
												k,knobtype[k],knob_order[k]);
							    /* knob k is a proper knob (with at least
								one complementary knob); knob k is also
								residue r, and part of this
								helix (helix'th helix of coiled coil c);
								but is k part of coiled coil c? Not
								necessarily, so check:

								confirm that knob k actually fits into a hole which is
								a member helix of this coiled coil; its possible for it
								to be fitting into a hole on a helix outside this coiled coil,
								eg Stat3B */

							    helix2 = 0;
							    while ((helix2 < coiled_coil_helices[c]) &&
									(coiled_coil[c][helix2] != helix_no[hole[k][0]]))
								helix2++;

							    if (coiled_coil[c][helix2] == helix_no[hole[k][0]])
								{
								/* Knob k and its hole are both part of helices which are in this
								   coiled coil. There are two ways that the extremities of the helices
								   (held in arrays coiled_coil_begin[][][] and coiled_coil_end[][][], each
								   indexed respectively by 1, the coiled coil ID; 2, the coiled coil's
								   helix number -NOT helix ID and 3, the extremity definition- read on)
								   can be determined.

								   They can be considered either to be the most N-terminal and most
								   C-terminal knobs (which fit into a hole in a helix which belongs to
								   this coiled coil) of each helix; this is stored
								   in coiled_coil_begin[][][0] and coiled_coil_end[][][0], where the
								   values are *the serial numbers of the KNOBS* (NOT the residues)

								   or the most N-terminal and most C-terminal *residues* which are part
								   of holes into which the above knobs fit - these *residue IDs* are
								   held in the arrays coiled_coil_begin[][][1] and coiled_coil_end[][][1]
								   and are determined by the extremity() function.

								   For now, only the former are determined (ie the most extreme knobs) -
								   the latter can be easily determined later, when all the knobs have
								   had their register assigned.

								   A complication is that which knob is the most N- or C-terminal cannot
								   be determined from the knob serial numbers - knobs are not necessarily
								   catalogued in the sequential order (because they are found in the
								   context of helix-helix pairs). Their residue serial numbers
								   (knob[] array) must be used, which is easy enough. The current knob is
								   residue serial number r, of course. */
/*printf("checkpoint 1\n");*/
								if ((coiled_coil_begin[c][helix][0] == -1) ||
								   (r < knob[coiled_coil_begin[c][helix][0]]))
									coiled_coil_begin[c][helix][0] = k;
								if ((coiled_coil_end[c][helix][0] == -1) ||
								   (r > knob[coiled_coil_end[c][helix][0]]))
									coiled_coil_end[c][helix][0] = k;



								/* knob k is in this coiled coil (c) if k
								is part of a daisy chain which is part of
								coiled coil c and the order of the coiled
								coil is > 2; or if other order of c is 2,
								and both k's helix, and the helix of k's
								hole, belong to coiled coil c */

								if (coiled_coil_helices[c] == 2)
									{/* NOT REDUNDANT BRACKET */
									/* its a 2-stranded coiled coil
									- just check that the helix of
									this knob's complementary partner
									is the other helix in this coiled
									coil. The helix of this knob is
									coiled_coil[c][helix], and 'helix'
									must be 0 or 1. Therefore the other
									helix is either 1 or 0 in the
									coiled_coil[c] array

									a complication is that it is possible
									for a knob to have more than one
									complementary knob - so check them
									all */

									if (setflag[flag_v])
										printf("\t\t\tassigning on 2-stranded basis\n");

									for (i = 0; i < n_compknob[k]; i++)
										if (helix_no[knob[compknob[k][i]]] ==
											coiled_coil[c][1-helix])
											{
											/* complementary knob to k found -
											assign register */

											side = -1;
											for (j = 1; j < 3; j++)
												if (hole[k][j] == knob[compknob[k][i]])
													{
													/* the hole helix is the
													helix2'th helix in this coiled
													coil */

/*printf("checkpoint 2A\n");*/
/*printf("before:c_h_begin[%d] = %d; c_h_end[%d] = %d\n",helix2,c_h_begin[helix2],helix2,c_h_end[helix2]);*/
													if ((c_h_begin[helix2] == -1)
												|| (knob[k] < knob[c_h_begin[helix2]]))
														{
														c_h_begin[helix2] = k;
														c_h_begin_partner[helix2]
														 = compknob[k][i];
														}
													if ((c_h_begin[helix2] == -1)
												     || (knob[k] > knob[c_h_end[helix2]]))
														{
														c_h_end[helix2] = k;
														c_h_end_partner[helix2]
														 = compknob[k][i];
														}

													check_extremes_of_hole(k,c);
/*printf(" after:c_h_begin[%d] = %d; c_h_end[%d] = %d\n",helix2,c_h_begin[helix2],helix2,c_h_end[helix2]);*/
/*printf("checkpoint 3A\n");*/

													side = j - 1;
													helix_orientation =
orientation_of_helices(helix_no[knob[k]],helix_no[hole[k][j]],orientation,orientation_first_helix,orientation_second_helix,pairs);

											/* NOTE 1. The register of the knob is determined
											from its complementary partner's position in the
											hole ('side' variable) = 0 if it is sequentially
											the first of the two residues forming the sides,
											= 1 if is the second; and from the relative
											orientation of the knob and hole helices (variable
											'helix_orientation':

												     helix_orientation
													0	1
												0	d	a
											side
												1	a	d

											so and XOR of the two variables gives 0 for d,
											1 for a; character 0 of string variable 'core'
											is "d", character 1 is "a" */
											
										tad_register[r][c] = core[side ^ helix_orientation];

											if (setflag[flag_v] || setflag[flag_l])
printf("coiled coil %2d; helix %2d; residue %4d is knob %3d;\n\t\tcomplementary knob (#%d = knob %2d) is side %d of hole; helix orientation = %sparallel => register = %c\n",
						c,helix,r,k,i,compknob[k][i],side,orientation_name[helix_orientation],tad_register[r][c]);
													}

											/* A check. One of the two sides of the hole
											in which k fits should be the complementary
											knob */

											if (side < 0)	{
							printf("could not find a complementary knob (#%d) to knob %d\n",
												i,k);
												exit(1);
												}


											}
									
									} /* end of if (coiled_coil_helices[c] = 2) REDUNDANT*/
								else	{
									/* its a 3 or more-stranded coils coil,
									so check all the daisy chains to find
									one which is in coiled coil c, and
									includes knob k */

									/* daisy is the serial number of the daisy chain of which
									   knob k is part, if any */
									daisy = -1;

									/* only bother checking the daisy chains for knob k if
									k is of the same order as the coiled coil; if its of order 2
									then it won't be in the list. Note that its possible for k to
									have the same order as the coiled coil, but not be in a daisy
									chain *in this coiled coil* - ie it participates in a cyclic
									arrangement outside of this coiled coil, but makes a pairwise
									interaction in this coiled coil */

									/* version 2.02: note also that its possible for a knob to
									have a lower order than the coiled coil, but a higher order than
									2, and to be in this coiled coil; in which case the same rules
									apply as if it had the same order as the coiled coil. Only
									known example of this so far is the 3-membered layer in the
									4-stranded 1sfc COMP coiled coil. So, the important thing is
									not that the order of the knob is the same as the order of the
									coiled coil; the important thing is that the order of the knob
									is not 2 */

									if (knob_order[k] > 2)
									    { /* NOT REDUNDANT BRACKET */

									    if (setflag[flag_v])
				printf("\t\t\tsame order as coiled coil - looking for daisy chains of which this knob is a member\n");

									    for (d = 0; d < daisy_chains; d++)
										if (daisy_chain_cc[d] == c)
											{
											/* daisy chain d is part of
											   coiled coil c 

											   i indexes the knobs in daisy
											   chain d*/
											i = 0;
											while ((daisy_chain[d][i] != -1)
												&& (daisy_chain[d][i] != k)
												&& (i < MAX_DAISIES))
												i++;

											if ((daisy_chain[d][i] == k) && (i < MAX_DAISIES))
												{
												/* knob k is a member of
												   daisy chain d */

												if (setflag[flag_v])
									printf("\t\t\t\t- knob %d is a member of daisy chain %d\n",k,d);

												/* one, and only one, of the
												two residues forming the sides
												of the hole into which k fits
												should also be in this
												daisy chain - which is it?

												side == 0 => the knob complementary
													to k is sequentially before
													the other side in hole[k],
													ie x+3 in the diagram below

													    x

													x+4	x+3

													   x+7


												side == 1 => the knob complementary
													to k is sequentially after
													the other side in hole[k]
													ie x+4 in the diagram above

												NB hole[0] == x; hole[1] == x+3;
												   hole[2] == x+4; hole[3] == x+7.
												x is the position in the amino acid
												sequence */


												side = -1;

												/* j indexes the hole residues */
/* make the other similar j loops into whiles */
												for (j = 1;
												/*while ((*/j < 3; j++)
/*) && (side == -1))*/
													{
													i = 0;
													while ((daisy_chain[d][i] != -1)
													   && (knob[daisy_chain[d][i]]
														 != hole[k][j])
													   && (i < MAX_DAISIES))
														i++;
/*printf("!!! daisy_chain[%d][%d] = %d\n",d,i,daisy_chain[d][i]);*/
													if (knob[daisy_chain[d][i]]
														 == hole[k][j])

														/* the jth hole residue
														is the complementary
														knob to k */
														{
														/* the hole helix is the
														helix2'th helix in this
														coiled coil */

/*printf("checkpoint 2B\n");*/
/*printf("before:c_h_begin[%d] = %d; c_h_end[%d] = %d\n",helix2,c_h_begin[helix2],helix2,c_h_end[helix2]);*/
													if ((c_h_begin[helix2] == -1) ||
												(knob[k] < knob[c_h_begin[helix2]]))
															{
														c_h_begin[helix2] = k;
											c_h_begin_partner[helix2] = daisy_chain[d][i];
															}
													if ((c_h_end[helix2] == -1) ||
												(knob[k] > knob[c_h_end[helix2]]))
															{
														c_h_end[helix2] = k;
														c_h_end_partner[helix2]
														 = daisy_chain[d][i];
															}
													check_extremes_of_hole(k,c);
/*printf(" after:c_h_begin[%d] = %d; c_h_end[%d] = %d\n",helix2,c_h_begin[helix2],helix2,c_h_end[helix2]);*/
/*printf("checkpoint 3B\n");*/

														side = j - 1;
														helix_orientation =
orientation_of_helices(helix_no[knob[k]],helix_no[hole[k][j]],orientation,orientation_first_helix,orientation_second_helix,pairs);

														/* see NOTE 1 for
														an explanation */

										tad_register[r][c] = core[side ^ helix_orientation];

												if (setflag[flag_v] || setflag[flag_l])
printf("coiled coil %2d; helix %2d; residue %4d is knob %3d;\n\t\tcomplementary knob %d is side %d of hole; helix orientation = %sparallel => register = %c\n",
					c,helix,r,k,daisy_chain[d][i],side,orientation_name[helix_orientation],tad_register[r][c]);
														}

													/*j++;*/
													}

												/* A check. Unless something is very
												wrong, one of the two sides of the hole
												in which k fits MUST be a knob in the
												same daisy chain */

												if (side < 0)	{
							printf("could not find a complementary knob to knob %d in daisy chain %d\n",
													k,d);
													exit(1);
													}

												daisy = d;



												} /* if (daisy_chain[d][i] == k) */

											} /* end of if (daisy_chain_cc[d] == c) */

									    if ((daisy == -1) && (setflag[flag_v]))
printf("\t\t\t - knob is of same order as coiled coil, but does not belong to any daisy chains constituting this coiled coil\n");

									    } /* end of if (knob_order[k] == coiled_coil_helices[c])
										NOT REDUNDANT BRACKET */


									else
									    {
									    /* knob k is not of the same order as the coiled coil c */
									    if (setflag[flag_v])
									printf("\t\t\t - not the same order as the coiled coil\n");
									    }

								    	if ((daisy == -1) || (knob_order[k] == 2))
									    {
									    /* knob k is either:
										of the same order as coiled coil c, but
										its participation in c is of a lower order (ie its
										higher-order nature is due to its belonging
										simultaneously to a different coiled coil)
										- assuming it is of order 2 *with respect to coiled
										coil c* then it is probably an e or g knob
										- for this to be true, its complementary knob should
										be of the same order as coiled coil c, and in a daisy
										chain, of the same order, which belongs to c

									    or:
										not of the same order as the coiled coil c */


									    if (setflag[flag_v])
		printf("\t\t\t\tlooking for complementary knobs which are members of daisy chains constituting this coiled coil\n");

									    side = -1;

									    for (l = 0; l < n_compknob[k]; l++)
										{
										/* check that k's lth compknob is in a daisy
										chain in this coiled coil */
									    	for (d = 0; d < daisy_chains; d++)
											if (daisy_chain_cc[d] == c)
												{
												/* daisy chain d is part of
												   coiled coil c 

										  		 i indexes the knobs in daisy
												   chain d*/
												i = 0;
												while ((daisy_chain[d][i] != -1)
												&& (daisy_chain[d][i] != 
													compknob[k][l])
												&& (i < MAX_DAISIES))
													i++;

												if (daisy_chain[d][i] == 
													compknob[k][l])
													{

													/* the l'th complementary
													knob of knob k is indeed
													in a daisy chain in this
													coiled coil, ie a *core*
													residue;
													assign register */

													if (setflag[flag_v])
			printf("\t\t\t\t\tcomplementary knob #%d ( = knob %d) of knob %d is in daisy chain %d\n", l,compknob[k][l],k,d);

													daisy = d;

													for (j = 1; j < 3; j++)
													if (hole[k][j] ==
													    knob[compknob[k][l]])
														{
														/* the hole helix is the
														helix2'th helix in this
														coiled coil */

/*printf("checkpoint 2C\n");*/
/*printf("before:c_h_begin[%d] = %d; c_h_end[%d] = %d\n",helix2,c_h_begin[helix2],helix2,c_h_end[helix2]);*/
													if ((c_h_begin[helix2] == -1)
												   || (knob[k] < knob[c_h_begin[helix2]]))
														{
														c_h_begin[helix2] = k;
														c_h_begin_partner[helix2]
														 = compknob[k][l];
														}
													if ((c_h_end[helix2] == -1)
												   || (knob[k] > knob[c_h_end[helix2]]))
														{
														c_h_end[helix2] = k;
														c_h_end_partner[helix2]
															= compknob[k][l];
														}
/*printf(" after:c_h_begin[%d] = %d; c_h_end[%d] = %d\n",helix2,c_h_begin[helix2],helix2,c_h_end[helix2]);*/
/*printf("checkpoint 3C\n");*/

														side = j - 1;
														helix_orientation =
orientation_of_helices(helix_no[knob[k]],helix_no[hole[k][j]],orientation,orientation_first_helix,orientation_second_helix,pairs);

														/* see NOTE 1 for
														an explanation, but
														substitute "g" for "d"
														and "e" for "a" - the
														string variable 'flank'
														= "ge" */

										tad_register[r][c] = flank[side ^ helix_orientation];
												if (setflag[flag_v] || setflag[flag_l])
	printf("coiled coil %2d; helix %2d; residue %4d is knob %3d;\n\t\tcomplementary knob (#%d = knob %2d) is side %d of hole; helix orientation = %sparallel => register = %c\n",
					c,helix,r,k,l,compknob[k][l],side,orientation_name[helix_orientation],tad_register[r][c]);
														}

											/* A check. One of the two sides of the hole
											in which k fits should be the complementary
											knob */

													if ((side < 0)
												&& (l == n_compknob[k] -1))
														{
						printf("could not find a complementary knob (#%d = knob %d) to knob %d\n",
														l,compknob[k][l],k);
														exit(1);
														}


													} /* end of if
														(daisy_chain[d][i]
													     == compknob[k][l]) */

												} /* end of if (daisy_chain_cc[d] == c)*/

										if (daisy == -1)
											{
											if (setflag[flag_v])
	printf("\t\t\t\t\tcomplementary knob #%d ( = knob %d) of knob %d is not in any daisy chains in this coiled coil\n",
											l,compknob[k][l],k);
											/* the situation here is that the knob k and
											its complementary knob form a *pairwise*
											interaction, but both are members of helices
											which belong to this coiled coil, which has
											a higher order (> 2). One option is to assign
											the register based on this pairwise interaction,
											as if it was part of a 2-stranded coiled coil
											(which arguably it is- perhaps at the end of
											a 3- or 4-stranded coiled coil, where the other
											helix or helices have 'fallen away' -like in
											Stat3B 1bg1 or colicin 1A 1cii). This will
											almost certainly result in the register of this
											part of the helices being out of synch with the
											rest of it.

											Or (as is currently actually done) - just leave it
											*/

											/* Whatever, a note needs to be made of this
											interaction in terms of measuring the extent of
											KiH packing along the length of the helix */

											for (i = 0; i < n_compknob[k]; i++)
												for (j = 1; j < 3; j++)
													if (hole[k][j] ==
														knob[compknob[k][i]])
													/* if the register *were* to be
													assigned based on this pairwise
													interaction, now is the time: by
													setting side to j - 1 */
													{
/*printf("checkpoint 2D\n");*/
/*printf("before:c_h_begin[%d] = %d; c_h_end[%d] = %d\n",helix2,c_h_begin[helix2],helix2,c_h_end[helix2]);*/
		 											if ((c_h_begin[helix2] == -1)
										   		|| (knob[k] < knob[c_h_begin[helix2]]))
														{
														c_h_begin[helix2] = k;
														c_h_begin_partner[helix2]
														 = compknob[k][i];
														}
													if ((c_h_end[helix2] == -1)
										   		|| (knob[k] > knob[c_h_end[helix2]]))
														{
														c_h_end[helix2] = k;
														c_h_end_partner[helix2]
															= compknob[k][i];
														}
													check_extremes_of_hole(k,c);
/*printf(" after:c_h_begin[%d] = %d; c_h_end[%d] = %d\n",helix2,c_h_begin[helix2],helix2,c_h_end[helix2]);*/
/*printf("checkpoint 3D\n");*/
													}


											}

										} /* end of l for loop */

									    } /* end of if ((daisy == -1) ||
											(knob_order[k] != coiled_coil_helices[c])) */

									} /* end of if (coiled_coil_helices[c] = 2) ELSE */

								} /* end of if (coiled_coil[c][helix2] == helix_no[hole[k][0]]) */
							    else
								{
								if (setflag[flag_v])
	printf("\t\t\t\tknob %d fits into a hole in a helix (%d) which is not part of this coiled coil\n",k,helix_no[hole[k][0]]);
								}

							    } /* end of if ((knobtype[k] > 2) && (knob[k] == r)) */
						} /* end of if (helix_no[r] == coiled_coil[c][helix]) */

				/* do a double check on the relative orientations of the helices (parallel or antiparallel).
				Necessarily, the residue which is the most C-terminal (end) knob will
				have a higher residue serial number than the most N-terminal (because residue serial numbers will be
				in the same order as the residue sequence numbers in the PDB file, because PDB files list residues in
				sequential order).

				Suppose that the helix = 0, and the 0th helix in coiled coil c is helix 10. Say that the next helix
				in c is helix 12. That is, coiled_coil[c][0] = 10, coiled_coil[c][1] = 12.
				Of all the knobs of helix 10 *that fit into a hole in helix 12*, suppose the most N-terminal is knob
				5, and the most C-terminal is knob 8. Say knob 5 is residue 100, and that knob 8 is residue 115.
				Now consider helix 12's knobs *which are complementary* to these two knobs. Say the complementary knob
				(on helix 12) to knob 5 is knob 9; and that the complement (on helix 12) of knob 8 is knob 14. Let's say
				that currently we are evaluating the 0th helix of coiled coil c, ie helix 10 (that is, the loop variable
				helix would currently equal 10).
				All this would be stored as:
				c_h_begin[1] = 5;  c_h_end[1] = 8;  c_h_begin_partner[1] = 9;  c_h_end_partner[1] = 14
				Now consider the residue serial numbers of knobs 9 and 14. Say knob 9 is residue 150 and knob 14 is
				residue 164. That is, the equivalent position on helix 12 to the first (most N-terminal) knob of
				helix 10 is residue 150, and the equivalent position on helix 12 to the last knob of helix 10 is
				residue 164. Conclusion: helix 12 goes in the same direction as helix 10, ie they are antiparallel.

				This second determination of orientation is worthwhile, because the inital attempt is done by vectors of
				alpha-carbon positions, which will usually work; but in cases of extreme packing angles between helices,
				this could fall down, especially if the helices are short. This would be rare in helices with knobs-into
				-holes interactions, but it is worth getting it right.

				Back to the example; imagine that coiled coil c has 3-strands, and the 3rd helix is helix 13, ie
				coiled_coil[c][2] = 13. Of the knobs on the current helix (helix 10), the first which fits into helix 13
				is knob 15 and the last is knob 18. Say that these knobs are respectively residues 103 and 112; and that
				the complementary knobs, on helix 13, to knobs 15 and 18 are knobs 21 and 20, respectively:
				c_h_begin[2] = 15;  c_h_end[2] = 18;  c_h_begin_partner[2] = 21;  c_h_end_partner[2] = 20
				Knob 21 turns out to be, say, residue 190 and knob 20 is residue 200.
				This would mean that helices 10 and 13 are antiparallel.
				*/

				for (helix2 = 0; helix2 < coiled_coil_helices[c]; helix2++)
					if (helix2 != helix)
						{
						if (setflag[flag_v])
printf("secondary evaluation of orientation of helix #%d (%d) v helix #%d (%d):\n\t#%d's complement to #%d's first (knob %d) is knob %d (residue %d = %d:%c, iCode='%c');\n\t#%d's complement to #%d's last (knob %d) is knob %d (residue %d = %d:%c, iCode='%c')\n",
						helix,
						coiled_coil[c][helix],
						helix2,
						coiled_coil[c][helix2],
						helix2,
						helix,
						c_h_begin[helix2],
						c_h_begin_partner[helix2],
						knob[c_h_begin_partner[helix2]],
						helix_residue_no[knob[c_h_begin_partner[helix2]]],
						helix_chain[helix_no[knob[c_h_begin_partner[helix2]]]],
						helix_residue_iCode[knob[c_h_begin_partner[helix2]]],
						helix2,
						helix,
						c_h_end[helix2],
						c_h_end_partner[helix2],
						knob[c_h_end_partner[helix2]],
						helix_residue_no[knob[c_h_end_partner[helix2]]],
						helix_chain[helix_no[knob[c_h_end_partner[helix2]]]],
						helix_residue_iCode[knob[c_h_end_partner[helix2]]]);

						/* if there is only one complementary knob on helix #helix2, then the most N-terminal
						and the most C-terminal are the same, so their diference cannot be used to determine
						orientation */
					
						if (knob[c_h_begin[helix2]] == knob[c_h_end[helix2]])
							{
							if (setflag[flag_v])
								printf("\t\tthe two are the same; orientation cannot be reevaluated\n");
							}
						else	{
							if (knob[c_h_begin_partner[helix2]] < knob[c_h_end_partner[helix2]])
								helix_orientation = 0;
							else	helix_orientation = 1;

							/* now compare this with the primary evaluation */
							if (helix_orientation == orientation_of_helices(coiled_coil[c][helix],
								coiled_coil[c][helix2],orientation,orientation_first_helix,
								orientation_second_helix,pairs))
								{
								if (setflag[flag_v])
									printf("\t\tresult %sparallel, agrees with primary evaluation\n",
										orientation_name[helix_orientation]);
								}
							else	{
								reorientate++;
								if (!setflag[flag_v])
printf("secondary evaluation of orientation of helix #%d (%d) v helix #%d (%d):\n\t#%d's complement to #%d's first (knob %d) is knob %d (residue %d = %d:%c, iCode='%c');\n\t#%d's complement to #%d's last (knob %d) is knob %d (residue %d = %d:%c, iCode='%c')\n",
								helix,
								coiled_coil[c][helix],
								helix2,coiled_coil[c][helix2],
								helix2,
								helix,
								c_h_begin[helix2],
								c_h_begin_partner[helix2],
								knob[c_h_begin_partner[helix2]],
								helix_residue_no[knob[c_h_begin_partner[helix2]]],
								helix_chain[helix_no[knob[c_h_begin_partner[helix2]]]],
								helix_residue_iCode[knob[c_h_begin_partner[helix2]]],
								helix2,
								helix,
								c_h_end[helix2],
								c_h_end_partner[helix2],
								knob[c_h_end_partner[helix2]],
								helix_residue_no[knob[c_h_end_partner[helix2]]],
								helix_chain[helix_no[knob[c_h_end_partner[helix2]]]],
								helix_residue_iCode[knob[c_h_end_partner[helix2]]]);

								printf("!!!!!!\t\tresult %sparallel, disagrees with primary evaluation!!!!!!\n\t\t\t- resetting orientation to %sparallel",
								orientation_name[helix_orientation],orientation_name[helix_orientation]);

								set_orientation_of_helices(coiled_coil[c][helix],
								coiled_coil[c][helix2],helix_orientation,orientation,
								orientation_first_helix,orientation_second_helix,pairs);
								}
							}
						
						}

				} /* end of helix for loop NOT REDUNDANT BRACKET*/

			/* 13-4-00
			This is a SECOND loop to determine the orientation (parallel or antiparallel)
			of each coiled coil. This is only re-done if there were any discrepancies in
			the helix-helix pair orientations (flagged by the variable reorientate)
			*/
			
			if (reorientate)
				{
				if (!setflag[flag_q])
				printf("!!!!!! %d helix pairs in coiled coils have had their orientation reassigned\n",
					reorientate);
				coiled_coil_orientation[c] = 0;
				for (helix1 = 0; helix1 < coiled_coil_helices[c] - 1; helix1++)
					{/* REDUNDANT BRACKET */
					for (helix2 = 1; helix2 < coiled_coil_helices[c]; helix2++)
						if (helix1 != helix2)
							{
							i = 0;
							while ((i < pairs) &&
								((orientation_first_helix[i] != coiled_coil[c][helix1]) ||
								 (orientation_second_helix[i] != coiled_coil[c][helix2])))
								i++;
							if ((orientation_first_helix[i] == coiled_coil[c][helix1]) &&
							    (orientation_second_helix[i] == coiled_coil[c][helix2]))
								/* this pair of helices are in contact */
								{
								if (orientation[i]) coiled_coil_orientation[c] = orientation[i];
								}
							}
					} /* end of SECOND orientation-determining loop REDUNDANT BRACKET*/
				if (!setflag[flag_q])
					printf("coiled coil %d now assigned as %sparallel\n",c,orientation_name[coiled_coil_orientation[c]]);
				}



			true_ccs++;

			coiled_coil_max_length[c] = 0;
			coiled_coil_mean_length[c] = 0.0;

			for (helix = 0; helix < coiled_coil_helices[c]; helix++)

				{

				/*
				The extremities of the knob-into-hole packing region of each helix within the context of this coiled
				coil (c) have already been determined: they are held in arrays coiled_coil_begin[c][helix][0] and
				coiled_coil_end[c][helix][0]. The values are the *knob IDs* of the most N- and C-terminal knobs in
				the helix.

				The more liberal extremes, defined by the most N- or C-terminal
				residues of each *hole* involved in coiled coil c will have mostly been determined already. These are
				basically the most extreme residues of the holes into which the knobs in coiled_coil_end[c][helix][0]
				fit. They are not exactly the same because in higher-order coiled coils, core knobs do not generally
				have pairwise-complementary knobs on another helix. That is, A knob fits in B hole, B knob fits in
				C hole, C knob fits in A hole in a 3-stranded. Say A is on helix 1, B on 2, C on 3. Knob A might be the
				most extreme knob on 1, but C might not be the most extreme knob on 3. It is possible that no hole on A
				corresponds to the most extreme knob on another helix. The way to get round this is to check each and
				every hole residue corresponding to any knobs on any helix in the coiled coil. This has been done during
				the loop which attempts to assign register to the knobs. 

				The 'hole' extremities are stored in arrays
				coiled_coil_begin[c][helix][1] and coiled_coil_end[c][helix][1]. */

				/* In case a 'hole' extremity has not been set, or a 'knob' extremity is more N- or C - terminal
				than a 'hole' extremity has not been set: */

				if ((coiled_coil_begin[c][helix][1] == -1) ||
					(knob[coiled_coil_begin[c][helix][0]] < coiled_coil_begin[c][helix][1]))
					coiled_coil_begin[c][helix][1] = knob[coiled_coil_begin[c][helix][0]];
				if ((coiled_coil_end[c][helix][1] == -1) ||
					(knob[coiled_coil_end[c][helix][0]] > coiled_coil_end[c][helix][1]))
					coiled_coil_end[c][helix][1] = knob[coiled_coil_end[c][helix][0]];

				if ((coiled_coil_end[c][helix][1] - coiled_coil_begin[c][helix][1]) > (coiled_coil_max_length[c] - 1))
					coiled_coil_max_length[c] = coiled_coil_end[c][helix][1] - coiled_coil_begin[c][helix][1] + 1;
				coiled_coil_mean_length[c] += coiled_coil_end[c][helix][1] - coiled_coil_begin[c][helix][1] + 1;
				} /* end of helix for loop */


			coiled_coil_mean_length[c] /= coiled_coil_helices[c];



			printf("\n\n%s%5.1f %1d coiled coil (%s) %d (%sparallel %d-stranded, length max %d mean %5.2f):\n",
				code,cutoff2,extend,roman[true_ccs],c,orientation_name[coiled_coil_orientation[c]]
				,coiled_coil_helices[c],coiled_coil_max_length[c],coiled_coil_mean_length[c]);

			/* define this coiled coil in the rasmol script file, if one has been requested */

			if (par[par_r] != NULL) fprintf(rasmol_file,"define coiled_coil%d",c);



			for (helix = 0; helix < coiled_coil_helices[c]; helix++)

				{


				/* define the stretch of this helix which contributes to this coiled coil,
				   in the rasmol script file if appropriate */

				if (par[par_r] != NULL)
					{
					if (helix) fprintf(rasmol_file,",");

					/* N.B. when writing the RasMol script file, insertion codes (iCode)
					are assumed to be null; not sure if RasMol deals with them correctly */

					fprintf(rasmol_file," %d-%d:%c",
					helix_residue_no[coiled_coil_begin[c][helix][1]],
					helix_residue_no[coiled_coil_end[c][helix][1]],
					helix_chain[coiled_coil[c][helix]]);
					if (helix == coiled_coil_helices[c] - 1) fprintf(rasmol_file,"\n");
					}


				/* reset the results strings */
				for (i = 0; i < 4; i++) strcpy(result[i],"");

				printf("\n\nassigning heptad to helix %d (%c) %d-%d:%c\n",coiled_coil[c][helix],
						alphabase[coiled_coil_helices[c]] + helix, helix_start[coiled_coil[c][helix]],
						helix_end[coiled_coil[c][helix]],helix_chain[coiled_coil[c][helix]]);
				printf("extent of coiled coil packing: %3d",
					helix_residue_no[coiled_coil_begin[c][helix][1]]);
				if (helix_residue_iCode[coiled_coil_begin[c][helix][1]] != ' ')
					printf("'%c'",helix_residue_iCode[coiled_coil_begin[c][helix][1]]);
				printf("-%3d",helix_residue_no[coiled_coil_end[c][helix][1]]);
				if (helix_residue_iCode[coiled_coil_end[c][helix][1]] != ' ')
					printf("'%c'",helix_residue_iCode[coiled_coil_end[c][helix][1]]);
				printf(":%c\n",helix_chain[coiled_coil[c][helix]]);

				previous_was_knob = ' ';
				current_is_knob = ' ';
				tad_index = -1;
				nonc_breaks = 0;
				region_length = 0;

				/* r cycles thru all the residues */

				for (r = 0; r < residue_index; r++)
					if (helix_no[r] == coiled_coil[c][helix])
						/* residue r is part of the helix'th helix
						   of coiled coil c */
						{

						/* fill in the gaps in the register assignment
						*/

						current_is_knob = tad_register[r][c];

						if ((r >= coiled_coil_begin[c][helix][1]) && (r <= coiled_coil_end[c][helix][1])
							&& (tad_register[r][c] == ' '))
							/* this residue currently has a blank register assignment but is in the coiled
							   coil */
							{
							if (previous_was_knob != ' ')
								tad_register[r][c] = relative_register(previous_was_knob,1);
							else if ((helix_no[r-1] == helix_no[r]) &&
									(tad_register[r-1][c] >= 'a') && (tad_register[r-1][c] <= 'g'))
								tad_register[r][c] = relative_register(tad_register[r-1][c],1);
							else	/* it must be before the first knob */
								{
								i = r + 1;
								while ((tad_register[i][c] == ' ') && (i < residue_index)) 
									i++;

								if (helix_no[i] != helix_no[r])
									{ printf("couldn't find first assigned knob\n"); exit(1);}
								tad_register[r][c] = 
								relative_register(tad_register[i][c],r - i);
								}
							}

						previous_was_knob = current_is_knob;

						strcat(result[0],amino_acid1[helix_residue_aacode[r]]);
						sprintf(tmpstr,"%c",tad_register[r][c]);
						strcat(result[1],tmpstr);

/* at this point, a register assignment has been made to residue r, if r is
   inside the coiled-coil region of the current helix of the current coiled
   coil (the helix'th helix of coiled coil c
*/
						
						if (tad_register[r][c] != ' ')
							{
							region_length++;
							if ((tad_register[r-1][c] == ' ') ||
								(tad_register[r][c] != tad_register[r-1][c] + 1))

								/* its the start of a new tad */
								{
								tadlength[++tad_index] = 0;

								/* the first tad can be < 7 residues in a canonical heptad
								i.e. if it's incomplete */

								if ((tad_index > 1) && (tadlength[tad_index-1] != 7))
									{

									/* its a non-canonical break */									
									nonc_breaks++;
									}
								}
							tadlength[tad_index]++;
							}

						i = 0; /* i here marks the highest knobtype (in any helix) which this residue has */
						l = 0; /* l marks the highest knobtype which this residue has in this coiled coil
							(its possible for a side chain to be > 1 knob, in holes in different
							helices; these might, in fact usually will, be in different coiled coils) */

						

						for (k = 0; k < knob_index; k++)
							if (knob[k] == r)
								{
								if (knobtype[k] > i) i = knobtype[k];
								helix1 = helix_no[hole[k][0]];
								helix2 = 0;
								while ((coiled_coil[c][helix2] != helix1) &&
									(helix2 < coiled_coil_helices[c]))
									helix2++;
								if ((coiled_coil[c][helix2] == helix1) && (knobtype[k] > l))
									{
									sprintf(tmpstr2,"%c",alphabase[coiled_coil_helices[c]]+ helix2);
									sprintf(tmpstr,"%d",knobtype[k]);
									l = knobtype[k];
									}
								}

						if (i)
							{
							if (l)
								{
								strcat(result[2],tmpstr2);
								strcat(result[3],tmpstr);
								}
							else	{
								strcat(result[2],"!");
								sprintf(tmpstr,"%d",i);
								strcat(result[3],tmpstr);
								}
							}
						else
							{
							strcat(result[2],"-");
							strcat(result[3],"-");
							}

						}

				/* print out the 4 results strings */
				for (i = 0; i < 4; i++)
					printf("%s %s\n",result_name[i],result[i]);

				/* print out tad-signature */

				printf("repeats  %2d non-canonical interrupts in %3d residues: ",
					nonc_breaks, region_length);
				for (i = 0; i <= tad_index; i++)
					{
					if (i) printf(",");
					printf("%d",tadlength[i]);
					}
				printf("\n");

				total_nonc_breaks += nonc_breaks;
				if (nonc_breaks && (region_length > longest_nonc)) longest_nonc = region_length;

				} /* end of helix for loop */

			/* if a rasmol script file is being created, add the
			the definitions of sets of residues of each heptad position
			for this particular coiled coil*/
			if (par[par_r] != NULL)
				{
		
				/* this adds the definitions of sets of residues of each heptad position;
				   i represents the register a..g */
				for (i = 0; i < 7; i++)
					{
					/* j counts the number of residues which have this register; r is the residue id */
					j = 0; for (r = 0; r < residue_index; r++) if (tad_register[r][c] == i+'a') j++;
					if (j > RASMOL_WRAP) k = 0; 
					fprintf(rasmol_file,"define register_%d%c",c,i+'a');
					if (j > RASMOL_WRAP) fprintf(rasmol_file,"_%d ",k);
					else fprintf(rasmol_file," (");
					l = 0;
					for (r = 0; r < residue_index; r++) if (tad_register[r][c] == i+'a')
						{
						if (l == RASMOL_WRAP)
							{
							l = 0; k++; fprintf(rasmol_file,"\ndefine register_%d%c_%d ",c,i+'a',k);
							}
						if (l++) fprintf(rasmol_file,",");

						/* N.B. when writing the RasMol script file, insertion codes (iCode)
						are assumed to be null; not sure if RasMol deals with them correctly */

						fprintf(rasmol_file,"%d:%c",helix_residue_no[r],helix_chain[helix_no[r]]);
						}
					if (j > RASMOL_WRAP)
						{
						fprintf(rasmol_file,"\ndefine register_%d%c (",c,i+'a');
						for (k = 0; k <= (j -1) / RASMOL_WRAP; k++)
							{
							if (k) fprintf(rasmol_file,",");
							fprintf(rasmol_file,"register_%d%c_%d",c,i+'a',k);
							}
						}
					/* see the notes below re user-friendly RasMol set names */
					fprintf(rasmol_file,") and (sidechain,*.ca)\ndefine register_%s%c register_%d%c\ndefine reg_%s%c register_%d%c\n",
						roman[true_ccs],i+'a',c,i+'a',roman[true_ccs],i+'a',c,i+'a');
					}		
				}

			/* add some more user-friendly definitions of RasMol sets.
			So far, each coiled coil has been defined by reference to its serial number (c). However, not all the coiled
			coils are 'true' coiled coils, so they won't show up. For example there might be two coiled coils in a structure,
			with IDs 4 and 6. The user can't be expected to know these numbers. So duplicate definitions are added here,
			using roman numerals */

			if (par[par_r] != NULL) 
				{
				fprintf(rasmol_file,"define coiled_coil_%s coiled_coil%d\ndefine cc_%s coiled_coil%d\necho cc_%s consists of helices ",
					roman[true_ccs],c,roman[true_ccs],c,roman[true_ccs]);
				for (helix = 0; helix < coiled_coil_helices[c]; helix++)
					{
					if (helix) fprintf(rasmol_file,",");
					fprintf(rasmol_file,"%d",coiled_coil[c][helix]);
					}
				fprintf(rasmol_file,"\n");
				}
			} /* end of if ((coiled_coil_subset[c] == -1) && ((coiled_coil_helices[c] > 2) || (coiled_coil_tally[c] > 1))) */

		} /* end of c for loop REDUNDANT BRACKET*/

	/* if a rasmol script file is being created, add the
	the definitions of sets of residues of each heptad position
	for *all* coiled coils*/
	if (par[par_r] != NULL)
		for (i = 0; i < 7; i++)
			{
			fprintf(rasmol_file,"define register_%c",i+'a');
			for (c = 0; c < true_ccs; c++)
				{
				if (c) fprintf(rasmol_file,",");
				fprintf(rasmol_file," reg_%s%c",roman[c+1],i+'a');
				}
			fprintf(rasmol_file,"\ndefine reg_%c register_%c\n",i+'a',i+'a');
			}		


	/* report if there have been any non-canonical tad signatures */

	if (total_nonc_breaks)
		printf("%s c %5.2f e %d REPEATS: %2d NON-CANONICAL TAD-INTERRUPTS (LONGEST MOTIF %3d RESIDUES)\n",
			 code, cutoff2, extend, total_nonc_breaks, longest_nonc);

	return true_ccs;

	} /* end of function find_register */

int terminal_orientation(int helix1, int helix2, int residue_index)
	{
	int i,/*t_orientation,*/ r, N1, C1, N2, C2;
	/* N1..C2 are the *Serial numbers* of the N-terminal residue of
	helix 1, etc ; because the helix_start and helix_end values are
	the residue numbers as they appear in the PDB */
	float vector1[3],vector2[3],magnitude1,magnitude2,result;

	N1 = -1; N2 = -1; C1 = -1; C2 = -1;

	for (r = 0; r < residue_index; r++)
		{
		if ((helix_residue_no[r] == helix_start[helix1]) && 
			(helix_residue_iCode[r] == helix_start_iCode[helix1]) &&
			(helix_no[r] == helix1))
			N1 = r;
		if ((helix_residue_no[r] == helix_end[helix1]) &&
			(helix_residue_iCode[r] == helix_end_iCode[helix1]) &&
			(helix_no[r] == helix1))
			C1 = r;
		if ((helix_residue_no[r] == helix_start[helix2]) &&
			(helix_residue_iCode[r] == helix_start_iCode[helix2]) &&
			(helix_no[r] == helix2))
			N2 = r;
		if ((helix_residue_no[r] == helix_end[helix2]) &&
			(helix_residue_iCode[r] == helix_end_iCode[helix2]) &&
			(helix_no[r] == helix2))
			C2 = r;
		}

	if ((N1 == -1) || (C1 == -1) || (N2 == -1) || (C2 == -1))
		{
		printf("could not find first and last residues of helices %d and %d (PDB numbers are %d..%d and %d..%d)\n",
			helix1,helix2,helix_start[helix1],helix_end[helix1],helix_start[helix2],helix_end[helix2]);
		exit(1);
		}

/*printf("N1 = %d;\tC1 = %d;\tN2 = %d;\tC2 = %d\n",N1,C1,N2,C2);*/

	magnitude1 = 0.0; magnitude2 = 0.0; result = 0.0;

	for (i = 0; i < 3; i++)
		{
		vector1[i] = coord[refatom0[C1]][i] - coord[refatom0[N1]][i];
		magnitude1 += vector1[i]*vector1[i];
		vector2[i] = coord[refatom0[C2]][i] - coord[refatom0[N2]][i];
		magnitude2 += vector2[i]*vector2[i];
/*printf("vector1[%d] = %8.3f, vector2[%d] = %8.3f\n",i,vector1[i],i,vector2[i]); */
		}

	magnitude1 = sqrt(magnitude1);
	magnitude2 = sqrt(magnitude2);

/*printf("magnitude1 = %8.3f, magnitude2 = %8.3f\n",magnitude1,magnitude2);*/

	for (i = 0; i < 3; i++) result +=vector1[i]*vector2[i];

	result = 180.0 * acosf(result/(magnitude1 * magnitude2)) / M_PI;

	if (!setflag[flag_q]) printf("\tangle between helices %2d and %2d is %8.3f\t",
		helix1,helix2,result);

	if (result < 90.0)
		return 0;
	else	return 1;
	}

int orientation_of_helices(int helix1, int helix2,
	int orientation[], int orientation_first_helix[], int orientation_second_helix[], int pairs)
	{
	/* NB helix1 and helix2 are *helix serial numbers*, not the serial numbers with the context of
	the coiled coil (because of course these two helices might simultaneously be in two coiled
	coils; its not a coiled-coil dependent thing */
	int i;
	if (helix1 > helix2)
		{
		i = helix1;
		helix1 = helix2;
		helix2 = i;
		}

	i = 0;
	while ((i < pairs) && ((orientation_first_helix[i] != helix1) || (orientation_second_helix[i] != helix2)))
		i++;

	if ((orientation_first_helix[i] != helix1) || (orientation_second_helix[i] != helix2))
		{
		printf("!!!! orientation of helix %d with respect to helix %d has not been determined !!!!\n",helix1,helix2);
		exit(1);
		}
	return orientation[i];
	}

void set_orientation_of_helices(int helix1, int helix2, int new_orientation,
	int orientation[], int orientation_first_helix[], int orientation_second_helix[], int pairs)
	{
	/* NB helix1 and helix2 are *helix serial numbers*, not the serial numbers with the context of
	the coiled coil (because of course these two helices might simultaneously be in two coiled
	coils; its not a coiled-coil dependent thing */
	int i;
	if (helix1 > helix2)
		{
		i = helix1;
		helix1 = helix2;
		helix2 = i;
		}

	i = 0;
	while ((i < pairs) && ((orientation_first_helix[i] != helix1) || (orientation_second_helix[i] != helix2)))
		i++;

	if ((orientation_first_helix[i] != helix1) || (orientation_second_helix[i] != helix2))
		{
		printf("!!!! set_orientation_of_helices: orientation of helix %d with respect to helix %d has not been determined !!!!\n",helix1,helix2);
		exit(1);
		}

	orientation[i] = new_orientation;

	}

void check_extremes_of_hole(int knob, int c /* the coiled coil ID */)
	{
	int h,i;
	/* given the ID of a knob, each of the four residues of the hole into which
	it fits are checked; if any are the most extreme (N- or C-terminal)
	residues on that helix (the one with the hole, not the knob) -evaluated by
	comparing with the current values of the coiled_coil_begin[][][1] and
	coiled_coil_end[][][1] arrays - then one (or both) of these arrays is
	updated appropriately */

	/* first determine which helix the hole corresponding to the knob is on,
	and where this comes in the coiled_coil[c] list */

	i = helix_no[hole[knob][0]];
	h = 0;
	while ((h < coiled_coil_helices[c]) && (coiled_coil[c][h] != i)) h++;
	if (coiled_coil[c][h] != i) { printf("check_extremes_of_hole: oops !\n"); exit (1);}

	for (i = 0; i < 4; i++)
		{
		if ((coiled_coil_begin[c][h][1] == -1) ||
			(hole[knob][i] < coiled_coil_begin[c][h][1]))
			coiled_coil_begin[c][h][1] = hole[knob][i];
		if ((coiled_coil_end[c][h][1] == -1) ||
			(hole[knob][i] > coiled_coil_end[c][h][1]))
			coiled_coil_end[c][h][1] = hole[knob][i];
		}

	}

char relative_register(char reg, int offset)
	{
	char new;
	new = reg + offset;
/*printf("relative_register(%c,%d) = %c = ",reg,offset,new);*/
	while (new > 'g') new -= 7;
	while (new < 'a') new += 7;
/*printf("%c\n",new);*/
	return new;
	}
