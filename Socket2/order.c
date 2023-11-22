/*

					SOCKET
					 v3.02

					order.c

					18-06-01
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

/*					order.c
					-------
*/


#include "socket.h"
#include <string.h>
#include <stdlib.h>
int determine_order(/*int residue_index*/)
	{
	/* determines the oligomerization state of the coiled coils */
	int i, j, k, l, m, helix[MAX_DAISIES];

	/* function uses these global variables:
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
	coiled_coils - the number of different coiled coils
	helix_order[MAX_HELICES] - the number of helices in the coiled coil

	and:		
	knob_index
	knob_threshold
	knob[]
	knobtype[]
	knob_order[]
	helix_residue_name[]
	helix_residue_no[]
	helix_chain[]
	helix_no[]
	daisy_chains
	daisy_chain[][]
	daisy_chain_cc[]
	*/

	coiled_coils = 0;

	/* first compile list of all pairwise coiled-coil interactions */

	for (i = 0; i < knob_index; i++)
		{
		if (!setflag[flag_q]) printf("knob %3d (residue %d = %s %d:%c iCode='%c') type %d order %d\n",
			i,knob[i],helix_residue_name[knob[i]],
			helix_residue_no[knob[i]],helix_chain[helix_no[knob[i]]],
			helix_residue_iCode[knob[i]],knobtype[i],knob_order[i]);

		if (knobtype[i] >= knob_threshold)
			{
			if (helix_no[knob[i]] > helix_no[knob[compknob[i][0]]])
				{
				helix[0] = helix_no[knob[compknob[i][0]]];
				helix[1] = helix_no[knob[i]];
				}
			else	{
				helix[0] = helix_no[knob[i]];
				helix[1] = helix_no[knob[compknob[i][0]]];
				}

			k = 0;
			/* l signals whether this combination is found in the list */
			l = 0;
			while ((k < coiled_coils) && (!l))
				{
				l = 1;
				for (m = 0; m < 2; m++) if (coiled_coil[k][m] != helix[m]) l = 0;
				if (l) coiled_coil_tally[k]++;
				k++;
				}
			if (!l)
				{
				if (coiled_coils == MAX_COILED_COILS)
					{ printf("maximum number of coiled coils (%d) exceeded\n\n",MAX_COILED_COILS); exit(1);}
				for (m = 0; m < 2; m++) coiled_coil[coiled_coils][m] = helix[m];
				coiled_coil_subset[coiled_coils] = -1;
				coiled_coil_helices[coiled_coils] = 2;
				coiled_coil_tally[coiled_coils] = 1;
				coiled_coils++;
				}
			}
			
		}
	if (!setflag[flag_q]) printf("\n");


	for (i = 0; i < coiled_coils; i++) coiled_coil_tally[i] /= 2;

	/* now compile list of all cyclical coiled-coil interactions (daisy chains)*/

	for (i = 0; i < daisy_chains; i++)
		{
		/* loop through all the daisy chains, compiling lists of each
		combination of helices - each combination is recorded as
		one entry in the coiled_coil[MAX_COILED_COILS] array; the
		number of different times this combination is found is recorded
		in coiled_coil_tally */

		j = 0;
		while ((daisy_chain[i][j] != -1) && (j < MAX_DAISIES))
			{
			/* the helix id of this knob in the daisy chain is
			helix_no[knob[daisy_chain[i][j]]] ; add it to
			the helix[] array (in the correct order) */
			k = 0;
			while ((k < j) && (helix[k] < helix_no[knob[daisy_chain[i][j]]])) k++;
			for (l = j - 1; l >= k; l--) helix[l+1] = helix[l];
			helix[k] = helix_no[knob[daisy_chain[i][j]]];
			j++;
			}
		if (setflag[flag_v]) { printf("daisy chain %d; helices",i);
		for (k = 0; k < j; k++) {printf("\t%d",helix[k]);}
		printf("\n");	}

		/* the list of j helices which are part of this daisy chain is now
		held in helix[] ; if this list is different to all those found so
		far, add it to the coiled_coil array */

		k = 0;
		/* l signals whether this combination is found in the list */
		l = 0;
		while ((k < coiled_coils) && (!l))
			{
			l = 1;
			/* printf("cf coiled coil %d:\t",k); */
			if (setflag[flag_v]) for (m = 0; m < j; m++) printf("%d v %d; ",coiled_coil[k][m],helix[m]);
			for (m = 0; m < j; m++) if (coiled_coil[k][m] != helix[m]) l = 0;
			if (l)
				{
				coiled_coil_tally[k]++;
				daisy_chain_cc[i] = k;
				}
			if (setflag[flag_v]) { if (l) printf("MATCHED\n"); else printf("different\n"); }
			k++;
			}
		if (!l)
			{
			if (coiled_coils == MAX_COILED_COILS)
				{ printf("maximum number of coiled coils (%d) exceeded\n\n",MAX_COILED_COILS); exit(1);}
			for (m = 0; m < j; m++) coiled_coil[coiled_coils][m] = helix[m];
			coiled_coil_subset[coiled_coils] = -1;
			coiled_coil_helices[coiled_coils] = j;
			coiled_coil_tally[coiled_coils] = 1;
			daisy_chain_cc[i] = coiled_coils;
			coiled_coils++;
			}
		}

	for (i = 0; i < coiled_coils; i++)
		{

		/* check each combination of helices (coiled coil) - it might be
		a subset of other coiled coils */
		for (k = 0; k < coiled_coils; k++)
			{
			if (i != k)
				{
				/* l signals that the helices in coiled_coil_helices[i]
				are all also present in array coiled_coil_helices[k] */
				l = 1;
				j = 0;
				while ((j < coiled_coil_helices[i]) && (l))
					{
					l = 0;
					m = 0;
					while ((m < coiled_coil_helices[k]) && (!l))
						if (coiled_coil[k][m++] == coiled_coil[i][j]) l = 1;
					j++;
					}
				if (l)	{
					/* coiled coil i is a subset of coiled coil k */
					if (coiled_coil_subset[k] == -1)
						coiled_coil_subset[i] = k;
					else coiled_coil_subset[i] = coiled_coil_subset[k];
					/* its just possible that other coiled coils have previously
					   been assigned as subsets of coiled coil i ; so they
					   must be found and their subsets changed to the new
					   coiled_coil_subset[i] 
					   eg coiled coil 1 consists of helices A, B
					      coiled coil 2 consists of helices A, B, C, D
					      coiled coil 3 consists of helices A, B, C
					      1 is a subset of 2 and 3; 3 is a subset of 2
					      (this arrangement is unlikely but does occur, eg 1sfc)
					   when assigning coiled coil 1, it will first
						have coiled_coil_subset[1] set to 2;
                                                later it will be set to 3.
                                           coiled coil 2 will remain with coiled_coil_subset[2] set to -1;
                                           coiled coil 3 will have coiled_coil_subset[3] set to 2.
                                           Therefore a second pass is needed to change coiled_coil_subset[1]
					   from 3 to coiled_coil_subset[3] = 2.
					   The point is, all values of coiled_coil_subset[] should be
					   either -1, or the ID of a coiled coil X whose coiled_coil_subset[X] = -1.
					   l is used to index this second pass. */
					for (l = 0; l < coiled_coils; l++)
						if (coiled_coil_subset[l] == i) coiled_coil_subset[l] = coiled_coil_subset[i];
					/* also, the coiled coil ID of any daisy chains (daisy_chain_cc) whose daisy_chain_cc
					is currently set to i, should also be changed to coiled_coil_subset[i];
					l is sued to index the daisy chains */
					for (l = 0; l < daisy_chains; l++)
						if (daisy_chain_cc[l] == i) daisy_chain_cc[l] = coiled_coil_subset[i];
					}
				}
			}

		printf("coiled coil %2d: %2d helices ",i,coiled_coil_helices[i]);
		for (j = 0; j < coiled_coil_helices[i]; j++)
			printf("%3d",coiled_coil[i][j]);
		printf("\tfrequency %d",coiled_coil_tally[i]);
/* !!!! CHOULD USE A CONSTANT INSTEAD OF '2' IN THE LINE BELOW*/
		if ((coiled_coil_tally[i] < 2) && (coiled_coil_helices[i] == 2))
			printf(" IGNORING");
		if (coiled_coil_subset[i] != -1)
			{
			printf(" (subset of coiled coil %2d)",coiled_coil_subset[i]);
			if (coiled_coil_helices[i] > 2)
				printf(" WARNING: THESE HELICES HAVE BOTH %d-STRANDED AND %d-STRANDED CHARACTERISTICS - THIS MAY MAKE THE REGISTER ASSIGNMENT UNRELIABLE",
					coiled_coil_helices[i],coiled_coil_helices[coiled_coil_subset[i]]);
			}
		printf("\n");
		}


	printf("\n");
/* assign an order (number of strands) to each helix = the oligomerization
	state of the coiled coil to which it belongs */
	for (i = 0; i < helix_index; i++)
		{
		helix_order[i] = 0;

		for (j = 0; j < coiled_coils; j++)
			if (coiled_coil_subset[j] == -1)
				{
				k = 0; while ((k < coiled_coil_helices[j]) && (coiled_coil[j][k] != i)) k++;

				if ((k < coiled_coil_helices[j]) && (coiled_coil_helices[j] > helix_order[i]))
					helix_order[i] = coiled_coil_helices[j];
				}

		if (helix_order[i] != 0) printf("helix %d is in a %d-stranded coiled coil\n",i,helix_order[i]);
		}
	if (!setflag[flag_q]) printf("\n");


/* report the coiled coils to which any daisy chains belong */
	for (i = 0; i < daisy_chains; i++)
		printf("daisy chain %2d is in coiled coil %2d\n",i,daisy_chain_cc[i]);

	return coiled_coils;
	}

void define_ras_coils()
	{
	int i,/*j,*/m;

	/* if a rasmol file is being created, define the set of coiled coils (the
		union of all coiled coils which arent subsets of other coiled coils */
	if (par[par_r] != NULL)
		{
		fprintf(rasmol_file,"define coiled_coils "); m = 0;
		for (i = 0; i < coiled_coils; i++) if ((coiled_coil_subset[i] == -1)
			&& ((coiled_coil_helices[i] > 2) || (coiled_coil_tally[i] >= 2)))
			{
			if (m++) fprintf(rasmol_file,",");
			fprintf(rasmol_file,"coiled_coil%d",i);
			}
		fprintf(rasmol_file,"\n");
		}
	}
