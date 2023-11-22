#include <stdlib.h>
/*

					SOCKET
					 v3.02

					knobs.c

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

/*					knobs.c
					-------

*/

#include "socket.h"

void find_knobs_and_holes(int residue_index)
	{
	int i,j,k,sub_index;
	enum boolean is_first;
/* this array of colours written to the rasmol script, if specified,
used to be local to find_knobs_and_holes(), but is now required
by another short routine (in main()) */
char *rasmol_colour[RASMOL_COLOURS] = {"red","green","blue","magenta",
					"cyan","yellow","orange","white"};


	char sub;

	knob_index = 0;

	for (i = 0; i < helix_index-1; i++)
		for (j = i+1; j < helix_index; j++)
			if (i != j)
				{
				if (setflag[flag_v] || setflag[flag_l]) printf("\nhelices %3d,%3d:\n",i,j);
				reset_contacts();
				measure_helix_pair(i,j,residue_index);
				best_kih(residue_index);
				report_kih(residue_index);
				}
	if (setflag[flag_debug])
		{
		printf("helix %d checked; DUMPING KNOBS\n",i); dumpknobs(-1);
		printf("all helices checked; DUMPING KNOBS\n");
		}

	check_complementarity();

	check_duplication();

	/* tally number of knobs (and holes) of each type; total is cumulative*/

	for (i = 0; i < helix_index; i++) for (j = 0; j < 5; j++)
		{	n_knobtype[i][j] = 0;	n_holetype[i][j] = 0; }


	for (k = 0; k < knob_index; k++)
		{
		for (j = 0; j <= knobtype[k]; j++)
			{
			n_knobtype[helix_no[knob[k]]][j]++;
			n_holetype[helix_no[hole[k][0]]][j]++;
			}
		}

	if ((setflag[flag_l]) || (setflag[flag_v])) for (i = 0; i < helix_index; i++)
			for (j = 0; j < 7; j++) printf("\thelix %d has %d type %d knobs and %d type %d holes\n",
				i,n_knobtype[i][j],j,n_holetype[i][j],j);
				


	/* Print out the list of knobs in each helix; if a rasmol script is
		being written, write the appropriate lines to it which define
		the set of knobs for each helix ('knobs0', 'knobs1' etc)
		and displays them */

	if (!setflag[flag_q]) printf("\nThese are the knobs and holes:\n\n");
	for (j = 0; j < helix_index; j++)
		{
		if (!setflag[flag_q]) printf("knobs in helix %d:\n",j);
		if (par[par_r] && n_knobtype[j][knob_threshold]) 
			{
			sub_index = 0;
/*printf("helix %d has %d knobs\n",j,n_knobs[j]);*/
			if (n_knobtype[j][knob_threshold] < RASMOL_WRAP)
				{sub = '\0'; }
			else {sub = 'a'; }

			
			}

		for (i = 0; i < knob_index; i++)
			if ((helix_no[knob[i]] == j) && (knobtype[i] >= knob_threshold) )
				{

				if (!setflag[flag_q]) printf("%d) %d (%s %d:%c, iCode='%c', helix %d) type %d\t",
					i,knob[i],helix_residue_name[knob[i]],helix_residue_no[knob[i]],
					helix_chain[helix_no[knob[i]]], helix_residue_iCode[knob[i]],
					helix_no[knob[i]],knobtype[i]);

				if (par[par_r] && n_knobs[j]) {

					/* check that the rasmol line defining this set is not too long-
						if it is, start a new line */

					if (!(sub_index % RASMOL_WRAP))
						{
						if (sub_index) fprintf(rasmol_file,") and (sidechain,*.ca)\n");

						if (sub == '\0')
							{
							fprintf(rasmol_file,"define knobs%d (",j);
							}
						else
							{
							fprintf(rasmol_file,"define knobs%d%c (",j,sub);
							sub++;
							}
						is_first = true;
						}

					sub_index++;

					if (!is_first) fprintf(rasmol_file,",");

					/* N.B. when writing the RasMol script file, insertion codes (iCode)
					are assumed to be null; not sure if RasMol deals with them correctly */

					fprintf(rasmol_file,"%d:%c",helix_residue_no[knob[i]],helix_chain[helix_no[knob[i]]]);

					is_first = false;

					}

				if (!setflag[flag_q]) printf("(hole: ");

				for (k = 0; k < 4; k++)
					{
					if ((k) && (!setflag[flag_q])) printf(",");
					if (!setflag[flag_q]) 
						printf(" %s %d:%c iCode='%c'",helix_residue_name[hole[i][k]],helix_residue_no[hole[i][k]],
					helix_chain[helix_no[hole[i][k]]],helix_residue_iCode[hole[i][k]]);
					}

				if (!setflag[flag_q]) printf(" helix %d) packing angle %8.3f\n",helix_no[hole[i][0]],angle[i]);
				}

		--sub_index;

		if (par[par_r] && n_knobtype[j][knob_threshold]) 
			{
			fprintf(rasmol_file,") and (sidechain,*.ca)\n"/*,j,n_knobs[j]*/);
			if (sub != '\0')
				{
				fprintf(rasmol_file,"define knobs%d (",j);
				for (i = 0; i <= (sub_index / RASMOL_WRAP); i++)
					{
					if (i) fprintf(rasmol_file,",");
					fprintf(rasmol_file,"knobs%d%c",j,'a'+i);
					}
				fprintf(rasmol_file,")\n");
				}
			}
		}

	/* Print out the list of holes in each helix; if a rasmol script is
		being written, write the appropriate lines to it which define
		the set of holes for each helix ('holes0', 'holes1' etc) */

	for (j = 0; j < helix_index; j++)
		{

		if (!setflag[flag_q]) printf("holes in helix %d:\n",j);
		if (par[par_r] && n_holetype[j][knob_threshold])
			{
			sub_index = 0;
			if ((n_hole_res[j] * 4) <= RASMOL_WRAP)
				{sub = '\0'; }
			else {sub = 'a'; }
			
			}

		for (i = 0; i < knob_index; i++)
			if ((helix_no[hole[i][0]] == j) && (knobtype[i] >= knob_threshold))
				{
				for (k = 0; k < 4; k++)
					{
					if ((k) && (!setflag[flag_q])) printf(",");
					if (!setflag[flag_q]) 
						printf(" %s %d:%c iCode='%c'",helix_residue_name[hole[i][k]],helix_residue_no[hole[i][k]],
					helix_chain[helix_no[hole[i][k]]], helix_residue_iCode[hole[i][k]]);
					if (par[par_r] && n_hole_res[j]) {

						/* check that the rasmol line defining this set is not too long-
							if it is, start a new line */
						if (!(sub_index % RASMOL_WRAP))
							{
							if (sub_index) fprintf(rasmol_file,") and (sidechain,*.ca)\n");

							if (sub == '\0')
								{
								fprintf(rasmol_file,"define holes%d (",j);
								}
							else
								{
								fprintf(rasmol_file,"define holes%d%c (",j,sub);
								sub++;
								}
							is_first = true;
							}
						sub_index++;

						if (!is_first) fprintf(rasmol_file,",");
						fprintf(rasmol_file,"%d:%c",helix_residue_no[hole[i][k]],helix_chain[helix_no[hole[i][k]]]);
						is_first = false;
						}
					}
				if (!setflag[flag_q]) 
					printf(" (knob: %d (%s %d:%c, helix %d))\n",knob[i],helix_residue_name[knob[i]],helix_residue_no[knob[i]],
				helix_chain[helix_no[knob[i]]],helix_no[knob[i]]);
				}
		if (par[par_r] && n_holetype[j][knob_threshold])
			{
			fprintf(rasmol_file,") and (sidechain,*.ca)\n");

			if (sub != '\0')
				{
				fprintf(rasmol_file,"define holes%d (",j);

				for (i = 0; i <= (--sub_index / RASMOL_WRAP); i++)
					{
					if (i) fprintf(rasmol_file,",");
					fprintf(rasmol_file,"holes%d%c",j,'a'+i);
					}

				fprintf(rasmol_file,")\n");
				}
			}
		}

	if (!setflag[flag_q]) printf("\n");


	/* Write to the rasmol script file the commands to highlight the helices and knobs */

	if (par[par_r] != NULL)
		{

		fprintf(rasmol_file,"restrict ");

		for (i = 0; i < helix_index; i++) 
			{
			if (i) fprintf(rasmol_file,",");
			fprintf(rasmol_file,"helix%d",i);
			}

		fprintf(rasmol_file,"\nstrands 1\ncolour chain\n");

	/* as well as printing out the RasMol commands, 
		add user-friendly 'echo' statements which state the extent and colour of
		each helix, and which helices belong to which coiled coil */

		fprintf(rasmol_file,"echo\necho *** HELICES ***\necho\n");

		for (i = 0; i < helix_index; i++)
			{
			j = i%RASMOL_COLOURS;

			/* N.B. when writing the RasMol script file, insertion codes (iCode)
			are assumed to be null; not sure if RasMol deals with them correctly;
			so, the start and end of each helix are referred to here only in terms
			of residue serial number in combination with chain identifier */

			fprintf(rasmol_file,"select helix%d\ncolour %s\necho helix%d (%d-%d:%c) is %s\n",
				i,rasmol_colour[j],i,helix_start[i],helix_end[i],helix_chain[i],rasmol_colour[j]);
			}

		is_first = true;

		j = 0;

		for (i = 0; i < helix_index; i++) 
			if (n_knobtype[i][knob_threshold]) 
				{
				j = 1;
				if (is_first) fprintf(rasmol_file,"select ");
				else fprintf(rasmol_file,",");
				fprintf(rasmol_file,"knobs%d",i);
				is_first = false;
				}
		if (j) fprintf(rasmol_file,"\nwireframe 100\n");
		}

	}


void reset_contacts()
	{ int i; for (i = 0; i < MAX_RESIDUES; i++) n_contacts[i] = 0;}

void best_kih(int residue_index)
/* If there are > 4 residues in contact with a side chain, then find the
	best (most 'hole-like') group of four;
	the first group of 4 putative hole residues with a 3,1,3 spacing
	(ie X--XX--X) will be used as the hole; if no such group exists then
	any group where the second and third residues are consecutive will
	be used.
	In either case, the contact[j][] array will be rearranged so that
	the 4 residues constituting the best group will occupy the first
	4 elements. This is because later functions which test the hole
	spacing pattern assume that there is an x,1,y spacing pattern (anything
	but the first 4 elements are ignored);
	groups which do not match this will be discarded*/
	{
	int pos[4], best[4], i,j, ok;

	for (j = 0; j < residue_index; j++)
		{
		if (n_contacts[j] > 4)
			{
			printf("residue %d (%s %d:%c) has > 4 contacts; looking for holes\n",
				j,helix_residue_name[j],helix_residue_no[j],helix_chain[helix_no[j]]);
			pos[0] = 0; best[0] = -1;
			i = 0;
			while (++i < 4) { best[i] = -1; pos[i] = pos[i-1] + 1;	}
			ok = 1;

			while (ok)
				{
				/*for (i = 0; i < 4; i++) printf("\t%d",pos[i]);
				printf("\n");*/

				/* There are N (= n_contacts[j]) positions to test; all possible
				groups of 4 of them are tested (until a 3,1,3 spacing is found)
				The current selection is held in pos[]. */

				if (contact[j][pos[2]] == contact[j][pos[1]] + 1)
					{
					/* save this combination */
					for (i = 0; i < 4; i++) best[i] = pos[i];
					/* now test it for 3,1,3 spacing */
					if ((contact[j][pos[1]] == contact[j][pos[0]] + 3) &&
						(contact[j][pos[3]] == contact[j][pos[2]] + 3))
						{
						ok = 0; break;
						}
					}

				i = 3;
				while (++pos[i] == n_contacts[j] + i - 3)
					{
					if (--i < 0) { ok = 0; break;}
					}

				if (i < 3)
					while (++i < 4) pos[i] = pos[i-1] + 1;
				}

			if (best[0] != -1)
				{

				printf("contacts: ");
				for (i = 0; i < n_contacts[j]; i++)
					printf("%d) %s %d:%c iCode='%c'\t",i, helix_residue_name[contact[j][i]],helix_residue_no[contact[j][i]],
						helix_chain[helix_no[contact[j][i]]], helix_residue_iCode[contact[j][i]]);

					printf("\n- includes hole:\n");
					for (i = 0; i < 4; i++)
						{
						contact[j][i] = contact[j][best[i]];
						printf("%s %d:%c iCode='%c'\t",helix_residue_name[contact[j][i]],helix_residue_no[contact[j][i]],
						helix_chain[helix_no[contact[j][i]]], helix_residue_iCode[contact[j][i]]);
						}

					printf("\n\n");
				}

			}
		}
	}

void report_kih(int residue_index)
	{
	int i,j;

	if (setflag[flag_v] || setflag[flag_l])
		{
		for (i = 0; i < residue_index; i++)
			if (n_contacts[i] > 3)
				{
				printf("knobs\t\t\t\tholes\n-----\t\t\t\t-----\n\n");
				break;
				}
		}

	for (i = 0; i < residue_index; i++)

		if (n_contacts[i] > 3)
			{
			/* this residue (i) is touching at least 4 other side chains
			(necessarily on the same helix); so residue (i) is a knob */

			n_knobs[helix_no[i]]++;
			if (setflag[flag_v] || setflag[flag_l])
				printf("%d (%s %d:%c iCode='%c', helix %d)\t",i,helix_residue_name[i],
			helix_residue_no[i],helix_chain[helix_no[i]],
			helix_residue_iCode[i], helix_no[i]);

			if (knob_index == MAX_KNOBS)
				{
				printf("Maximum number of knobs (%d) exceeded\n",MAX_KNOBS);
				exit(1);
				}

			for (j = 0; j < n_contacts[i]; j++)
				{
				if (setflag[flag_v] || setflag[flag_l])
					printf("\t%d (%s %d:%c iCode='%c', helix %d) ",
					contact[i][j],
					helix_residue_name[contact[i][j]],
					helix_residue_no[contact[i][j]],
					helix_chain[helix_no[contact[i][j]]],
					helix_residue_iCode[contact[i][j]],
					helix_no[contact[i][j]]);

				hole[knob_index][j] = contact[i][j];
				hole_distance[knob_index][j] = contact_distance[i][j];
				}

			knobtype[knob_index] = 0;
			knob[knob_index] = i;
			n_hole_res[helix_no[contact[i][0]]]++;

			if (setflag[flag_debug])
				printf("n_hole_res[%d] is now %d\n\n",helix_no[contact[i][0]],n_hole_res[helix_no[contact[i][0]]]);

/* do the determining of 'knobtype' here... */

			/* numbering of the 4 hole residues should follow this pattern: x,x+n,x+n+1,x+n+1+m
			where n and m are either 3 or 4.
			in canonical heptad repeats, 'a' and 'd' knobs have hole residues: x,x+3,x+4,x+7 (n = m = 3)
			- numbering can be checked in terms of internal numbers, or actual resSeq number;
			latter is not always reliable, because of discontinuities in the sequence of
			residue serial numbers, or inserted residues (which will have non-null iCodes)*/

			/* any holes where the 2nd and 3rd residues are not consecutive are disregarded: */

			if (contact[i][2] == (contact[i][1] + 1))
				{
				knobtype[knob_index]++;

				/* give a warning if its not a x,x+3,x+4,x+7 spacing */

				if ((contact[i][1] != contact[i][0] + 3) || (contact[i][3] != contact[i][2] + 3))
					{printf("!!!!odd-knob\n");	}

				/* assuming that both the XXX cutoff and insertion-cutoff are 7.0A, then:
					the 4 hole side chains must all be < 7.0A from the knob side chain
					  (as defined by centre of mass);
					if the knob is in, and not across, the hole, then the end of the knob
					  side chain would usually be expected to be even closer to, but certainly
					no more than 7.0A from, the 4 centres-of-mass of the hole side chains;
					it is possible for the end of the k side chain to be > 7.0A from one or
					  two of the h side-chain-centres; but if the k is in the h, then this
					  increase in distance must be matched by a very nearly identical decrease
					  in distance from one or two of the other h side chain centres;
					therefore, the *mean* distance - between the k end and each of the 4 h centres-
					  is necessarily < 7.0A, if the k is in the h;
					this mean distance is calculated by the measure_k_end_h_CA() function. */

				if (measure_k_end_h_CA(knob_index) < 7.0) /* XXX get rid of this naughty hardcoded constant */
					{
					knobtype[knob_index]++;
					}

				else if (setflag[flag_v] || setflag[flag_l]) printf("knob not in hole");
				angle[knob_index] = packing_angle(i,contact[i][1],contact[i][2]);
				}
			if (setflag[flag_v] || setflag[flag_l]) printf("\n");
			if (knobtype[knob_index]) knob_index++;
			}
	}

void check_complementarity()
	{
	int i,j,k,l,m,n,comp;
	for (i = 0; i < knob_index; i++) n_compknob[i] = 0;

	daisy_chains = 0;
	for (i = 0; i < MAX_DAISY_CHAINS; i++)
		for (j = 0; j < MAX_DAISIES; j++) daisy_chain[i][j] = -1;

	for (i = 0; i < knob_index; i++)
		{
		if (setflag[flag_v] || setflag[flag_l])
			printf("checking knob\t%d (%s\t%d:%c iCode='%c', helix\t%d)",i,helix_residue_name[knob[i]],
			helix_residue_no[knob[i]], helix_chain[helix_no[knob[i]]],
			helix_residue_iCode[knob[i]], helix_no[knob[i]]);

		knob_order[i] = -1;

		for (j = 1; j < 3; j++)
			{
			comp = complementary(i,j);

			if (comp != -1)
				{
				if (n_compknob[i] == 3) { printf("!!this knob has > 3 **complementary** knobs - something is seriously wrong (even 3 is pretty unbelievable) - have you used a stupidly large cutoff? These are the complementary knobs:\n");
					for (j = 0; j < 3; j++) printf("%d, ",compknob[i][j]); printf("%d\n",comp);
					exit(1);
					}
				compknob[i][n_compknob[i]++] = comp;
				knobtype[i] += 2;
				if ((knobtype[i] > 4) && (setflag[flag_v] || setflag[flag_l])) printf(" DOUBLE type %d",knobtype[i]);

				knob_order[i] = 2;

				if (setflag[flag_v] || setflag[flag_l]) printf(" complementary with knob %d",comp);
				}
			}
		if (setflag[flag_v] || setflag[flag_l]) printf("\n");

		for (n = 1; n > -2; n -= 2)
			{

			if (daisy_chains == MAX_DAISY_CHAINS)
				{
				printf("Maximum number of daisy chains (%d) exceeded\n",MAX_DAISY_CHAINS);
				exit(1);
				}

			k = check_daisy_chain(i,0,daisy_chain[daisy_chains],n);

			if (k >= 0)
				{
				knob_order[i] = k;
			/* its possible that the daisy chain thats just been found is an identical copy of a previously
			discovered chain, but with the knobs listed in a different order. The following checks this and
			if its the case, the latest chain found is discarded. NB this is distinct from the case of a
			knob being part of two different chains, which is possible eg in gp41-like arrangements where
			the central 3 helices each simultaneously belong to two different 3-strand coiled coils.

			The first thing to do is sort the knobs belonging to the chain, which makes it easier to compare them */

				/* start of sort routine */

				l = 1;
				while (l)
					{
					l = 0;
					for (j = 1; j < k; j++)
						if (daisy_chain[daisy_chains][j] < daisy_chain[daisy_chains][j-1])
							{
							l = daisy_chain[daisy_chains][j];
							daisy_chain[daisy_chains][j] = daisy_chain[daisy_chains][j-1];
							daisy_chain[daisy_chains][j-1] = l;
							l = 1;
							}
					}
				/* end of sort routine */


				k = 0; /* k flags whether or not an identical daisy-chain has been found */
				j = 0; while ((j < daisy_chains) && (!k))
					{
					k = 1; l = 0;
					while ((daisy_chain[daisy_chains][l] != -1) && (k) && (l < MAX_DAISIES))
						{
						if (daisy_chain[daisy_chains][l] != daisy_chain[j][l]) k = 0;
						l++;
						}
					j++;
					}
				if (k)	{for (j = 0; j < MAX_DAISIES; j++) daisy_chain[daisy_chains][j] = -1; }
				else daisy_chains++;
				}

			else for (j = 0; j < MAX_DAISIES; j++)
				daisy_chain[daisy_chains][j] = -1;

			}

		}



	/* only report 'daisy-chains' (closed, non pairwise complementary, cycles of knob-into-hole interactions),
		when all the daisy-chains have been collected (done in the above i=0..knob_index loop) */


	if (setflag[flag_v])
		for (i = 0; i < daisy_chains; i++)
			{
			for (j = 0; j < MAX_DAISIES; j++)
				{
				printf("daisy_chain[%2d][%2d] = %2d\n",i,j,daisy_chain[i][j]);
				}
			}

	/* first list by daisy chain */
	for (i = 0; i < daisy_chains; i++)
		{
		printf("daisy chain %2d : knobs ",i);
		j = -1; while ((daisy_chain[i][++j] != -1) && (j < MAX_DAISIES))
			printf("%2d (helix %2d)\t",daisy_chain[i][j],helix_no[knob[daisy_chain[i][j]]]);
		printf("\n");
		}

	/* first list by knob */
	for (i = 0; i < knob_index; i++)
		{
		l = -1; 
		for (j = 0; j < daisy_chains; j++)
			{
/*printf("searching daisy-chain %d for knob %d\n",j,i);*/
			k = 0; while ((daisy_chain[j][k] != -1) && (k < MAX_DAISIES)) 
				{
				if (daisy_chain[j][k] == i)
					{
/*printf("daisy_chain[%d][%d] = %d\n",j,k,daisy_chain[j][k]);*/
					if (l == -1)
						{
						if (!setflag[flag_q])
						printf("knob %3d (residue %d = %s %d:%c iCode='%c')",i,knob[i],helix_residue_name[knob[i]],
						helix_residue_no[knob[i]],helix_chain[helix_no[knob[i]]],helix_residue_iCode[knob[i]]);
						if (knobtype[i] < 3) {
							/* turn any knobs of type 1 or 2 into proper knobs (3 or 4) */
							knobtype[i] += 2;
							/* this while loop puts the knobs (apart from knob i itself)
							into the compknob array for knob i (while theres room) */
							}
						m = 0;
						while ((daisy_chain[j][m] != -1) && (m < MAX_DAISIES))
							{
							if ((daisy_chain[j][m] != i) && (n_compknob[i] < MAX_COMPKNOBS))
								compknob[i][n_compknob[i]++] = daisy_chain[j][m];
							m++;
							}
						}
					l = j;
					}
				k++;
				}
			if (l == j) {
				if (!setflag[flag_q]) {printf(" forms a %d-knob cycle with knobs ",k);
				k = 0; while ((daisy_chain[j][k] != -1) && (k < MAX_DAISIES)) {printf("%4d",daisy_chain[j][k]); k++;}
				printf("; ");}
				}
			}
		if ((l != -1) && (!setflag[flag_q])) printf("\n");

		/* it is possible that a knob in a daisy chain has ended up with an order of -2.
		For example, Leu 305:A in Stat3B 1bg1. In this structure there are knobs which
		are simultaneously part of 2 daisy chains:

		Leu 158:A is knob 1
		Gln 280:A is knob 15
		Leu 221:A is knob 17
		Leu 305:A is knob 24
		Ile 281:A is knob 28

		These form two 3-membered daisy chains
		 i) 17 -> 15 -> 1
		ii) 17 -> 28 -> 24

		Note that knobs 15 and 28 are consecutive residues; they form
		the two sides of the hole into which knob 17 fits.

		When daisy chains are being traced by the check_daisy_chains()
		function, then the next knob found after knob 17 will always
		be knob 15, because it is sequentially first; daisy chain (i)
		will be traced, terminating back at knob 17 (this happens
		for example when knob 17 is checked).

		The only case where this won't happen is when checking knob 28.
		The trace will find knob 24, then 17, and will then terminate
		because one of the two sides of the hole into which knob 17
		fits is the first knob in the chain (28).

		The problem arises when knob 24 is checked. The trace will find
		knob 17, then knob 15, and not knob 28. So the trace enters
		daisy chain (i):
		24 -> 17 -> 15 -> 1

		- this loop closes back on knob 17 of course, which is the second
		knob in the chain, not the first. That is, this 4-membered loop
		does not constitute a closed chain itself, so the function
		check_daisy_chains() returns a null value (-2) when checking
		knob 24. So this knob won't be given its correct order, which
		is 3.

		The way round this is to go through the list of knobs, and if
		any still have a null order (-1) then to look for it in the list
		of daisy chains. If it appears in any daisy chain, then it
		should have its order set to the size of that daisy chain.

		It is just conceivable, but very unlikely, that a knob could be
		part of 2 different daisy chains of different order

		        / 6
		       /  ^
		      L   |
		1 -> 2 -> 5
		^    | 	 		 i) 1 -> 2 -> 4 -> 3
		|    v  
		3 <- 4  		ii) 2 -> 5 -> 6


		Suppose that again, 4 and 5 are the two sides of the hole into
		which knob 2 fits, and that 5 is sequentially first. Only when
		checking knob 4 will daisy chain (i) be traced. When checking
		knob 2, (ii) will be traced, giving knob 2 an order of 3.
		Arguably its order should be the highest of any daisy chain of
		which it is part, ie 4. However, the SOCKET program does not
		really legislate for knobs simultaneously being part of daisy
		chains of different order, and the find_register() function
		might not like it. All cases of knobs belonging to more than
		one daisy chain so far seen involve only 3-stranded arrangements.

		Therefore, the above will be assumed NOT to happen, which means
		that a check which finds the *highest-order* daisy chain
		involving a knob, is unnecessary. Instead, the search is for
		*any* daisy chain of which an as yet un-ordered knob is part.


		NB the above problem can also affect knobs which already have an order
		of > 0. Eg in 1aik.mmol, there is this situation:
		(knob numbers) 7 -> 11 -> 45 -> 21 -> 11....
		knob 7 = 642:C, 11 = 562:N, 45 = 562:A, 11 = 562:D
		However * knob 7 has already been given an order of 2 *, so this
		doesn't get overwritten with -2. Therefore, *all* knobs should be
		checked; their knob-order will only be changed if the order of
		the daisy chain they are in is greater than their current order
		 */

			for (j = 0; j < daisy_chains; j++)
				{
				k = 0;
				while ((daisy_chain[j][k] != -1)
					&& (daisy_chain[j][k] != i)
					&& (k < MAX_DAISIES))
					k++;
				if (daisy_chain[j][k] == i)
					/* knob i is in daisy chain j; what order is daisy chain j? */
					{
					k = 0; while ((daisy_chain[j][k] != -1) && (k < MAX_DAISIES)) k++;
					if (k > knob_order[i])
						{
						if (!setflag[flag_q])
							printf("order of knob %d was %d; changing to %d\n",i,knob_order[i],k);
						knob_order[i] = k;
						}
					}
				}
		}

	}


int complementary(int knobno, int holeresno)
	{
	int i,done;
	done = -2;
	i = 0;
	/*printf("doing complementary(%d,%d)\n",knobno,holeresno);*/
	while ((done == -2) && (i < knob_index))
		{
/*printf("bigloop i = %d\n",i);*/
		while ((i < knob_index) && (knob[i] != hole[knobno][holeresno])) i++; 
		/* if the below is true, then the holeresno'th residue of the hole
			which contains knob knobno isnt a knob at all */
		if /*(knob[i] != hole[knobno][holeresno])*/ (i == knob_index) done = -1;
		else
			{
			/* so the hole residue is also a knob, but is it a knob with
			knob knobno as either the 2nd or 3rd residue of the hole? */
			if ((hole[i][1] == knob[knobno]) || (hole[i][2] == knob[knobno]))
				{
				done = i;
				}
			else i++;
			}
		}
	if (done == -2) done++;
	return done;
	}

