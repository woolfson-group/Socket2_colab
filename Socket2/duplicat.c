#include <stdlib.h>
/*

					SOCKET
					 v3.02

					duplicat.c

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

/*					duplicat.c
					----------
*/


#include "socket.h"
#include <string.h>
void check_duplication()
	{
	int i,j,copies,checked[MAX_KNOBS];

	/* the checked[] array is so that each knob appears only once in the
	duplicates table at most*/
	for (i = 0; i < knob_index; i++) checked[i] = 0;

	/* initialize duplicates table (global) */
	for (i = 0; i < MAX_DUPLICATES; i++)
		for (j = 0; j < 3; j++) duplicate_knobs[i][j] = -1;

	/* global */
	n_duplicate_knobs = 0;


	for (i = 0; i < knob_index - 1; i++)
		{
		/* for each set of duplicates, this is how many there are
		(should be no more than 2 */
		copies = 0;

		for (j = i+1; j < knob_index; j++)
			if ((!checked[j]) && (knob[i] == knob[j]))
				{
				if (!copies) {
					if (copies == 3)
						{
						printf("too many knob duplicates: knobs %d, %d, %d, %d are all residue %d (%s %d:%c, iCode='%c') ; only 3 knob duplicates can be stored; is your specified cutoff too high?\n\n",
duplicate_knobs[n_duplicate_knobs][0],duplicate_knobs[n_duplicate_knobs][1],duplicate_knobs[n_duplicate_knobs][2],j,knob[i],helix_residue_name[knob[i]],helix_residue_no[knob[i]],helix_chain[helix_no[knob[i]]],helix_residue_iCode[knob[i]]);
					exit(1);
					}
					printf("duplicate knobs: %3d",i);
					duplicate_knobs[n_duplicate_knobs][copies++] = i;
					}
				printf(",%3d",j);
				duplicate_knobs[n_duplicate_knobs][copies++] = j;
				checked[j] = 1;
				}
		if (copies) {printf("\n"); n_duplicate_knobs++;}
		}

	if (n_duplicate_knobs > 0) printf("sets of duplicate knobs:\n");
	for (i = 0; i < n_duplicate_knobs; i++)
		{
		printf("%2d)",i);
		j = 0; while ((j < 3) && (duplicate_knobs[i][j] != -1)) printf(" %3d",duplicate_knobs[i][j++]);
		j = duplicate_knobs[i][0];
		printf("\tare all residue %d (%s %d:%c, iCode='%c')\n",knob[j],helix_residue_name[knob[j]],helix_residue_no[knob[j]],helix_chain[helix_no[knob[j]]],helix_residue_iCode[knob[i]]);
		}
	}
