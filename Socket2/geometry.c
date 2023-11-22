/*

					SOCKET
					 v3.03

					geometry.c

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

/*					geometry.c
					----------

*/

#include "socket.h"

float packing_angle(int knobres, int holeres1, int holeres2)
/* returns the angle between these two vectors:
	 i) C-alpha to C-beta of knobres
	ii) C-alpha of holeres1 to C-alpha of holeres2
 - knobres, holeres1 and holres2 are the residue keys */
	{
	float A[3],B[3],magnitudeA,magnitudeB,result;
	int i;
	/* A is knobres C-alpha...C-beta vector */
	magnitudeA = 0; magnitudeB = 0; result = 0;
	for (i = 0; i < 3; i++)
		{
		A[i] = coord[refatom3[knobres]][i] - coord[refatom0[knobres]][i];
		magnitudeA += A[i]*A[i];
		B[i] = coord[refatom0[holeres1]][i] - coord[refatom0[holeres2]][i];
		magnitudeB += B[i]*B[i];
/*printf("A[%d] = %8.3f, B[%d] = %8.3f\n",i,A[i],i,B[i]);*/
		}
	magnitudeA = sqrt(magnitudeA);
	magnitudeB = sqrt(magnitudeB);
/*printf("magnitudeA = %8.3f, magnitudeB = %8.3f\n",magnitudeA,magnitudeB);*/
	for (i = 0; i < 3; i++) result +=A[i]*B[i];
	return 180.0 * acosf(result/(magnitudeA * magnitudeB)) / M_PI;
	}

void determine_centre_of_mass(int atom_index)
	{
	int i, j, last_residue_index, n_atoms;
	last_residue_index = -1;
	n_atoms = 0;
	for (i = 0; i < atom_index; i++)
		if (atom_res[i] != last_residue_index)
			{
			if (strcmp(atom_name[i]," N  ") && strcmp(atom_name[i]," C  ") && strcmp(atom_name[i]," O  ") && 
((setflag[flag_i]) || strcmp(atom_name[i]," CA "))  && (setflag[flag_i] || (strcmp(atom_name[i]," HA ") && 
strcmp(atom_name[i],"1HA ") && strcmp(atom_name[i],"2HA ") ) ) ) 
				{
				if (last_residue_index != -1)
					{for (j = 0; j < 3; j++) refatom2[last_residue_index][j] /= n_atoms;
					if (setflag[flag_debug]) printf("centre of mass: %8.3f , %8.3f , %8.3f\n\n",
refatom2[last_residue_index][0],refatom2[last_residue_index][1],refatom2[last_residue_index][2]);}
				if (setflag[flag_debug]) printf("first atom of this residue (%d) is atom %d \"%s\"\n",atom_res[i],i,atom_name[i]);
				last_residue_index = atom_res[i];
				n_atoms = 1;
				for (j = 0; j < 3; j++) refatom2[last_residue_index][j] = coord[i][j];
				}
			}
		else if (strcmp(atom_name[i]," N  ") && strcmp(atom_name[i]," C  ") && strcmp(atom_name[i]," O  ") && 
((setflag[flag_i]) || strcmp(atom_name[i]," CA "))  && (setflag[flag_i] || (strcmp(atom_name[i]," HA ") && 
strcmp(atom_name[i],"1HA ") && strcmp(atom_name[i],"2HA ") ) ) ) 
			{
			if (setflag[flag_debug]) printf("\tadding atom %d \"%s\" to residue (%d)\n",i,atom_name[i],atom_res[i]);
			n_atoms++;
			for (j = 0; j < 3; j++) refatom2[last_residue_index][j] += coord[i][j];
			}
	for (j = 0; j < 3; j++) refatom2[last_residue_index][j] /= n_atoms;
	if (setflag[flag_debug]) printf("centre of mass: %8.3f , %8.3f , %8.3f\n",
refatom2[last_residue_index][0],refatom2[last_residue_index][1],refatom2[last_residue_index][2]);
	}

void determine_end(int residue_index, int atom_index)
	{
	int i, j, k, n_atoms;

	for (i = 0; i < residue_index; i++)
		{
		if (setflag[flag_debug]) printf("finding end of residue %d\n",i);

		n_atoms = 0;

		for (k = 0; k < 3; k++) refatom1B[i][k] = 0.0;

		for (j = 0; j < 2; j++)
			if (refatom1[j][i] != -1)
				{
				n_atoms++;
				for (k = 0; k < 3; k++) refatom1B[i][k] += coord[refatom1[j][i]][k];
				}

		if (!n_atoms)
				{
				printf("residue %d (%s %d:%c iCode='%c') has no 'end' atoms\n\n",
					i, helix_residue_name[i],helix_residue_no[i],
				helix_chain[helix_no[i]], helix_residue_iCode[i]); 

				for (k = 0; k < atom_index; k++)
					{
					if (setflag[flag_debug])
						printf("atom_res[%d] (of %d) is %d, vs residue %d\n",
						k,atom_index,atom_res[k],i);
					if (atom_res[k] == i) refatom1[0][i] = k;
					if (atom_res[k] > i) break;
					}

				for (k = 0; k < 3; k++) refatom1B[i][k] = coord[refatom1[0][i]][k];

				n_atoms++;

				printf("Using atom %d (%s %d%c iCode='%c' %s) as end atom for this residue\n", refatom1[0][i],
				helix_residue_name[i],helix_residue_no[i],helix_chain[helix_no[i]],
				helix_residue_iCode[i],atom_name[refatom1[0][i]]);

				}

		if (n_atoms > 1) for (k = 0; k < 3; k++) refatom1B[i][k] /= n_atoms;

		if (setflag[flag_debug]) printf("residue %d (%s %d:%c iCode='%c') end (pseudo)atom: %8.3f , %8.3f , %8.3f\n", 
		i, helix_residue_name[i],helix_residue_no[i], helix_chain[helix_no[i]],
		helix_residue_iCode[i], refatom1B[i][0], refatom1B[i][1], refatom1B[i][2]);
		}
	}

void measure_helix_pair(int helix1, int helix2, int residue_index)
	{
	int i,j;
	if (setflag[flag_v] || setflag[flag_l])
		printf("- contacts between helix %3d v helix %3d\n\n",
			helix1,helix2);
	for (i = 0; i < residue_index; i++)
		if (helix_no[i] == helix1)
			for (j = 0; j < residue_index; j++)
				if (helix_no[j] == helix2)
					measure_residue_pair(i,j);
	if (setflag[flag_v] || setflag[flag_l])
		printf("- done contacts helix %3d v helix %3d\n\n",
			helix1,helix2);
	}

void measure_residue_pair(int res1, int res2)
	{
	float CA_distance, centre_distance, end_distance;
	if (setflag[flag_v] || setflag[flag_l])
		CA_distance = measure_CA_distance(res1,res2);
	centre_distance = measure_centre_distance(res1,res2);
	if (setflag[flag_v] || setflag[flag_l])
			end_distance = measure_end_distance_B(res1,res2);
	if (centre_distance < cutoff2)
		{
		if (setflag[flag_v] || setflag[flag_l])
			{
			printf("\n%4d (%s %5d:%c, iCode='%c', helix %3d),%4d (%s %5d:%c, iCode='%c', helix %3d): %8.3f, %8.3f, %8.3f \n",
				res1, helix_residue_name[res1], helix_residue_no[res1],
				helix_chain[helix_no[res1]], helix_residue_iCode[res1], helix_no[res1], res2,
				helix_residue_name[res2], helix_residue_no[res2],
				helix_chain[helix_no[res2]], helix_residue_iCode[res2], helix_no[res2], CA_distance,
				centre_distance, end_distance);

			measure_end_distance(res1,res2);
			}
		add_contact(res1,res2,centre_distance); add_contact(res2,res1,centre_distance);
		}
	}

void add_contact(int res1, int res2, float distance)
	{
	if (++n_contacts[res1] > 4)
		printf("%d sidechains in contact with residue %4d (%s %5d:%c, iCode='%c' helix %3d)\n",
			n_contacts[res1], res1, helix_residue_name[res1],
			helix_residue_no[res1], helix_chain[helix_no[res1]],
			helix_residue_iCode[res1], helix_no[res1]);
	
	contact[res1][n_contacts[res1]-1] = res2;
	contact_distance[res1][n_contacts[res1]-1] = distance;
	}

float measure_CA_distance(int res1, int res2)
	{
	return distance(coord[refatom0[res1]],coord[refatom0[res2]]);
	}

void measure_end_distance(int res1, int res2)
	{
	int i,j;
	for (i = 0; i < 2; i++)
		if (refatom1[i][res1] != -1)
			for (j = 0; j < 2; j++)
				if ((refatom1[j][res2] != -1) /* && (setflag[flag_v] || setflag[flag_l])*/ )
					printf("\t(%d,%d): %8.3f\n", i,j,distance(coord[refatom1[i][res1]],coord[refatom1[j][res2]]));
	}

float measure_end_distance_B(int res1, int res2)
	{
	return distance(refatom1B[res1],refatom1B[res2]);
	}

float measure_centre_distance(int res1, int res2)
	{
	if (null_refatom2(res1) || null_refatom2(res2))
		return 99999.9;
	else	return distance(refatom2[res1],refatom2[res2]);
	}

float measure_k_end_h_CA(int knobid)
	{
	int i;
	float d;
	d = 0;
	for (i = 0; i < 4; i++)
		{
		d += distance(refatom1B[knob[knobid]],coord[refatom0[hole[knobid][i]]]);
		}
	return d/4.0;
	}

float distance(float coords1[], float coords2[])
	{
	int i;
	float a;
	a = 0.0;
	for (i = 0; i < 3; i++) a += pow(coords1[i] - coords2[i],2);
	return sqrt(a);
	}

int null_refatom2(int resno)
	{ /* XXX change these hardcoded constants */
	/* if atom with serial no atomno has all 3 coordinates == 9999.99, then
	a true result is returned; otherwise false */
	int i,j;
	j = 1;
	for (i = 0; i < 3; i++) if (refatom2[resno][i] < 99999.0) { j = 0;  break;}
	return j;
	}
