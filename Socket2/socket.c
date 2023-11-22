/*

					SOCKET
					 v3.02

					socket.c

					17-10-01
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

/*					socket.c
					--------
*/

#include "socket.h"

#include "argument.h"

#include "aminoa2.h"

/* GNU c compiler complains if main() does not return type int */
int main(int argc, char *argv[])
	{
	int i,/*j,k, last_residue_index,*/ residue_index,
		 atom_index;

/*	char lastchain; */
/*	enum amino_acid ;*/

	printf(SOCKET_TITLE);

	statchar();

	initialize(argc,argv);

	check_files();

	if (extend)
		pre_parse_dssp(&helix_index, helix_start, helix_start_iCode,
		helix_end, helix_end_iCode, helix_chain, extend);

/* read in the alpha-helical residues from the DSSP file */
	residue_index = read_helical_dssp(extend);

	/* next line should now be redundant */
	/*if (extend)
		prune_extended_helices(residue_index, helix_index, helix_start, helix_end, helix_chain);*/

	atom_index = read_helical_pdb();

	determine_centre_of_mass(atom_index);
	determine_end(residue_index,atom_index);

	if (setflag[flag_v])
		{
		printf("\nThese are the reference coordinates for each residue:\n\n");
		for (i = 0; i < residue_index; i++)
			printf("\tresidue %d) end: %8.3f,%8.3f,%8.3f; centre of volume: %8.3f,%8.3f,%8.3f\n",
			i,refatom1B[i][0],refatom1B[i][1],refatom1B[i][2],refatom2[i][0],refatom2[i][1],refatom2[i][2]);
		printf("\n");
		}

	find_knobs_and_holes(residue_index);

	/*if ((par[par_d] != NULL) || (par[par_o] != NULL))
		write_files(residue_index);*/

	/* determine_order returns the number of coiled coils as specified
	by any helix-helix kih interactions */
	n_total_ccs = determine_order(residue_index);
	if (n_total_ccs)
		{
		/* find_register returns the number of the above coiled coils
		which are 'real', after weeding out those which dont have enough
		complementary knobs for example */
		n_true_ccs = find_register(residue_index);

		if (n_true_ccs)
			{
			if (par[par_r] != NULL)
				{
				define_ras_coils();
				fprintf(rasmol_file,"\nselect not coiled_coils\nstrands 1\nselect coiled_coils\nribbon 300\n");
				}
			printf("%s c %5.2f e %d result %d COILED COILS PRESENT", code, cutoff2, extend, n_true_ccs);
				if (n_total_ccs - n_true_ccs)
				printf(" (+ %d helix groups are either pairs with too few complementary knob in hole interactions or are subsets of larger coiled coils)\n",
				n_total_ccs - n_true_ccs);
				else printf("\n");
			}
		else printf("%s c %5.2f e %d result NO COILED COILS (but %d helix pairs have a single complementary knob in hole interaction)\n",
			code, cutoff2, extend, n_total_ccs);
		}

	else printf("%s c %5.2f e %d result NO COILED COILS\n",code, cutoff2, extend);

	/* the 'long' (-o) and 'summary' (-d) output files are created here, if specified;
		N.B. the rasmol script is more complicated and is written to by several
		different subroutines, including main() - see above - , 
		determine_order() and find_register() */

	if ((par[par_d] != NULL) || (par[par_o] != NULL))
		write_files(residue_index);


	printf("Finished\n");
	return(0);
	/* end of main() */
	}

void initialize(int argc, char *argv[])
	{
	int i,j,l;
	char temp[4];

	aa_map();

	set_flagsNpars(argc,argv);

	if (par[par_c] == NULL) cutoff2 = DEFAULT_CUTOFF2;
	else sscanf(par[par_c],"%f",&cutoff2);

	printf("using cutoff of %4.1f Angstroms for centre of mass distances\n\n",cutoff2);

	/* check for a user-specified refatom0 type (default is REFATOM0) */
	if (par[par_t] == NULL) strcpy(refatom0type,REFATOM0);
	else
		{
		strcpy(refatom0type,par[par_t]);
		/* the atom name field is 4 characters long:
			digit (or blank), letter, 
			letter (or blank), digit (or blank) */
		l = strlen(refatom0type);
		if (l > 4)
			{
			printf("\nYou specified an invalid reference atom type: %s\n\n",par[par_t]);
			exit(1);
			}
		if (l < 4)
			{
			strcpy(temp,refatom0type);
			strcpy(refatom0type,"");
			j = 0;
			for (i = 0; i < l; i++)
				{
				if (!(isalnum(temp[i]) || isspace(temp[i]))) {
				printf("\nYou specified an invalid reference atom type: %s\n\n",par[par_t]);
				exit(1);
				}
				if (temp[i] != ' ')
					{
					if (isalpha(temp[i]) && (j == 0)) refatom0type[j++] = ' ';
					refatom0type[j++] = temp[i];
					}
				if (j > 4) break;
				}
			while (strlen(refatom0type) < 4) strcat(refatom0type," ");
			}
		for (i = 0; i < 5; i++) refatom0type[i] = toupper(refatom0type[i]);
		}

	/*printf("Reference atom type is \"%s\"\n\n",refatom0type); */
	if (strcmp(refatom0type," CA ")) printf("- packing geometry cannot be calculated (requires C alpha atoms)\n");
	

	for (i = 0; i < MAX_RESIDUES; i++)
		{
		refatom0[i] = -1;
		refatom1[0][i] = -1;
		refatom1[1][i] = -1;
		for (j = 0; j < 3; j++) refatom2[i][j] = 99999.9;
		}

	for (i = 0; i < MAX_HELICES; i++) {n_knobs[i] = 0; n_hole_res[i] = 0;}

	if (setflag[flag_debug])
		{
		for (i = 0; i < FLAGS; i++)
			if (setflag[i]) printf("flag \"%s\" is set\n",flagname[i]);
		for (i = 0; i < PARS; i++)
			printf("par \"%s\" = \"%s\"\n",parname[i],par[i]);
		}

	if (par[par_k] == NULL) knob_threshold = DEFAULT_KNOB_THRESHOLD;
	else sscanf(par[par_k],"%i",&knob_threshold);

	if ((knob_threshold < 0) || (knob_threshold > 4)) {
		printf("minimum knob threshold must be between 0 and 4\n\n");
		exit(1);
		}
	cutoff_warning = ' ';

	if (par[par_e] == NULL) extend = 0;
	else sscanf(par[par_e],"%i",&extend);
	if (extend) printf("Helices defined by DSSP file will be extended by %d residues at each end\n",extend);


	for (i = 0; i < MAX_SEQCHAINS; i++)
		{
		seqchain[i] = NULL_SEQCHAIN;
		}

	for (i = 0; i < MAX_COILED_COILS; i++)
		{
		coiled_coil_tally[i] = 0;
		coiled_coil_helices[i] = 0;
		coiled_coil_subset[i] = -2;
		for (j = 0; j <  MAX_HELICES_PER_COIL; j++)
			{
			coiled_coil[i][j] = -1;
			for (l = 0; l < 2; l++)
				{
				coiled_coil_begin[i][j][l] = -1;
				coiled_coil_end[i][j][l] = -1;
				}
			}
		}

	} /* end of initialize */

void check_files()
	{
	char pdb_file_name[MAX_FILE_NAME_LENGTH],
		data_outfile_name[MAX_FILE_NAME_LENGTH],
		long_outfile_name[MAX_FILE_NAME_LENGTH],
		dssp_file_name[MAX_FILE_NAME_LENGTH],
		helix_register_file_name[MAX_FILE_NAME_LENGTH];

	if (par[par_f] == NULL)
		{
		printf("no PDB file was specified\n\n");
		exit(1);
		}

	strcpy(pdb_file_name,par[par_f]);
	if (!setflag[flag_p]) printf("attempting to open \"%s\"\n",pdb_file_name);
	if ( (pdb_file = fopen(pdb_file_name,"r") ) == NULL)
		{
		if (!setflag[flag_p]) printf("Failed to open PDB file \"%s\"\n\n",pdb_file_name);
		exit(1);
		}
	if (!setflag[flag_p]) printf("opened \"%s\" (input)\n",pdb_file_name);

	if (setflag[flag_p]) strcpy(code,"");
	else	get_code(pdb_file_name,code);

	if (par[par_s] == NULL)
		{
		printf("no secondary structure (DSSP) file was specified\n\n");
		exit(1);
		}
	strcpy(dssp_file_name,par[par_s]);

	if ( (dssp_file = fopen(dssp_file_name,"r") ) == NULL)
		{
		if (!setflag[flag_p]) printf("Failed to open DSSP file \"%s\"\n\n",dssp_file_name);
		exit(1);
		}
	if (!setflag[flag_p]) printf("opened \"%s\" (input)\n",dssp_file_name);

	if (par[par_r] != NULL)
		{
		if ( (rasmol_file = fopen(par[par_r],"w") ) == NULL)
			{
			if (!setflag[flag_p]) printf("Failed to open rasmol script file \"%s\"\n\n",par[par_r]);
			exit(1);
			}
			if (!setflag[flag_p]) printf("opened \"%s\" (output)\n",par[par_r]);
			fprintf(rasmol_file,"# RasMol script, created by ");
			fprintf(rasmol_file,SOCKET_TITLE);
			fprintf(rasmol_file,"echo\necho RasMol script, created by:\necho ");
			fprintf(rasmol_file,SOCKET_TITLE);
			fprintf(rasmol_file,"echo\nselect\nwireframe off\nspacefill off\nribbon off\ndots off\nbackbone off\nstrands off\n");
		}

	if (par[par_d] != NULL)
		{  
		strcpy(data_outfile_name,par[par_d]);
		if ( (data_outfile = fopen(data_outfile_name,"w") ) == NULL)
			{
			if (!setflag[flag_p]) printf("Failed to open output file \"%s\"\n\n",data_outfile_name);
			exit(1);
			}
		}

	if (par[par_o] != NULL)
		{  
		strcpy(long_outfile_name,par[par_o]);

		if ( (long_outfile = fopen(long_outfile_name,"w") ) == NULL)
			{
			if (!setflag[flag_p]) printf("Failed to open output file \"%s\"\n\n",long_outfile_name);
			exit(1);
			}
		if (!setflag[flag_p]) printf("opened \"%s\" (output)\n",long_outfile_name);
		}

	if (par[par_w] != NULL)
		{  
		strcpy(helix_register_file_name,par[par_w]);

		if ( (helix_register_file = fopen(helix_register_file_name,"w") ) == NULL)
			{
			if (!setflag[flag_p]) printf("Failed to open helix register file \"%s\"\n\n",helix_register_file_name);
			exit(1);
			}
		if (!setflag[flag_p]) printf("opened \"%s\" (output)\n",helix_register_file_name);
		}

	
	}





void write_files(int residue_index)
	{
	int i,j,k,l,m,n;
	char knob_pattern[3][MAX_KNOBS_PER_HELIX*2],s[20];
	static char* orientation_id = {"pa"};

	if (par[par_d] != NULL) fprintf(data_outfile,SOCKET_TITLE);
	if (par[par_o] != NULL) fprintf(long_outfile,SOCKET_TITLE);


	for (i = 0; i < helix_index; i++)
		{

		for (j = 0; j < 3; j++) strcpy(knob_pattern[j],"");

		if (par[par_d] != NULL)
			{
			fprintf(data_outfile,"%s\t%d (%c)\t%d(iCode='%c')..%d(iCode='%c')\t%4.1f%c\t%d",
				par[par_f],i,helix_chain[i],
				helix_start[i],helix_start_iCode[i],
				helix_end[i],helix_end_iCode[i],
				cutoff2,cutoff_warning,n_knobs[i]);

			for (j = 0; j < 7; j++) fprintf(data_outfile,",%d",n_knobtype[i][j]);
			}

		if (par[par_o] != NULL)
			{
			fprintf(long_outfile,"%s helix\t%d (chain %c)\t%d(iCode='%c')..%d(iCode='%c')\tcutoff %4.1f%c\t%d knobs",
				code,i,helix_chain[i],
				helix_start[i],helix_start_iCode[i],
				helix_end[i],helix_end_iCode[i],
				cutoff2,cutoff_warning,n_knobs[i]);
			for (j = 0; j < 7; j++) fprintf(long_outfile,", %d type %d",n_knobtype[i][j],j);
			fprintf(long_outfile,"\n");
			}

		l = 0;

		for (j = 0; j < residue_index; j++)
			{
			if (helix_no[j] == i) 
				{
				if (par[par_o] != NULL)
					fprintf(long_outfile,"%s %s%5d:%c iCode='%c'",code,helix_residue_name[j],helix_residue_no[j],
					helix_chain[helix_no[j]],helix_residue_iCode[j]);

				/* print out ALL the register assignments for this residue (there can be more than one; some
				residues can simultaneously belong to 2 coiled coils) */

				for (n = 0; n < MAX_COILED_COILS; n++)
					{
					if ((coiled_coil_subset[n] == -1) && 
						(tad_register[j][n] != ' ') && (tad_register[j][n] != 0))
						fprintf(long_outfile,"R%c[%d%c]",tad_register[j][n],
							coiled_coil_helices[n],orientation_id[coiled_coil_orientation[n]]);
					}

				if (l == MAX_KNOBS_PER_HELIX)
					{
					printf("Too many knobs (%d) in this helix (%d)\n",l,i);
					exit(1);
					}
				k = 0;

				while ((knob[k] != j) && (k < knob_index)) k++;

				if (k == knob_index)
					{
					knob_pattern[0][l] = '-';
					/*knob_pattern[1][l] = '-';*/
					}
				else
					{
					knob_pattern[0][l] = '0' + knobtype[k];
					/*knob_pattern[1][l] = '-';*/

					if (par[par_o] != NULL) {
						/* prints the knobtype, helix id and packing angle */
						fprintf(long_outfile," T%1d H%3d:%8.3f ; ",
							knobtype[k],helix_no[hole[k][0]],angle[k]);
						/* prints the list of 4 hole residues to the file */

						for (n = 0; n < MAX_COILED_COILS; n++)
							{
							if ((coiled_coil_subset[n] == -1) &&
								(tad_register[hole[k][1]][n] != ' ') && (tad_register[j][n] != 0))
								{
								fprintf(long_outfile,"hole (");
								for (m = 0; m < 4; m++)
									fprintf(long_outfile,"%c",tad_register[hole[k][m]][n]);
								fprintf(long_outfile,") ");
								}
							}

						fprintf(long_outfile,"chain %c: ",helix_chain[helix_no[hole[k][0]]]);
                                                /* JW 13-7-6 amended print statement so that it prints the CoV separation
						   of the knob and each of the 4 hole sidechains */
						for (m = 0; m < 4; m++)
							fprintf(long_outfile," (%d) %s%5d'%c' %8.3f",m,
								helix_residue_name[hole[k][m]],
								helix_residue_no[hole[k][m]],
								helix_residue_iCode[hole[k][m]],
								hole_distance[k][m]);
						/* calculate and print out the hole dimensions, ie the lengths of its four
							sides: distances h0-h1 , h0-h2, h1-h3, h2-h3 , where h0, h1 , h2, h3 are the
							centres of volume respectively of the four hole residues in serial
							order (usually, if h0 is residue x, then h1 is x+3, h2 is x+4, h3 is x+7) */
						fprintf(long_outfile,"; sides 0-1:%8.3f, 0-2:%8.3f, 1-3:%8.3f, 2-3:%8.3f",
							measure_centre_distance(hole[k][0],hole[k][1]),
							measure_centre_distance(hole[k][0],hole[k][2]),
							measure_centre_distance(hole[k][1],hole[k][3]),
							measure_centre_distance(hole[k][2],hole[k][3]));
						}
					}

				for (k = 0; k < knob_index; k++)
					{
					/*if (strcmp(knob_pattern[2],"")) strcat(knob_pattern[2],",");
					strcat(knob_pattern[2],"-");*/
					if (knob[k] == j)
						{
						if (strcmp(knob_pattern[2],"")) strcat(knob_pattern[2],",");
						sprintf(s,"%d",helix_no[hole[k][0]]);
						strcat(knob_pattern[2],s);
						}
					}

				if (par[par_o] != NULL) fprintf(long_outfile,"\n");
				l++;
/*knob_pattern[0][l] = '\0';
printf("checkpoint 4 i = %d, j = %d, l = %d, \"%s\", \"%s\"\n",i,j,l,knob_pattern[0],knob_pattern[2]);
*/
				}
			}

		for (j = 0; j < 2; j++) knob_pattern[j][l] = '\0';
/*		if (par[par_d] != NULL) fprintf(data_outfile,"\n");*/
		if (par[par_d] != NULL)
			{
			/* prints the type of each knob */
			fprintf(data_outfile," %s ",knob_pattern[0]);
			fprintf(data_outfile," %s ",knob_pattern[2]);

			/* prints the helix containing the hole in which each knob fits */
			/*for (j = 0; j < n_knobs[i]; j++)
				fprintf(data_outfile,"\t%d",helix_no[hole[j][0]]);*/
			fprintf(data_outfile,"\n");
			}
		}
	}

/* function get_code: filename[] is input (the name of a user-specified file,
which may or may not include a path specification) and code[] is effectively
the result of this function, ie a 4-letter PDB-style code. Such a 4-letter
code is returned only if the file name (ignoring the path) is of the form
either: \d\w{3} of pdb\d\w{3}
- otherwise, the full file-name (minus path) is returned, which is treated
as a pseudo-PDB 'code' */

void get_code(char filename[], char code[])
	{
	/* i marks the position in a string;
	j marks the position in a string, relative to i;
	k is a boolean (ought to be redundant) specifying
		whether a character is alphanumeric;
	l holds the length of various strings;
	m is a boolean specifying that a 4-letter string
		 matches the regular expression \d\w{3} */
	int i,j,k,l,m;

	char no_path[MAX_FILE_NAME_LENGTH], tmpstr[5];

	l = strlen(filename);

	/* ignore anything in the file name before the '/', if there is one */
	i = l;
	while ((i >= 0) && (filename[i] != '/')) i--;
	j = 0;
	while (i++ < l) no_path[j++] = filename[i];

	/* at this point no_path holds the name of the file, minus any path */

	/* copy no_path into code, translating all upper to lower case */
	for (i = 0; i < strlen(no_path); i++)
		{
		if ((no_path[i] >= 'A') && (no_path[i] <= 'Z'))
			code[i] = no_path[i] -'A' + 'a';
		else code[i] = no_path[i];
		}

	/*find the position of a 4 letter \d\w{3} code (if any) */
	i = -1;
	m = 0;
	l = strlen(code);
	while ((!m) && (++i < l - 3) )
		{
		if ((code[i] >= '0') && (code[i] <= '9'))
			{
			m = 1;
			k = 0;
			for (j = 1; j < 4; j++)
				{
				if ((code[i+j] >= '0') && (code[i+j] <= '9')) k = 1; /* messy */
 				else if ((code[i+j] >= 'a') && (code[i+j] <= 'z')) k = 1;
				if (!k) m = 0;
				}
			/* if now m == 1, its a \d\w{3} code starting at position i */
			}
		}

	if ((m) && ((!i) || (i == 3)) )
		{
		substring(code,0,3,tmpstr);

		if ((i == 3) && (strcmp(tmpstr,"pdb"))) strcpy(code,no_path);
		else	{
			substring(code,i+4,4,tmpstr);

			if ((strcmp(tmpstr,".ent")) && (strcmp(tmpstr,".pdb"))) strcpy(code,no_path);
			else {substring(code,i,4,code); 
				printf("file name implies standard PDB entry, code %s\n",code);}
			}
		}
	else strcpy(code,no_path);

	}

void substring(char parent[], int start, int length, char sub[])
	{
	int i,l; char tmpstr[132];
	l = strlen(parent);

	/* uses temporary string tmpstr in case parent and sub are the same string */
	strcpy(tmpstr,"");
	for (i = 0; i < length; i++)
		if (i < l) tmpstr[i] = parent[start+i];
	if (start + i > l) tmpstr[start + i - l] = '\0';
	else tmpstr[i] = '\0';
	strcpy(sub,tmpstr);
	}



void dumpknobs(int aknob)
	{
	int i,j;
	for (i = 0; i < knob_index; i++)
		{
		if ((i == aknob) || (aknob == -1))
		printf("knob %3d: knob[%d]=%3d hole[%d][0]=%3d hole[%d][1]=%3d hole[%d][2]=%3d hole[%d][3]=%3d knobtype[%d]=%d knob_order[%d]=%d angle[%d]=%8.3f n_compknob[%d]=%d ",
i,i,knob[i],i,hole[i][0],i,hole[i][1],i,hole[i][2],i,hole[i][3],i,knobtype[i],i,knob_order[i],i,angle[i],i,n_compknob[i]);
	for (j = 0; j < n_compknob[i]; j++) printf("compknob[%d][%d] = %d; ",i,j,compknob[i][j]);
		printf("\n");
		}
	}
