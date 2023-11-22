/*

					SOCKET
					v3.02

					read.c

					25-10-01
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

/*					read.c
					------

5 functions:

	int	read_helical_dssp(int extend)

	int	within_helix(int resno, char iCode, char chain, int helix_start[], char helix_start_iCode[],
		int helix_end[], char helix_end_iCode[], char helix_chain[], int n_helices)

	void	pre_parse_dssp(int *helix_index, int helix_start[], char helix_start_iCode[],
		int helix_end[], char helix_end_iCode[], char helix_chain[], int extend)

	void	prune_extended_helices(int residue_index, int helix_index, int helix_start[], int helix_end[], char helix_chain[])
		(should be redundant)

	int	read_helical_pdb()
*/

#include "socket.h"
#include <stdlib.h>
int read_helical_dssp(int extend /*, int join*/)
	{
	int i,/*j,*/ lastresidue, helical_residues, residue_index, helix_id;
	char chainID,/*altLoc,*/ iCode,lastiCode,ch,aacode,lastchainID;

/* read in the alpha-helical residues from the DSSP file */

	if (setflag[flag_debug])
		{
		printf("\nReading DSSP");
		if (extend) { printf("(pre-parsed helix_index = %d)",helix_index); }
		printf("\n");
		}

	/* first find the beginning of the residue data */
	while ((fgets(record_type,7,dssp_file) != NULL) && (strcmp(record_type,"  #  R")))
		{
		if (setflag[flag_debug]) printf("\"%s\"\n",record_type);
		fgets(textstring,MAX_LINE_WIDTH,dssp_file);
		if (setflag[flag_debug])
			{
			printf("DSSP line: \"%s\"\n",textstring);
			}
		}

	if (strcmp(record_type,"  #  R") == 0) printf("Found beginning of residue data\n");
	fgets(textstring,MAX_LINE_WIDTH,dssp_file);
	if (setflag[flag_debug]) printf("DSSP line: \"%s\"\n",textstring);

	residue_index = 0;
	if (!extend) helix_index = 0;
	is_helical = false;

/*for (i = 0; i < helix_index; i++)
			printf("\thelix %d (%d-%d:%c)\n",i,helix_start[i],helix_end[i],helix_chain[i]);XXX*/

	while (fgets(textstring,18,dssp_file) != NULL)
		{
		if (setflag[flag_debug]) printf("DSSP line: \"%s\"",textstring);
		if ((textstring[13] == '!') || (textstring[14] == '*'))
			{fgets(textstring,MAX_LINE_WIDTH,dssp_file); printf("\n"); continue;}
		sscanf(textstring,"%*d%5d%c%c%*c%c%*c%*c%c",&i,&iCode,&chainID,&aacode,&ch);
		if (setflag[flag_debug]) printf("\tread residue data for %d %c, iCode='%c'\n",i,chainID,iCode);

		/* check for out-of-sequence residues */

		if (	(chainID == lastchainID) &&
					(	(i < lastresidue)	||
						(	(i == lastresidue) && (iCode < lastiCode)	)
					)
			)
			printf("!!! WARNING: NON-INCREMENTAL SEQUENCE: %d:%c, iCode='%c' PRECEDES %d:%c, iCode='%c'\n",
				lastresidue,lastchainID,lastiCode,i,chainID,iCode);

		 /* if extend is non-zero, then the helix positions will have already been determined;
		    so check that the current residue is in one of them */
		if (extend) helix_id =
			within_helix(i,iCode,chainID,helix_start,helix_start_iCode,
							helix_end,helix_end_iCode,helix_chain,helix_index);

		if ((ch == 'H') || (extend && (helix_id != -1)) )
			{
			if (!is_helical)
				{
				if (!extend)
					{
					if (helix_index == MAX_HELICES)
						{
						printf("maximum number of helices (%d) exceeded\n",MAX_HELICES);
						exit(1);
						}

					/* this is *helix-specific* data
						(this residue is the first in the helix) */
					helix_start[helix_index] = i;
					helix_start_iCode[helix_index] = iCode;
					helix_chain[helix_index] = chainID;
					if (setflag[flag_debug]) printf("start of new helix (%d)\n",helix_index);
					}
				else	{if (setflag[flag_debug]) printf("start of new helix (%d)\n",helix_id);}
				is_helical = true;
				}
			if (residue_index == MAX_RESIDUES)
				{
				printf("Maximum number of alpha-helix residues (%d) exceeded\n",MAX_RESIDUES);
				exit(1);
				}

			/* this is *residue-specific* data */

			if (extend) helix_no[residue_index] = helix_id;
				else helix_no[residue_index] = helix_index;
			helix_residue_no[residue_index] = i;
			helix_residue_iCode[residue_index] = iCode;

			/* a lower case aacode indicates a Cystine bridge- DSSP labels each bridge
				starting from 'a' */
			if islower(aacode) aacode = 'C';

			helix_residue_aacode[residue_index] = map_alpha_to_amino_acid[aacode-'A'];
			lastresidue = i;
			lastiCode = iCode;
			lastchainID = chainID; /* so that out-of-sequence residues are spotted */
			if (setflag[flag_debug]) printf("residue %d: helix_no=\t%d; helix_residue_no=\t%d; helix_residue_iCode=\t%c; helix_residue_aacode=\t%d\n",
				residue_index,helix_no[residue_index],helix_residue_no[residue_index],helix_residue_iCode[residue_index],helix_residue_aacode[residue_index]);
			residue_index++;
			}
		else if (is_helical)
			{
			if (!extend)
				{
				helix_end[helix_index] = lastresidue;
				helix_end_iCode[helix_index] = lastiCode;
				if (setflag[flag_debug]) printf("end of helix (%d)\n",helix_index);
				helix_index++;
				}
			else if (setflag[flag_debug]) printf("end of helix\n");
			is_helical = false;
			}

		fgets(textstring,MAX_LINE_WIDTH,dssp_file);
		}

	helical_residues = residue_index;

	printf("There are %d alpha-helical residues in this structure\n\n",helical_residues);

	if (!helical_residues) {printf("%s c %5.2f e %d result NO COILED COILS\nFinished\n",code,cutoff2,extend); exit(1);}

	if (setflag[flag_v])	printf("These are the %d helices:\n\n\thlx# res iCd  res iCd\tch\n\n",helix_index);

	for (i = 0; i < helix_index; i++)
		{
		if (setflag[flag_v]) printf("\t%3d) %4d %c - %4d %c\t%c\n",
			i,helix_start[i],helix_start_iCode[i],
			helix_end[i],helix_end_iCode[i],helix_chain[i]);
		if (par[par_r] != NULL)
			fprintf(rasmol_file,"define helix%d %d-%d:%c\ndefine h%d %d-%d:%c\n",
			i,helix_start[i],helix_end[i],helix_chain[i],i,helix_start[i],helix_end[i],helix_chain[i]);
		}


	return residue_index;
	/* end of read_helical_dssp */
	}

int within_helix(int resno, char iCode, char chain, int helix_start[], char helix_start_iCode[],
	int helix_end[], char helix_end_iCode[], char helix_chain[], int n_helices)
	{

	int h, helix_match;
	helix_match = -1;
	h = 0;

	if (setflag[flag_debug])
		{
		printf("checking query residue: %d:%c, iCode='%c' :\n", resno, chain, iCode);
		}

	if (setflag[flag_debug])
		{
		printf("\t\tversus helix %d (chain %c, %d[iCode='%c']..%d[iCode='%c'])\n",
				h,helix_chain[h],helix_start[h],helix_start_iCode[h],
				helix_end[h],helix_end_iCode[h]);
		}

	while ( ((chain != helix_chain[h]) || (resno < helix_start[h]) || (resno > helix_end[h]) )
		&& (h < n_helices -1) )
		{
		h++;
		if (setflag[flag_debug])
			{
			printf("\t\tversus helix %d (chain %c, %d[iCode='%c']..%d[iCode='%c'])\n",
					h,helix_chain[h],helix_start[h],helix_start_iCode[h],
					helix_end[h],helix_end_iCode[h]);
			}
		}

	if ((h < n_helices) && (chain == helix_chain[h])
		&& (resno >= helix_start[h]) && (resno <= helix_end[h]))
		{
		/* the residue is within the bounds of the helix (at least, numerically in terms
		of the residue number; see below), and is in the same chain;
		however, we're still not quite there; complications arise where residues have
		non-null iCodes. This is only an issue if the query residue has the SAME residue
		number as either the first or last residue of the helix; so lets check it: */

		if (setflag[flag_debug])
			{
			printf("\t\tversus helix %d (chain %c, %d[iCode='%c']..%d[iCode='%c']) MATCH\n",
					h,helix_chain[h],helix_start[h],helix_start_iCode[h],
					helix_end[h],helix_end_iCode[h]);
			}

		if (resno == helix_start[h])
			{
			/*
			The following assumes that the query residue and the N-terminal residue of the
			helix have the same residue number.

			We assume the following:
			in the PDB file, for a given residue number i, the residue (if any) with a null
			iCode always precedes residue(s) with non-null iCode(s) (if any). Also, assume
			that the insertion codes are always in an ascending sequence.
			Example sequence:

				250
				250A
				250B
				251
				 ^ ^
   	              	 i iCode

			- so, if the query residue and the start of the helix are identical in terms
			of residue number and insertion code, then the residue is within the helix;
			if they are non-identical, then the residue is within the helix only if the
			query residue's iCode is 'greater' (in ASCII terms) than that of the N-terminus.
			(a space ' ' has the lowest ASCII code of any non-control character);
			UNLESS this is an unusual case where the START and END of the helix ARE THE
			SAME (these do occur, e.g. 1a2c); in which case, an extra check must be made-
			the query residue's iCode cannot in that case be greater than that of the
			end residue */

			/* e.g. 1a2c in more detail:

				this is part of the DSSP file. The helix spans 14C..14G of chain L.
				If 14H of chain L is compared v 14C, then the chain and residue number
				are both identical. 'H' > 'C', which would conclude that 14H:L is in
				the helix (wrong) - unless the extra check is made. 14A and 14B will
				never be compared with the end of the helix, so there will not be a
				failure in the 'reverse' case (they will always get compared with
				the start first, which will match, so the 'if (resno == helix_end[h])'
				never gets tried

			   20   12 L L      < -
			   21   13 L E        -
			   22   14 L D        -
			   23   14AL K  S    S+
			   24   14BL T  S >> S+
			   25   14CL E  H 3> S+
			   26   14DL R  H 3> S+
			   27   14EL E  H <4 S+
			   28   14FL L  H >X S+
			   29   14GL L  H ><>S+
			   30   14HL E  T 3<5S+
			   31   14IL S  T <45S+
			   32   14JL Y  T <<5S-
			   33   14KL I  T   5S-
			   34   14LL D  S   <S-
			   35   14ML G
			   36   15 L R
			   37        !*
			*/

			if (setflag[flag_debug])
				{
				printf("residue number %d is same as start of chain; comparing iCodes ('%c' v '%c')\n",
						resno,iCode,helix_start_iCode[h]);
				}


			if (	(iCode >= helix_start_iCode[h]) /* this alone suffices in most cases */
				&&	(	(helix_end[h] != helix_start[h]) /* see the above re unusual cases */
					||	(iCode <= helix_end_iCode[h])
					)
				)
				helix_match = h;
			}

		else if (resno == helix_end[h])
			{
			/*
			The following assumes that the query residue and the C-terminal residue of the
			helix have the same residue number.

			- see example above

			- so, if the query residue and the end of the helix are identical in terms
			of residue number and insertion code, then the residue is within the helix;
			if they are non-identical, then the residue is within the helix only if the
			query residue's iCode is 'less' (in ASCII terms) than that of the C-terminus.
			(a space ' ' has the lowest ASCII code of any non-control character). */

			if (setflag[flag_debug])
				{
				printf("residue number %d is same as end of chain; comparing iCodes ('%c' v '%c')\n",
						resno,iCode,helix_end_iCode[h]);
				}

			if (iCode <= helix_end_iCode[h])	helix_match = h;
			}

		else
			{
			/* the query residue equals neither the N-terminal nor C-terminal residue
			of the helix; so it must be somewhere in between, in which case we don't
			care what the insertion codes are, if any */
			helix_match = h;
			}


		if ((helix_match != -1) && (setflag[flag_v] || setflag[flag_debug]))
			printf("%d:%c iCode='%c' lies between %d:%c iCode='%c' and %d:%c iCode='%c' (helix %d)\n",
				resno,chain,iCode,
				helix_start[h],helix_chain[h],helix_start_iCode[h],
				helix_end[h],helix_chain[h],helix_end_iCode[h],h);
		}

	if ((helix_match == -1)  && (setflag[flag_debug]))
		{
		printf("NOT IN A HELIX\n");
		}

	return helix_match;
	}



void pre_parse_dssp(int *helix_index, int helix_start[], char helix_start_iCode[],
	int helix_end[], char helix_end_iCode[], char helix_chain[], int extend)
	{
	int h,i,j,r, merged, lastresidue, all_index, residue_index, chain_start[63],
	chain_end[63],

		/* linked list of residues (all of them, irrespective of secondary
		structure) - needed to perform helix-extension, because you can't
		rely on uniform sequential numbering of residues; deletions or
		insertions (which will have non-null iCode fields) cause problems, so
		these 3 arrays keep track of whats next to what */

		all_residue_i[MAX_RESIDUES_ALL],
		helix_start_index[MAX_RESIDUES], /* index for a helical residue, specifying */
		helix_end_index[MAX_RESIDUES]; /* index for a helical residue, specifying */

	char	all_residue_iCode[MAX_RESIDUES_ALL],
		all_residue_chainID[MAX_RESIDUES_ALL],

		chainID, lastchainID, iCode, lastiCode, ch, aacode;


	/* initialize the chain_start and chain_end arrays */
	for (i = 0; i < 63; i++)
		{
		chain_start[i] = 9999; chain_end[i] = -9999;
		}


	if (setflag[flag_debug] || setflag[flag_l] || setflag[flag_v]) printf("Pre-parsing DSSP file\n");
/* read in the alpha-helical residues from the DSSP file */
	/* first find the beginning of the residue data */
	while ((fgets(record_type,7,dssp_file) != NULL) && (strcmp(record_type,"  #  R")))
		{
		if (setflag[flag_debug]) printf("\"%s\"\n",record_type);
		fgets(textstring,MAX_LINE_WIDTH,dssp_file);
		if (setflag[flag_debug])
			{
			printf("DSSP line: \"%s\"\n",textstring);
			}
		}

	if (strcmp(record_type,"  #  R") == 0) printf("Found beginning of residue data\n");
	fgets(textstring,MAX_LINE_WIDTH,dssp_file);
	if (setflag[flag_debug]) printf("DSSP line: \"%s\"\n",textstring);

	residue_index = 0; /* the number of *alpha-helix* residues */
	all_index = 0;	/* the number of residues */
	*helix_index = 0;
	is_helical = false;

	while (fgets(textstring,18,dssp_file) != NULL)
		{
		if (setflag[flag_debug]) printf("DSSP line: \"%s\"",textstring);
		if ((textstring[13] == '!') || (textstring[14] == '*'))
			{fgets(textstring,MAX_LINE_WIDTH,dssp_file); printf("\n"); continue;}
		sscanf(textstring,"%*d%5d%c%c%*c%c%*c%*c%c",&i,&iCode,&chainID,&aacode,&ch);
		if (setflag[flag_debug]) printf("\tread residue data for %d %c %c\n",i,iCode,chainID);

		/* check for out-of-sequence residues */

		if (	(chainID == lastchainID) &&
					(	(i < lastresidue)	||
						(	(i == lastresidue) && (iCode < lastiCode)	)
					)
			)
			printf("!!! WARNING: NON-INCREMENTAL SEQUENCE: %d:%c, iCode='%c' PRECEDES %d:%c, iCode='%c'\n",
				lastresidue,lastchainID,lastiCode,i,chainID,iCode);

		/* store attributes of this residue */

		all_residue_i[all_index] = i;				/* residue number (as in PDB)*/
		all_residue_iCode[all_index] = iCode;		/* residue insertion code */
		all_residue_chainID[all_index] = chainID;	/* chain identifier */

		/* if the helices are to be extended, the first and last residues of each chain must
			be known, to specify limits of the extension */

		if ( (chainID != ' ') && ((chainID < 'A') || (chainID > 'Z')) &&
			((chainID < 'a') || (chainID > 'z')) &&
			((chainID < '0') || (chainID > '9')))
			{ printf("PDB file has unexpected chain identifier: '%c':\n%s\n",
			chainID,textstring); exit(1);}

		if (chainID == ' ') j = 0;
		else if ((chainID >= 'A') && (chainID <= 'Z')) j = chainID - 'A' + 1;
		else if ((chainID >= 'a') && (chainID <= 'z')) j = chainID - 'a' + 27;
		else j = chainID - '0' + 53;


		if (i < chain_start[j]) chain_start[j] = i;
		else if (i > chain_end[j]) chain_end[j] = i;


		/* is it a helical (DSSP) residue? */

		if ((ch == 'H'))
			{
			if (!is_helical)
				{
				is_helical = true;
				if (setflag[flag_debug]) printf("start of new helix (%d)\n",*helix_index);
				/*helix_start[*helix_index] = i;
				helix_start_iCode[*helix_index] = iCode;
				helix_chain[*helix_index] = chainID;*/
				helix_start_index[*helix_index] = all_index; /* points to residue in master list */
				}
			if (residue_index == MAX_RESIDUES)
				{
				printf("Maximum number of alpha-helix residues (%d) exceeded\n",MAX_RESIDUES);
				exit(1);
				}
			lastresidue = i;
			lastiCode = iCode;
			lastchainID = chainID; /* so that out-of-sequence residues are spotted */
			}
		else if (is_helical)
			{
			if (setflag[flag_debug]) printf("end of helix (%d)\n",*helix_index);
			helix_end_index[*helix_index] = all_index - 1; /* points to residue in master list */
			/*helix_end_iCode[*helix_index] = lastiCode;
			helix_end[(*helix_index)++] = lastresidue;*/
			(*helix_index)++;
			is_helical = false;
			}

		all_index++;

		fgets(textstring,MAX_LINE_WIDTH,dssp_file);
		}
	if (setflag[flag_debug] || setflag[flag_l] || setflag[flag_v]) printf("Pre-parsing DSSP completed\n");

	/* the start[] and end[] of each helix (as defined by the DSSP file) has now
	   been determined; these are now extended by extend residues at each end */

	for (h = 0; h < *helix_index; h++)
		{
		if (helix_chain[h] == ' ') j = 0;
		else if ((helix_chain[h] >= 'A') && (helix_chain[h] <= 'Z')) j = helix_chain[h] - 'A' + 1;
		else if ((helix_chain[h] >= 'a') && (helix_chain[h] <= 'z')) j = helix_chain[h] - 'a' + 27;
		else j = helix_chain[h] - '0' + 53;


		/* extend the N-terminal end of the helix, 1 residue at a time; it is done
		this way, to go as far to the end of the same chain as possible */

		for (r = 0; r < extend; r++)
			{
			/* check that there is a residue r-residues N-terminal to the first
				residue of this helix h */
			if (helix_start_index[h] > 0)
				/* check that the next residue is in the same chain */
				if (all_residue_chainID[helix_start_index[h] - 1] ==
					all_residue_chainID[helix_start_index[h]])
					{
					/* update info on the first residue of helix h */
					helix_start_index[h]--;
					}
			}

/*printf("helix_start[1] = %d;\n",helix_start[1]);
printf("helix_start[%d] = %d; chain_start[%d] = %d\n",h,helix_start[h],j,chain_start[j]);*/

		/* this check ought now to be redundant */
		/*if (helix_start[h] < chain_start[j]) helix_start[h] = chain_start[j];*/

/* DO THE EXTENDING THE SAME WAY AS ABOVE */

		/* extend the C-terminal end of the helix, 1 residue at a time; it is done
		this way, to go as far to the end of the same chain as possible */

		for (r = 0; r < extend; r++)
			{
			/* check that there is a residue r-residues N-terminal to the first
				residue of this helix h */
			if (helix_end_index[h] < all_index - 1)
				/* check that the next residue is in the same chain */
				if (all_residue_chainID[helix_end_index[h] + 1] ==
					all_residue_chainID[helix_end_index[h]])
					{
					/* update info on the first residue of helix h */
					helix_end_index[h]++;
					}
			}


		/*if (setflag[flag_v]) printf("helix_start[%d] = %d; helix_end[%d] = %d\n",
				h,helix_start[h],h,helix_end[h]);*/
		}


	/* which may mean that some of them merge, necessitating a new edited list of helices */

	/* merged specifies how many helix-pairs have been joined */
	merged = 0;

	/* h is the index of the helix */
	h = 0;
	/* i is the residue number */
	while (h < *helix_index -1)
		{
		/* need only compare each helix with the subsequent helix in the list;
		specifically, C-terminus of helix h with N-terminus of helix h+1;
		only merge them if they are overlapping/touching AND part of the same
		chain, of course */
		if ( (all_residue_chainID[helix_end_index[h]] ==
			  all_residue_chainID[helix_start_index[h+1]]) &&
			 (helix_start_index[h+1] - helix_end_index[h] < 2) )
			{
			if (setflag[flag_q] == false)
				printf("joining helices %d (%d%c-%d%c:%c) and %d (%d%c-%d%c:%c)\n",
					h,
					all_residue_i[helix_start_index[h]],
					all_residue_iCode[helix_start_index[h]],
					all_residue_i[helix_end_index[h]],
					all_residue_iCode[helix_end_index[h]],
					all_residue_chainID[helix_start_index[h]],
					h+1,
					all_residue_i[helix_start_index[h+1]],
					all_residue_iCode[helix_start_index[h+1]],
					all_residue_i[helix_end_index[h+1]],
					all_residue_iCode[helix_end_index[h+1]],
					all_residue_chainID[helix_start_index[h+1]]);
				helix_end_index[h] = helix_end_index[h+1];
			merged++;
			for (j = h+1; j < *helix_index - 1; j++)
				{
				helix_start_index[j] = helix_start_index[j+1];
				helix_end_index[j] = helix_end_index[j+1];
				}
			(*helix_index)--;
			}
		else h++;
		}
	if (merged && (setflag[flag_q] == false)) {
		printf("\nThe new list of helices:\n\n\thlx# res iCd  res iCd\tch\n\n");
		for (h = 0; h < *helix_index; h++)
			printf("\t%3d) %4d %c - %4d %c\t%c\n",
				h,
				all_residue_i[helix_start_index[h]],
				all_residue_iCode[helix_start_index[h]],
				all_residue_i[helix_end_index[h]],
				all_residue_iCode[helix_end_index[h]],
				all_residue_chainID[helix_start_index[h]]);
		printf("\n");
		}

	/* the first and last residues (in the all_residue_.. list) have now
 	been determined for each helix; so, set the global variables specifying
	these N- and C-termini */

	for (h = 0; h < *helix_index; h++)
		{
		helix_start[h] = all_residue_i[helix_start_index[h]];
		helix_start_iCode[h] = all_residue_iCode[helix_start_index[h]];
		helix_end[h] = all_residue_i[helix_end_index[h]];
		helix_end_iCode[h] = all_residue_iCode[helix_end_index[h]];
		helix_chain[h] = all_residue_chainID[helix_start_index[h]];
		}



	rewind(dssp_file);
	} /* end of pre_parse dssp */


/* this function should now be redundant XXX */
void prune_extended_helices(int residue_index, int helix_index, int helix_start[], int helix_end[], char helix_chain[])
	{
	int found,h,r;
	/* If the helices have been extended, it is possible that either the N-terminal of the first helix, or the C-terminal of the
	last, goes beyond the beginning or end of the amino acid sequence. Eg if the sequence went from 1-100, and the first and
	last helices (unextended) were 2-17 and 85-99, then an extension of two residues would cause both to go over the end.
	helix_start[0] would be set to zero for example, but there would be no residue with this actual sequence number. This will
	cause problems later on (ie the residue corresponding to the helix_start and helix_end of a particular helix h not actually
	existing. This function checks that all the residues at the beginning/end of each helix actually exist, and if they don't,
	then helix_start[h] is increased, or helix_end[h] decreased, until they do. */
	for (h = 0; h < helix_index; h++)
		{
		found = 0;
		while (!found)
			{
			r = 0;
			while ((r < residue_index) &&
				((helix_residue_no[r] != helix_start[h]) || (helix_chain[helix_no[r]] != helix_chain[h])))
				{
				printf("residue (%d): %d:%c, iCode='%c' v %d:%c\n",
					r, helix_residue_no[r], helix_chain[helix_no[r]], helix_residue_iCode[r],
					helix_start[h], helix_chain[h]);
				r++;
				}
			if ((helix_residue_no[r] == helix_start[h]) && (helix_chain[helix_no[r]] == helix_chain[h]))
				found = 1;
			else
				{
				printf("pruning start of helix %d from %d:%c to %d:%c\n",h,helix_start[h],
					helix_chain[h],helix_start[h]+1,helix_chain[h]);
				helix_start[h] += 1;
				}
			}
		found = 0;
		while (!found)
			{
			r = 0;
			while ((r < residue_index) &&
				((helix_residue_no[r] != helix_end[h]) || (helix_chain[helix_no[r]] != helix_chain[h])))
				r++;
			if ((helix_residue_no[r] == helix_end[h]) && (helix_chain[helix_no[r]] == helix_chain[h]))
				found = 1;
			else
				{
				printf("pruning end of helix %d from %d:%c to %d:%c\n",h,helix_end[h],
					helix_chain[h],helix_end[h]-1,helix_chain[h]);
				helix_end[h] -= 1;
				}
			}
		}
	} /* end of prune_extended_helices */

int read_helical_pdb()
	{
	int lastresidue ,atom_index, last_residue_index, serial, resSeq,
		residue_index, i,j, null_heterogen, mapped_residue_aacode;
	char name[5] = "XXXX", resName[4] = "XXX", segID[5] = "XXXX",
		stdRes[4] = "XXX", lastchain,chainID,altLoc,iCode,lastiCode,aacode;

	float x,y,z,occupancy,tempFactor;

	lastresidue = -9999;
	lastchain = '\0';
	lastiCode = '\0';
	is_helical = false;
	residue_index = 0;
	atom_index = 0;
	last_residue_index = 0;

	while (fgets(record_type,7,pdb_file) != NULL)
		{
/*		printf("\"%s\"\n",record_type); */
		fgets(textstring,MAX_LINE_WIDTH,pdb_file);

		/* check MODRES records */

		if (!strcmp(record_type,"MODRES"))
			{
			sscanf(textstring,"%*c%*c%*c%*c%*c%*c%c%c%c%*c%*c%*c%*c%*c%*c%*c%*c%*c%c%c%c",
				&resName[0],&resName[1],&resName[2],
				&stdRes[0],&stdRes[1],&stdRes[2]);

			for (i = 1; i < 3; i++)
				if (isupper(resName[i])) resName[i] = tolower(resName[i]);

			/* check it against the list of known modified residues which
			should be treated as standard residues for the purposes of SOCKET */

			for (i = 0; i < n_heterogens; i++)
				{
				if (setflag[flag_debug])
					printf("identifying MODRES residue: \"%s\" v \"%s\"",
								resName,heterogen3[i]);

				if (strcmp(resName,heterogen3[i]) == 0)
					{
					if (!setflag[flag_q])
						printf("MODRES record specifies %s; will treat as %s\n",
								resName,amino_acid3[map_alpha3_to_amino_acid(resName)]);

					for (i = 1; i < 3; i++)
						if (isupper(stdRes[i])) stdRes[i] = tolower(stdRes[i]);

					if (strcmp(stdRes,amino_acid3[map_alpha3_to_amino_acid(resName)]))
						printf(" - !!! but MODRES record says its a modified %s !!!\n",
							stdRes);

					break;
					}

				if (setflag[flag_debug])
					printf("\n");

				}

			if (i >= n_heterogens)
				{
				printf("MODRES record specifies previously unlisted residue \"%s\" (%s);\n\t- will treat %s as %s\n",
					resName, textstring+23, resName, stdRes);

				/* add this new heterogen-type residue to the list */

				/* - unless the list is full */

				if (n_heterogens == HETEROGENS_MAX)
					{
					printf("Maxmimum number of heterogen-type residues (%d) exceeded\n",
						HETEROGENS_MAX);
					exit(1);
					}

				map_heterogen_no_to_amino_acid[n_heterogens] =
					map_alpha3_to_amino_acid(stdRes);


				/* XXX printf("&(heterogen3[n_heterogens]) = %d\nheterogen3[n_heterogens]=%d\n",&(heterogen3[n_heterogens]),heterogen3[n_heterogens]);
printf("heterogen3[0] = \"%s\"; *(heterogen3[0]) = \"%c\"\n",heterogen3[0],*(heterogen3[0]));
/* *(heterogen3[0]) = 'F'; */
				strcpy(heterogen3[n_heterogens],resName); /* pointers */

				printf("%s added to list; treating as amino acid (%d), %s\n",
					heterogen3[n_heterogens], map_heterogen_no_to_amino_acid[n_heterogens],
					amino_acid3[map_heterogen_no_to_amino_acid[n_heterogens]]);

				n_heterogens++;

				}

			}

		else

		/* only process ATOM and HETATM records */
		if (!(strcmp(record_type,"ATOM  ") && strcmp(record_type,"HETATM")) )
			{
			/* the sscanf line is unwieldy, as the name, resName and
			segID strings are read a character at a time; this is
			because sscanf does not read spaces into strings, and
			these strings will usually contain spaces; it could be
			done with several sscanf/fgets statements instead */
			sscanf(textstring,"%d%*c%c%c%c%c%c%c%c%c%*c%c%4d%c%*c%*c%*c%f%f%f%f%f%*c%*c%*c%*c%*c%*c%c%c%c%c",
			/*                  |   |______|   |___|    |  | iCode      x y z | |                   |_____|
			                serial    name    resName   | resSeq      occupancy  tempFactor          segID
			                                         chainID

					ATOM    394  NE2 HIS A  48A     37.416  27.800  50.108  1.00 48.89      201L 520 */

				&serial,
				&name[0],&name[1],&name[2],&name[3],
				&altLoc,
				&resName[0],&resName[1],&resName[2],
				&chainID,
				&resSeq,
				&iCode,
				&x,&y,&z,
				&occupancy,
				&tempFactor,
				&segID[0],&segID[1],&segID[2],&segID[3]);

			if (setflag[flag_debug])
				printf("serial=%d, name=\"%s\", altLoc='%c', resName=\"%s\", chainID='%c', resSeq=%d, iCode='%c', x=%8.3f, y=%8.3f, z=%8.3f, occupancy=%6.2f, tempFactor=%6.2f, segID=\"%s\"\n",serial,name,altLoc,resName,chainID,resSeq,iCode,x,y,z,occupancy,tempFactor,segID);


			/* 17-3-1
				Watch out for modified amino acid residues like selenocysteine and
				selenomethionine; these	correspond to standard amino acids but are
				invariably specified in HETATM records not ATOM, and the one-
				letter code in the DSSP file is the same as the corresponding
				standard amino acid.

				This kludge looks for any such non-standard amino acids (in array
				heterogen3[]) and simply substitutes the 3-letter code of any it
				finds, with the corresponding standard residue's 3-letter code,
				before proceeding as normal. BOLLOX theres no kludge. XXX

				However, if any HETATM records are encountered which are not
				part of the polypeptide chain (i.e. are not listed in heterogen3[])
				then they should be completely ignored.
			*/

				null_heterogen = true;

				if (strcmp(record_type,"HETATM") == 0)
					{


					/* heterogen3[] elements are stored in Xxx format;
					this kludge allows the strcmp() below to work */
					for (i = 1; i < 3; i++)
						if (isupper(resName[i])) resName[i] = tolower(resName[i]);

					if (setflag[flag_debug])
						printf("HETATM record: resName changed to \"%s\"\n",resName);

					for (i = 0; i < n_heterogens; i++)
						{
						if (setflag[flag_debug])
							printf("identifying HETATM residue: \"%s\" v (%3d)\"%s\"",
								resName,i,heterogen3[i]);

						if (strcmp(resName,heterogen3[i]) == 0)
							{
							if (setflag[flag_debug]) printf(" MATCH\n");
							/* change the residue name (messy) */
							/* NO NEED, as map_alpha3_to_amino_acid(resName) comes XXX
							up with the right answer */
							/*strcpy(resName,amino_acid3[map_heterogen_no_to_amino_acid[i]]);*/
							null_heterogen = false;
							break;
							}

						if (setflag[flag_debug])
							printf("\n");

						}


					if (null_heterogen)
						{
						/* its not a recognized amino-acid heterogen, so
						check that its a recognized solvent; in any case,
						it will be ignored */

							/* note the following check is done only if -u has
							been specified */

						if (setflag[flag_u])	{
							for (i = 0; i < NON_AA_HETEROGENS; i++)
								{
/*printf("i=%d, \"%s\"\n",i,solvent_heterogen3[i]);*/
								if (strcmp(resName,solvent_heterogen3[i]) == 0)
									break;
								}

							if (i == NON_AA_HETEROGENS)
								printf("Unidentified heterogen: %s, %d:%c, iCode='%c'\n",
								resName, resSeq, chainID, iCode);
							}


						continue; /* abort this HETATM record */
						}
					}

				if (chainID == ' ') j = 0;
				else if ((chainID >= 'A') && (chainID <= 'Z')) j = chainID - 'A' + 1;
				else if ((chainID >= 'a') && (chainID <= 'z')) j = chainID - 'a' + 27;
				else j = chainID - '0' + 53;
				if (seqchain[j] == NULL_SEQCHAIN)	{
				seqchain[j] = resSeq;
				printf("chain %c starts at residue %4d, iCode='%c'\n",chainID,resSeq,iCode);
				}

			/* ignore hydrogen atoms unless the -a command-line flag has been used */

			if ( (name[1] == 'H') &&  !setflag[flag_a] ) continue;

			/* is this atom in the same residue as the last one read in,
			or is it the first atom of a new residue? */

			if ((resSeq == lastresidue) && (chainID == lastchain) && (iCode == lastiCode))
				{ /* only bother with residues which are already known to be in
					an alpha-helix */

				if (is_helical) {
					if (atom_index == MAX_ATOMS)
						{
						printf("Maximum no of atoms (%d) exceeded\n",MAX_ATOMS);
						exit(1);
						}
					if (strcmp(name,refatom0type) == 0) refatom0[last_residue_index] = atom_index;
					if (strcmp(name,REFATOM3) == 0) refatom3[last_residue_index] = atom_index;
					aacode = map_alpha3_to_amino_acid(resName);
					for (j = 0; j < 2; j++)
						if (strcmp(name,refatomtype[j][aacode]) == 0)
							{
							refatom1[j][last_residue_index] = atom_index;
							break;
							}
					atom_no[atom_index] = serial;
					strcpy(atom_name[atom_index],name);
					atom_res[atom_index] = last_residue_index;
					coord[atom_index][0] = x;
					coord[atom_index][1] = y;
					coord[atom_index++][2] = z;
					/*printf("%d,%d,%d,%s,%8.3f,%8.3f,%8.3f,%c,%4d\n",atom_index,last_residue_index,atom_no[atom_index-1],atom_name[atom_index-1],x,y,z, chainID,resSeq );*/
					if ((strcmp(atom_name[atom_index-1]," CA ")==0)&&(strcmp(resName,"GLY"))==0)
					{
						atom_no[atom_index]=serial-2;
						strcpy(atom_name[atom_index],"CB");
						atom_res[atom_index] = last_residue_index;
						coord[atom_index][0] = x+1.126;
						coord[atom_index][1] = y+0.872;
						coord[atom_index++][2] = z+0.512;
						/*printf("%s,%s\n",atom_name[atom_index-1],resName);
						printf("%d,%d,%d,%s,%8.3f,%8.3f,%8.3f,%c,%4d\n",atom_index,last_residue_index,atom_no[atom_index-1],atom_name[atom_index-1],x,y,z, chainID,resSeq );*/
					}

					}
				} /* end if ((resSeq == lastresidue) && (chainID == lastchain) && (iCode == lastiCode))*/
			else
				{
				/* its the first atom of a new residue */
				/* check that this residue is in an alpha-helix */


				i = within_helix(resSeq, iCode, chainID, helix_start, helix_start_iCode,
						helix_end, helix_end_iCode, helix_chain, helix_index);

					if (i != -1)
						{
						/* this residue falls between the start and end of a known helix
						(as specified in DSSP in conjunction with any helix-extension) */

						is_helical = true;

						/* check that the new residue name (3-letter) corresponds with the
						1-letter code specified in the DSSP */

						mapped_residue_aacode = map_alpha3_to_amino_acid(resName);

						/* mapped_residue_aacode is the numeric code for an amino acid
						residue as translated from a 3-letter string; see function
						map_alpha3_to_amino_acid() in aminoa2.h; see also aminoa1.h
						and statchar.c */

						if (helix_residue_aacode[residue_index] != mapped_residue_aacode)
								{
								if (helix_residue_aacode[residue_index]
									|| setflag[flag_v] || setflag[flag_debug])
								printf("DSSP file doesnt match PDB file: residue (%d) %d:%c iCode='%c'\n\tDSSP %d (\"%s\") v PDB %d (\"%s\"; original \"%s\")\n",
								residue_index, resSeq, chainID, iCode,
								helix_residue_aacode[residue_index],
								amino_acid3[helix_residue_aacode[residue_index]],
								mapped_residue_aacode,
								amino_acid3[mapped_residue_aacode],
								resName);

/* the amino acid previously read in from the DSSP file (represented by helix_residue_aacode[residue_index]) is not the same as the
one just read in (resName) from the PDB file. This indicates that the DSSP and PDB files don't match. Another possibility is that
the amino acid in the DSSP is an 'X' residue, which is stored in the PDB file as a series of HETATM records instead of ATOM. In this
case, the HETATM records will have been ignored, so (assuming there are x consecutive such HETATM residues) the current residue
read from PDB will correspond to residue_index + x . Therefore incrementing residue_index by x will get things back in synch. This is
acheived by adding one to the value of residue_index until helix_residue_aacode[residue_index] is not 0 (= 'X'). */

/* XXX !!!! SHOULD PROBABLY BOMB OUT HERE !!!! */

/*								while (helix_residue_aacode[residue_index] == 0)
									{
									printf("WARNING: residue %d (%d %c) in the DSSP file is 'X' and does not match the PDB file; a possible reason is that this residue is stored as HETATM rather than ATOM records; please check the PDB file\n",
									residue_index++,resSeq,helix_chain[helix_no[residue_index]]);
									}
XXX delete this crap*/

								/* a mismatch between an XXX (code == 0) and a 'special' heterogen such
								as Cse or Mse is allowed */

								if (!(helix_residue_aacode[residue_index] || null_heterogen))
									{
									if (helix_residue_aacode[residue_index]
										|| setflag[flag_v] || setflag[flag_debug])
									printf("\t- allowing (DSSP files list %s as 'X')\n",resName);
									}

								/* now residue_index is equivalent to a non-X residue; if there is
								still a mismatch, bomb out */

								else if (helix_residue_aacode[residue_index] !=
									map_alpha3_to_amino_acid(resName) )
									{
									printf("DSSP file doesnt match PDB file:\n");
/*printf("%d\n",helix_residue_aacode[residue_index]);*/
									printf("\tDSSP: residue %d) %d iCode='%c' chain %c is %s (%s)\n",
										residue_index,resSeq,
										helix_residue_iCode[residue_index],
										helix_chain[helix_no[residue_index]],
										amino_acid1[helix_residue_aacode[residue_index]],
										amino_acid3[helix_residue_aacode[residue_index]]);

									printf("\tPDB:  residue %d) %d iCode='%c' chain %c is %s (%s)\n",
										residue_index,resSeq,iCode,
										helix_chain[helix_no[residue_index]],
										amino_acid1[map_alpha3_to_amino_acid(resName)],resName);
									exit(1);
									}
								}

							strcpy(helix_residue_name[residue_index++],resName);

							/*helix_residue_no[residue_index] = resSeq;
							helix_residue_chain[residue_index++] = chainID; these already set by
							read_helical_dssp*/

							if (atom_index == MAX_ATOMS)
								{
								printf("Maximum no of atoms (%d) exceeded\n",MAX_ATOMS);
								exit(1);
								}

							atom_no[atom_index] = serial;
							strcpy(atom_name[atom_index],name);

							for (last_residue_index = 0;
								last_residue_index < residue_index; last_residue_index++)
								if ((helix_residue_no[last_residue_index] == resSeq)
									&& (helix_residue_iCode[last_residue_index] == iCode)
									&& (helix_chain[helix_no[last_residue_index]] == chainID))
									{
									atom_res[atom_index] = last_residue_index;
									break;
									}
								else if (last_residue_index == residue_index -1)
									{
									printf("oops- couldnt find this residue (resSeq=%d,chainID='%c') in the list read from the DSSP file\n",
									resSeq,chainID);
									exit(1);
									}
							if (strcmp(name,refatom0type) == 0) refatom0[last_residue_index] = atom_index;
							if (strcmp(name,REFATOM3) == 0) refatom3[last_residue_index] = atom_index;
							aacode = map_alpha3_to_amino_acid(resName);
							for (j = 0; j < 2; j++)
								if (strcmp(name,refatomtype[j][aacode]) == 0)
									{
									refatom1[j][last_residue_index] = atom_index;
									break;
									}
							coord[atom_index][0] = x;
							coord[atom_index][1] = y;
							coord[atom_index++][2] = z;
							/*break;*/
						} /* end if (i == -1)*/
					else is_helical = false;
					lastresidue = resSeq;
					lastchain = chainID;
					lastiCode = iCode;

				} /* end else of  ((resSeq == lastresidue) && (chainID == lastchain) && (iCode == lastiCode)) */
			}
		else if (strcmp(record_type,"ENDMDL") == 0)
			{
			printf("ENDMDL card found: implies this is PDB file contains multiple NMR models; all but the first will be ignored\n\n");
			break;
			}

		}

	if (setflag[flag_v])
		{
		printf("These are the alpha-helical residues:\n\n");
		for (i = 0; i < residue_index; i++)
			printf("\t%d) %s %d iCode='%c' %c (helix %d)\n",i,helix_residue_name[i], helix_residue_no[i],
			helix_residue_iCode[i], helix_chain[helix_no[i]], helix_no[i]);

		printf("\nThese are the atoms in the above residues:\n\n");
		for (i = 0; i < atom_index; i++)
			printf("\t%d) %s %d %8.3f %8.3f %8.3f  residue %d (%s %d iCode='%c' %c)\n",i, atom_name[i], atom_no[i], coord[i][0], coord[i][1], coord[i][2], atom_res[i],helix_residue_name[atom_res[i]], helix_residue_no[atom_res[i]], helix_residue_iCode[atom_res[i]], helix_chain[helix_no[atom_res[i]]]);

		printf("\nThese are the reference atoms for each residue:\n\n");
		for (i = 0; i < residue_index; i++)
			printf("\tresidue %d) CA: atom %d; CB: atom %d; end1: atom %d; end2: atom %d\n",
			i,refatom0[i],refatom3[i],refatom1[0][i],refatom1[1][i]);
		}

	return atom_index;
	/* end of function read_helical_pdb */
	}
