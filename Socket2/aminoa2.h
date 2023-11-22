/*

					SOCKET
					 v3.02

					aminoa2.h

					16-10-01
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

/*					aminoa2.h
					---------
*/

/* JW 4-5-98

	Function to initialize integer array map_alpha_to_amino_acid[],
	which converts an ascii character (from 'A' to 'Z') to the
	correct value as specified by enum amino_acid (see aminoa2.h).
	The reason that everything is not done directly by one array (of
	26 elements, corresponding to 'A' to 'Z' , with null values for
	'J', 'O', etc.) is that the declaration of the type enum amino_acid
	should ideally be flexible - for instance you might want to group
	hydrophobic side chains together, to make things easier to process.
	This way, the amino acids can be declared in any order in the
	enum amino_acid {...} line.

	The function below should be called near the beginning of your code.
	Or before you want to actually use the enum amino_acid type anyway.
*/

void aa_map()
	{
	int i;
	char st[2];
	enum amino_acid aa;
	/* remember char is really an integer type */
	for (i = 0; i <26; i++)
		{
		map_alpha_to_amino_acid[i] = 0;
		for (aa = X; aa < AMINO_ACIDS; aa++)
			{
			if ((i + 'A') == amino_acid1[aa][0])
				{
				map_alpha_to_amino_acid[i] = aa;
				break;
				}
			}
		}
	if (setflag[flag_debug]) for (i = 0; i <26; i++)
		{
		st[0] = i + 'A'; st[1] = '\0';
		printf("char \"%s\" maps to %d, ie \"%s\" , \"%s\"\n", st,map_alpha_to_amino_acid[i],amino_acid1[map_alpha_to_amino_acid[i]],amino_acid3[map_alpha_to_amino_acid[i]]);
		}
	}
	
/*	Function to convert a single 3-letter string to type enum
	amino_acid. If the 3-letter string is unrecognized, amino
	acid 'X' is returned
*/

int map_alpha3_to_amino_acid(char aa_string[])
	{
	char ch3_string[4] = "   ";
	int i,j;
	j = 0;
	/* ignore any leading spaces etc */
	for (i = 0; i < strlen(aa_string); i++)
		if (isalpha(aa_string[i])) ch3_string[j++] = aa_string[i];
	if (islower(ch3_string[0])) ch3_string[0] = toupper(ch3_string[0]);
	for (i = 1; i < 3; i++)
		if (isupper(ch3_string[i])) ch3_string[i] = tolower(ch3_string[i]);

	/* now check against the array of 3-letter amino acid strings */
	j = 0;
	for (i = 0; i < AMINO_ACIDS; i++)
		if (strcmp(amino_acid3[i],ch3_string) == 0)
			{ j = i; break;}

	/* if no match, do the same for the list of recognized 3-letter heterogen strings
		- these will still map to one of the normal amino acids */
	if (j == 0)	for (i = 0; i < n_heterogens; i++)
				if (strcmp(heterogen3[i],ch3_string) == 0)
					{ j = map_heterogen_no_to_amino_acid[i]; break;}
	return j;
			
	}

/*
	Function to convert a string of (single-letter) amino acid codes
	to an array of type enum amino_acid (see aminoa2.h)
*/

int aa_interpret(char sequence_string[], int sequence[], int code_1_or_3)
	{
	int i,alpha,length,residues;
	char st[1];
	length = strlen(sequence_string);
	if (code_1_or_3 != 3) {
		residues = -1;
		for (i = 0; i < length; i++)
			{ residues++;
/* 10 is the ascii value of a carriage return to mark the end of a line */
			if ((sequence_string[i] == 10) || (sequence_string[i] == '*') || (sequence_string[i] == '/')) break;
			if (islower(sequence_string[i])) { sequence_string[i] = toupper(sequence_string[i]);}
			alpha = sequence_string[i] - 'A';
			if ((alpha < 0) || (alpha > 25))
				{ st[0] = sequence_string[i]; 
				printf("sequence has non alphabetic character (\"%s\", chr %d) at position %d\n\n",st,sequence_string[i],i+1);
				exit(1);
				}
			sequence[i] = map_alpha_to_amino_acid[alpha];
			} 
		}
	return residues;
	}
