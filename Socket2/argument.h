#include <string.h>
#include <stdlib.h>
/*

					SOCKET
					 v3.02

					argument.h

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

/*					argument.h
					----------
*/

/*JW 1-5-98*/
/* This function is used to set 'flags' (which are either on or off) and
	'parameters' (which have a value) used by the main function of a
	C program. The user executing the program can specify which flags
	are set, and the parameter values, like this:
		> cprogname -a -blah -c 20 -val thingy
	- the above would set the flags called 'a' and 'blah' to true, and
	would give the parameter 'c' the value '20', and 'val' the value
	'thingy' .
	How these flags and parameters are referenced from within your program
	is partly a matter of choice, but a way of doing it is described
	below.
	Currently, all parameter *values* are strings (arrays of characters).
	If you want to use them in a different context, then you must convert
	them yourself in your C program using the atoi(), atof() or atol()
	functions for example.
	Later versions of this routine may allow specification of the desired
	type of each parameter, to do the conversion automatically.

	To use this function, you should include the following in the main
	program	file:
		#define MAXQUALNAMELENGTH <m>
		#define FLAGS <f>
		#define PARS <p>
	-where m is the maximum character length of the qualifiers that are used;
	f is the number of different flags recognized by the program;
	p is the number of different parameters recognized by the program.

	In addition, you need to specify your flag an parameter names as
	enumerated types in the main program file, like this:
		enum flags {flag_e,flag_v};
		enum pars {par_f,par_l,par_p,par_s};
	- for example. This just means that 'flag_e' is equivalent to 0,
	'flag_v' to 1, 'par_f' to 0, 'par_l' to 2 etc .

		enum boolean {false,true};


*/

int set_flagsNpars(int argc, char *argv[])
	{
	/*enum boolean {false,true};GLOBAL*/
	enum boolean found;

	char qualifier[MAXQUALNAMELENGTH];
	int i,j;
	i = 1;
	if (setflag[flag_debug]) printf("starting loop\n");
	while (i < argc)
		{ if (setflag[flag_debug]) printf("looping\n");
		found = false;
		for (j = 0; j < FLAGS; j++)
			{
			if (found) continue;
			*(qualifier) = '-';
			*(qualifier+1) = '\0';
/*			printf("\"%s\"\n",qualifier);*/
			strcat(qualifier,flagname[j]);
/*			printf("\"%s\"\n",qualifier);*/
			if (strcmp(qualifier,argv[i]) == 0)
				{
				setflag[j] = true; 
				found = true;
				}
			}
		for (j = 0; j < PARS; j++)
			{
			if (found) continue;
			*(qualifier) = '-';
			*(qualifier+1) = '\0';
/*			printf("\"%s\"\n",qualifier);*/
			strcat(qualifier,parname[j]);
/*			printf("\"%s\"\n",qualifier);*/
			if ((strcmp(qualifier,argv[i]) == 0) && (i < argc - 1))
				{
/*				printf("hey: \"%s\"\n",argv[++i]);*/
				par[j] = argv[++i]; 

				found = true;
				}				
			}
		if (setflag[flag_debug]) printf("i = %d\n",i);
		i++;
/*		printf("finishing this loop (i now set to:%d)\n",i); */
		}


		return 0;
}
