/*

					SOCKET
					 v3.02

					daisies.c

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

/*					daisies.c
					---------
*/


#include "socket.h"

/* this subroutine is intended to be called like this:

order_of_knob_k = check_daisy_chain(k,0,daisy_chain_2Darray[n],direction);

-where daisy_chain_2Darray is a two dimensional array, each
element of which corresponds to a single daisy chain and is an array
storing the serial numbers of all the knobs which are part of that
daisy chain. So, n is the daisy chain serial number;
daisy_chain_2Darray[n] is a 1D array storing the knob IDs of the
members of the nth daisy chain in the list.
When the above call is made, daisy_chain_2Darray[n] should be
blank (initialized to -1 values).

The daisy_chain_2Darray[n] of the call is therefore equivalent
to the daisy_chain array in the check_daisy_chain function

direction is either 1 or -1, determining in what order to check
the 2 sides of a hole (see note 8-12-99 below)

The check_daisy_chain function updates the daisy_chain array
by calling itself recursively.

When called, the function does the following things:

	it adds knob k to the next free element of the
		daisy_chain array; the current size of the
		array is held in the variable order, which is
		zero when the function is first called externally

	it increments the size of the daisy_chain array by one

	it checks the two residues which form the two sides of the
		hole into which knob k fits.
		(All the knobs are held in the global array knob;
		the global variable knob_index is how many
		knobs there are; hole residues are stored in
		global array hole[][])
		These two residues are the second and third of the
		four residues in the hole array, ie hole[k][1] and
		hole[k][2].

		If either of these two hole residues are themselves
		knobs, which *fit into holes in a helix OTHER than
		the one of which knob k is part*, then this indicates
		that knobs k and the new knob (call it k2) form a link
		in a 'chain' of knobs.

		That is, knob k fits into a hole, the sides of which
		are k2 and h2 (h2 is a side of k's hole which is not a
		knob, or at least not in this knob-chain), and k2 fits
		into a hole whose sides are k3 and h3; and neither
		k3 nor h3 are knob k. This sort of arrangement usually
		only occurs in cyclic (closed) chains of knobs.

		(If either k3 or h3 were in fact k, then this would
		simply represent a pairwise complementary interaction,
		ie order 2. This function is not concerned with these.)

		If k2 turns out to be first knob in the daisy_chain
		array (this can only happen in the first *recursive* call
		of the function or later), then the circle is complete,
		and the function terminates by returning the order. Eg
		the chain is: k -> k2 -> k3 -> k4 -> k -> k2 -> k3 -> k4...
		ie the order is 4 and the elements of the array are:
		daisy_chain[0] = k; daisy_chain[1] = k2;
		daisy_chain[2] = k3; daisy_chain[3] = k4. 

		If k2 is in fact any other knob in the daisy_chain
		array, then the chain is 'partially cyclic' ie it goes
		k -> k2 -> k3 -> k4 -> k2 -> k3 -> k4 -> k2 -> k3 -> k4...
		- so k1 leads into the chain but is not actually part
		of the cycle. If this is the case, the function terminates
		with a negative (null) value - there is no point in dealing
		with the daisy chain here, because the cyclic portion will
		be dealt with elsewhere (when the function is called
		externally to check knob k2).

		If k2 is not already in the daisy_chain array, then the
		function calles itself recursively, to check k2 in the
		same way that k has been checked:
		check_daisy_chain(k2,order,daisy_chain)

	if neither of the two residues forming the sides of k's hole are
		themselves knobs which fit into a hole in a helix other
		than that of which k is part, then the function terminates
		with a negative (null) value.


Therefore, the recursion will end when either a negative (null) value, or
the order of a complete, cyclic daisy chain, is returned.

*/

int check_daisy_chain(int thisknob, int order, int daisy_chain[], int direction)
	{

	int g,h,i,j,result;

	result = -1;

	daisy_chain[order++] = thisknob;

	/* f = (x+1+ (d-1)*(2x-1)/2)*/

	for (g = 0; g < 2; g++)
		{
		h = g + 1 + ((direction - 1) * (2*g - 1)/2);
		for (i = 0; i < knob_index; i++)
			{
			if (knob[i] == hole[thisknob][h])
				{

				/* check that knob[i]'s hole is not in the
				helix of previousknob */
				if (helix_no[hole[i][0]] != 
					helix_no[knob[thisknob]])
					{

					j = 0;
					while ((j < order) && (hole[i][1] != knob[daisy_chain[j]]) 
						&& (hole[i][2] != knob[daisy_chain[j]])) j++;

					if ((hole[i][1] == knob[daisy_chain[j]]) || (hole[i][2] == knob[daisy_chain[j]]))
						{

						daisy_chain[order++] = i;
						if (j) result = -2;
						else result = order;
						}
					else result = check_daisy_chain(i,order,daisy_chain,direction);
					}
				}
			if (result != -1) break;
			}
		if (result != -1) break;
		}
	return result;
	}


/*
eg 1aik at cutoff 8.2

residue	1	Q548:N		helix_no[1] = 0
residue	4	Q551:N		helix_no[4] = 0
residue	5	Q552:N		helix_no[5] = 0
residue	8	Q555:N		helix_no[8] = 0
residue	65	I548:A		helix_no[65] = 2
residue	68	Q551:A		helix_no[68] = 2
residue	69	Q552:A		helix_no[69] = 2
residue	72	L555:A		helix_no[72] = 2
residue	129	I548:D		helix_no[129] = 4
residue	132	Q551:D		helix_no[132] = 4
residue	133	Q552:D		helix_no[133] = 4
residue	136	L555:D		helix_no[136] = 4

knob[17] = 5	hole[17][0] = 129	hole[17][1] = 132	hole[17][2] = 133	hole[17][3] = 136
knob[14] = 69	hole[14][0] =   1	hole[14][1] =   4	hole[14][2] =   5	hole[14][3] =   8
knob[52] = 133	hole[52][0] =  65	hole[52][1] =  68	hole[52][2] =  69	hole[52][3] =  72

suppose that residue 132 is also a knob (99) in a hole formed by residues of helix 0; ie knob[99] = 132;

order = check_daisy_chain(17,17,0,daisy_chain[daisy_chains]);

(level 0)
firstknob = 17;
thisknob = 17;
result = -1;
h = 1;
i = 0;
i = 1;
i = ..
i = 99
knob[99] == hole[17][1]
but helix_no[hole[99][0]] == helix_no[knob[17]] == 0
i = ..
h = 2;
i = ..
i = 52
knob[52] == hole[17][2] == 133
helix_no[hole[52][0]] == 2 != helix_no[knob[17]]
hole[52][1] == 68 != knob[17]
hole[52][2] == 69 != knob[17]
result = check_daisy_chain(17,52,1);

	(level 1)
	firstknob = 17;
	thisknob = 52;
	order = 1;
	result = -1;
	h = 1;
	i = 0;
	i = 1;
	i = ..
	h = 2;
	i = ..
	i = 14;
	knob[14] == hole[52][2] == 69
	helix_no[hole[14][0]] == 0 != helix_no[knob[52]]
	hole[14][1] == 4 != knob[17]
	hole[14][2] == knob[17] == 5
	result = 1

(level 0)
result = 1
*/

/*

8-12-99

Problem with finding all daisy chains in 'complex' coiled coils
(ie central coiled coils with helices packed on the outside forming
outer coiled coils). Its possible to omit some of the core daisy
chains, resulting in incorrect register assignment, because the
core side chains appear to be taking part *only* in peripheral
KiH interactions, so don't get labelled as 'a' or 'd'.

E.g. structure 1ebo

helices 1, 5 and 9 (chains A,B and C respectively) form a 3-stranded
parallel coiled coil (coiled coil iii)
around the outside of this bundle are helices 3, 7 and 11 (also
chains A, B and C), giving these coiled coils:

helices 1, 9, 11 (antiparallel) (coiled coil iv)
helices 1, 3, 5 (antiparallel) (coiled coil ii)
helices 5, 7, 9 (antiparallel) (coiled coil v)

Consider this KiH layer:

Leu 68:A (Residue #41, Knob #13)
Leu 68:B (Residue #124, Knob #48)
Leu 68:C (Residue #213, Knob #25)

41 (LEU 68:A, helix 1)
hole:
	120 (THR 64:B, helix 5)
        123 (ALA 67:B, helix 5)
        124 (LEU 68:B, helix 5)
	127 (PHE 71:B, helix 5) 

124 (LEU 68:B, helix 5)
hole:
	209 (THR 64:C, helix 9)
        212 (ALA 67:C, helix 9)
        213 (LEU 68:C, helix 9)
	216 (PHE 71:C, helix 9)

213 (LEU 68:C, helix 9)
hole:
         37 (THR 64:A, helix 1)
	 40 (ALA 67:A, helix 1)
	 41 (LEU 68:A, helix 1)
	 44 (PHE 71:A, helix 1)


Each of these core residues (Leu 68) is also involved in
an outer, antiparallel coiled coil:

Ala 67:A (Residue #40, Knob #30)
Ile 122:C (Residue #249, Knob #61)

40 (ALA 67:A, helix 1)
hole:
	245 (ILE 118:C, helix 11)
	248 (LYS 121:C, helix 11)
	249 (ILE 122:C, helix 11)
	252 (ILE 125:C, helix 11) 

249 (ILE 122:C, helix 11)
hole:
	210 (THR 65:C, helix 9)
	213 (LEU 68:C, helix 9)
        214 (GLN 69:C, helix 9)
	217 (LEU 72:C, helix 9)


There is also:
Lys 121:C (Residue #248, Knob #32)
248 (LYS 121:C, helix 11)
hole:
	36 (GLU 63:A, helix 1)
	39 (GLN 66:A, helix 1)
	40 (ALA 67:A, helix 1)
	43 (LEU 70:A, helix 1)


So, what happens when checking knob 25 (Leu 68:C) ?

order = check_daisy_chains(25,0,daisy_chain[daisy_chains]);

(level 0)

thisknob = 25;
order = 0;
result = -1;

daisy_chain[0] = 25;
order = 1;

h = 1;
i = 0;
i = ..
i = 30

knob[30] == hole[25][1] == 40
helix_no[hole[30][0]] == 11
helix_no[knob[25]] == helix_no[213] == 9

j = 0

hole[30][1] == 248
hole[30][2] == 249
knob[daisy_chain[0]] == knob[25] == 213

j = 1;
j >= order;

knob[daisy_chain[1]] == knob[-1] != 248 ; != 249

result = check_daisy_chain(30,1,daisy_chain);

	(level 1)

	thisknob = 30;
	order = 1;
	result = -1;

	daisy_chain[1] = 30;
	order = 2;

	h = 1;
	i = 0;
	i = ..
	i = 32

	knob[32] == hole[30][1] == 248
	helix_no[hole[32][0]] == 1
	helix_no[knob[30]] == helix_no[40] == 1

	i = ..

	h = 2
	i = 0;
	i = ..
	i = 61

	knob[61] == hole[30][2] == 249
	helix_no[hole[61][0]] == 11
	helix_no[knob[30]] == helix_no[40] == 1

	j = 0

	hole[61][1] == 213
	hole[61][2] == 214
	knob[daisy_chain[0]] == knob[25] == 213

	hole[61][2] ==  knob[daisy_chain[0]] == 213

	daisy_chain[2] = 61;

	order = 3;
	j == 0;
	result = 3;

	return 3;

(level 0)
result = 3;

return 3;


- gives order == 3.

Therefore, knob 25 (residue 213) becomes part of the outer, 
antiparallel daisy chain only. This will also happen to its
counterparts on chains A and B, ie knobs 13 and 48.

The solution would be to force the routine to begin by
checking the *second* side of the hole into which knob 25
fits (instead of the third) - ie the third hole residue,
h = 2, not h = 1. If h = 2, then the side residue would be
Leu 68:A, and so long as this in turn checked h = 2 of
68:A's hole first, the central core daisy chain would be
found.

A solution could therefore be to first call check_daisy_chain
(to determine the order of a knob) with h = 1 checked first,
and then do a separate call with h = 2 checked first. However,
what happens with the recursive calls of check_daisy_chain?
The function cannot call itself twice (one for h = 1 first,
one for h = 2 first), because each sequence of recursive calls
is adding to a single daisy chain; the second call would
overwrite the results of the first.

Therefore, it would be problematic to examine every possibility
- eg with a four-stranded coiled coil, the recursion would be 3
deep, so there would be 8 paths to follow to be sure of finding
all daisy chains.

The best solution is to call
check_daisy_chain(i,0,daisy_chain[daisy_chains])
*twice* from the subroutine check_complementarity() (currently
the source is in socket2.06.c) ; the first would check h = 1
before h = 2; the second h = 2 before h = 1. This will require
another variable, which will result in ascending h when 1, and
descending h when -1 . All recursive self-calls of
check_daisy_chain will have the same direction as the calling
function.
*/
