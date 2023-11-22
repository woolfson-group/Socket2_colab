/*

					SOCKET
					 v3.02

					statchar.c

					24-10-01
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

/*					statchar.c
					----------
*/


#include "preproc.h"
#include "global.h"
#include <stdlib.h>
#include <string.h>
/* statchar: initializes the values of static character arrays specifying
	command-line flags and parameters, and 1- and 3-letter codes for
	amino acid residues
*/


void statchar()
	{
	int i;
	/* the names of the available flags, as specified by the user on the command-
	line when preceded directly by '-' */

	static char *local_flagname[FLAGS] = {"debug","a","b","i","l","p","q","u","v"};


	/* the names of the available parameters, as specified by the user on the
	command-line when preceded directly by '-' */

	static char *local_parname[PARS] = {"c","d","e","f","j","k","o","r","s","t","w"};



	/* The amino acid codes in the enum type amino_acid can be specified in
		any order (you might want to change it so that they are grouped
		by property eg hydrophobicity). Type enum amino_acid is mapped
		from the corresponding ascii characters by the int array
		map_alpha_to_amino_acid[]. If you do change the order, make sure
		that amino_acid1[] and amino_acid3[] are in the same order...

		MAKE SURE THAT 'X' THE FIRST CODE THOUGH.

		The other thing is, make sure you have a line calling the function
		'aa_map()' in your code.

	*/

	static char *local_amino_acid1[AMINO_ACIDS] = {"X","A","C","D","E","F","G","H","I","K","L",
				"M","N","P","Q","R","S","T","V","W","Y"};

	static char *local_amino_acid3[AMINO_ACIDS] = {"XXX", "Ala", "Cys", "Asp", "Glu", "Phe",
				"Gly", "His", "Ile", "Lys", "Leu", 
				"Met", "Asn", "Pro", "Gln", "Arg", 
				"Ser", "Thr", "Val", "Trp", "Tyr"};

	/* 17-3-1
	Regard these heterogen3[] residue types as equivalent to a normal amino acid
	residue, but recorded as a HETATM record, not an ATOM
	*/

	static char *local_heterogen3[HETEROGENS_MAX] = {

		/*   0 */	"Cse", /* modified CYS selenocysteine */
		/*   1 */	"Mse", /* modified MET selenomethionine */
		/*   2 */	"Mis", /* modified SER monoisopropylphosphorylserine */
		/*   3 */	"Tys", /* modified TYR sulphonated Tyrosine */
		/*   4 */	"Glz", /* modified GLY amino-acetaldehyde */
		/*   5 */	"Aba", /* modified ASN alpha-aminobutyric acid	*/
		/*   6 */	"Agm", /* modified ARG post-translational modification	*/
		/*   7 */	"Aib", /* modified ALA alpha-aminoisobutyric acid	*/
		/*   8 */	"Ala", /* modified ALA amidated terminal carboxy group	*/
		/*   9 */	"Ar2", /* modified ARG arginyl-benzothiazole-6-carboxylic acid	*/
		/*  10 */	"Asn", /* modified ASN glycosylation site	*/
		/*  11 */	"Cas", /* modified CYS s-(dimethylarsenic)cysteine	*/
		/*  12 */	"Ccs", /* modified CYS s carboxymethylated	*/
		/*  13 */	"Cea", /* modified CYS cysteine sulfenic acid	*/
		/*  14 */	"Cgu", /* modified GLU gamma-carboxy-glutamic acid	*/
		/*  15 */	"Clb", /* modified SER covalently bound inhibitor	*/
		/*  16 */	"Cld", /* modified SER covalently bound inhibitor	*/
		/*  17 */	"Cme", /* modified CYS s-(thioethylhydroxy)cystine	*/
		/*  18 */	"Csb", /* modified CYS lead bound to sg	*/
		/*  19 */	"Csd", /* modified CYS modified cysteine	*/
		/*  20 */	"Cso", /* modified CYS s-hydroxycystine	*/
		/*  21 */	"Css", /* modified CYS thiocystine	*/
		/*  22 */	"Csw", /* modified CYS double oxidized cysteine	*/
		/*  23 */	"Cxm", /* modified MET carboxylated methionine	*/
		/*  24 */	"Cyg", /* modified CYS glutamyl-s-cysteine	*/
		/*  25 */	"Cys", /* modified CYS disulfide bond to glutathione	*/
		/*  26 */	"Dar", /* modified ARG d-arginine	*/
		/*  27 */	"Dgl", /* modified GLU d-glutamic acid	*/
		/*  28 */	"Dil", /* modified ILE d-isoleucine	*/
		/*  29 */	"Div", /* modified VAL d-isovaline	*/
		/*  30 */	"Doh", /* modified ASP hydroxylation of cb atom	*/
		/*  31 */	"Dpn", /* modified PHE d-phenylalanine	*/
		/*  32 */	"Ehp", /* modified PHE 3-hydroxyphenylalanine	*/
		/*  33 */	"Fgl", /* modified GLY post-translational oxidation	*/
		/*  34 */	"Ftr", /* modified TRP fluorotryptophane	*/
		/*  35 */	"Gl3", /* modified GLY post-translational modification	*/
		/*  36 */	"Glh", /* modified GLU inhibited by dccd	*/
		/*  37 */	"Hip", /* modified HIS phosphorylated histidine	*/
		/*  38 */	"Htr", /* modified TRP beta-hydroxytryptophane	*/
		/*  39 */	"Hyp", /* modified PRO 4-hydroxyproline	*/
		/*  40 */	"Iil", /* modified ILE allo-isoleucine	*/
		/*  41 */	"Kcx", /* modified LYS carbamylated lysine	*/
		/*  42 */	"Llp", /* modified LYS lysine-pyridoxal-5'-phosphate	*/
		/*  43 */	"M3l", /* modified LYS methylated lysine	*/
		/*  44 */	"Mgn", /* modified GLN 2-methyl-glutamine	*/
		/*  45 */	"Mho", /* modified MET s-oxymethionine	*/
		/*  46 */	"Mhs", /* modified HIS post-translational modification	*/
		/*  47 */	"Mly", /* modified LYS methylated lysine	*/
		/*  48 */	"Mse", /* modified MET selenomethionine	*/
		/*  49 */	"Mty", /* modified TYR meta-tyrosine	*/
		/*  50 */	"Ncb", /* modified ALA chemical modification	*/
		/*  51 */	"Nep", /* modified HIS n1-phosphonohistidine	*/
		/*  52 */	"Nle", /* modified LEU norleucine	*/
		/*  53 */	"Oas", /* modified SER o-acetylserine	*/
		/*  54 */	"Ocs", /* modified CYS cysteinesulfonic acid	*/
		/*  55 */	"Ocy", /* modified CYS hydroxyethylcysteine	*/
		/*  56 */	"Phe", /* modified PHE d-phenylalanine	*/
		/*  57 */	"Phl", /* modified PHE l-phenylalaninol	*/
		/*  58 */	"Prs", /* modified PRO thioproline	*/
		/*  59 */	"Sbd", /* modified SER covalently bound inhibitor	*/
		/*  60 */	"Sbl", /* modified SER covalently bound inhibitor	*/
		/*  61 */	"Sch", /* modified CYS thio methylated cysteine	*/
		/*  62 */	"Seb", /* modified SER o-benzylsulfonyl-serine	*/
		/*  63 */	"Sep", /* modified SER phosphoserine	*/
		/*  64 */	"Smc", /* modified CYS post-translational modification	*/
		/*  65 */	"Snc", /* modified CYS s-nitroso-cysteine	*/
		/*  66 */	"Soc", /* modified CYS dioxyselenocysteine	*/
		/*  67 */	"Sty", /* modified TYR tyrosine-o-sulphonic acid	*/
		/*  68 */	"Sva", /* modified SER serine vanadate	*/
		/*  69 */	"Trf", /* modified TRP ne1-formal-tryptophan	*/
		/*  70 */	"Trn", /* modified TRP aza-tryptophan	*/
		/*  71 */	"Tyi", /* modified TYR 3,5-diiodotryrosine	*/
		/*  72 */	"Tys", /* modified TYR sulfonated tyrosine	*/
		/*  73 */	"Yof", /* modified TYR 3-fluorotyrosine	*/
		/*  74 */	"XYZ", /* blank - can be filled if any	*/
		/*  75 */	"XYZ", /* unexpected MODRES records are	*/
		/*  76 */	"   ", /* encountered; up to 10 unknown	*/
		/*  77 */	"   ", /* types of modified residue can	*/
		/*  78 */	"   ", /* be handled, and added to those*/
		/*  79 */	"   ", /* listed above; the number of	*/
		/*  80 */	"   ", /* types of modified residue is	*/
		/*  81 */	"   ", /* therefore dynamic, and held in*/
		/*  82 */	"   ", /* variable:	  heterogen_modres	*/
		/*  83 */	"   " /* 	*/

	};

	/* the following array specifies to which normal amino acid each heterogen
	is equivalent; e.g. CSE is selenocysteine, equivalent to amino acid 2 (i.e.
	cysteine); note the GLZ is recorded as 'X' by DSSP, unlike CSE,MSE,MIS,TYS
	*/

	int local_map_heterogen_no_to_amino_acid[HETEROGENS_MAX] = {
		/*   0 */	2,	/* modified CYS selenocysteine */
		/*   1 */	11,	/* modified MET selenomethionine */
		/*   2 */	16,	/* modified SER monoisopropylphosphorylserine */
		/*   3 */	20,	/* modified TYR sulphonated Tyrosine */
		/*   4 */	6,	/* modified GLY amino-acetaldehyde */
		/*   5 */	12,	/* modified ASN alpha-aminobutyric acid	*/
		/*   6 */	15,	/* modified ARG post-translational modification	*/
		/*   7 */	1,	/* modified ALA alpha-aminoisobutyric acid	*/
		/*   8 */	1,	/* modified ALA amidated terminal carboxy group	*/
		/*   9 */	15,	/* modified ARG arginyl-benzothiazole-6-carboxylic acid	*/
		/*  10 */	12,	/* modified ASN glycosylation site	*/
		/*  11 */	2,	/* modified CYS s-(dimethylarsenic)cysteine	*/
		/*  12 */	2,	/* modified CYS s carboxymethylated	*/
		/*  13 */	2,	/* modified CYS cysteine sulfenic acid	*/
		/*  14 */	4,	/* modified GLU gamma-carboxy-glutamic acid	*/
		/*  15 */	16,	/* modified SER covalently bound inhibitor	*/
		/*  16 */	16,	/* modified SER covalently bound inhibitor	*/
		/*  17 */	2,	/* modified CYS s-(thioethylhydroxy)cystine	*/
		/*  18 */	2,	/* modified CYS lead bound to sg	*/
		/*  19 */	2,	/* modified CYS modified cysteine	*/
		/*  20 */	2,	/* modified CYS s-hydroxycystine	*/
		/*  21 */	2,	/* modified CYS thiocystine	*/
		/*  22 */	2,	/* modified CYS double oxidized cysteine	*/
		/*  23 */	11,	/* modified MET carboxylated methionine	*/
		/*  24 */	2,	/* modified CYS glutamyl-s-cysteine	*/
		/*  25 */	2,	/* modified CYS disulfide bond to glutathione	*/
		/*  26 */	15,	/* modified ARG d-arginine	*/
		/*  27 */	4,	/* modified GLU d-glutamic acid	*/
		/*  28 */	8,	/* modified ILE d-isoleucine	*/
		/*  29 */	18,	/* modified VAL d-isovaline	*/
		/*  30 */	3,	/* modified ASP hydroxylation of cb atom	*/
		/*  31 */	5,	/* modified PHE d-phenylalanine	*/
		/*  32 */	5,	/* modified PHE 3-hydroxyphenylalanine	*/
		/*  33 */	6,	/* modified GLY post-translational oxidation	*/
		/*  34 */	19,	/* modified TRP fluorotryptophane	*/
		/*  35 */	6,	/* modified GLY post-translational modification	*/
		/*  36 */	4,	/* modified GLU inhibited by dccd	*/
		/*  37 */	7,	/* modified HIS phosphorylated histidine	*/
		/*  38 */	19,	/* modified TRP beta-hydroxytryptophane	*/
		/*  39 */	13,	/* modified PRO 4-hydroxyproline	*/
		/*  40 */	8,	/* modified ILE allo-isoleucine	*/
		/*  41 */	9,	/* modified LYS carbamylated lysine	*/
		/*  42 */	9,	/* modified LYS lysine-pyridoxal-5'-phosphate	*/
		/*  43 */	9,	/* modified LYS methylated lysine	*/
		/*  44 */	14,	/* modified GLN 2-methyl-glutamine	*/
		/*  45 */	11,	/* modified MET s-oxymethionine	*/
		/*  46 */	7,	/* modified HIS post-translational modification	*/
		/*  47 */	9,	/* modified LYS methylated lysine	*/
		/*  48 */	11,	/* modified MET selenomethionine	*/
		/*  49 */	20,	/* modified TYR meta-tyrosine	*/
		/*  50 */	1,	/* modified ALA chemical modification	*/
		/*  51 */	7,	/* modified HIS n1-phosphonohistidine	*/
		/*  52 */	10,	/* modified LEU norleucine	*/
		/*  53 */	16,	/* modified SER o-acetylserine	*/
		/*  54 */	2,	/* modified CYS cysteinesulfonic acid	*/
		/*  55 */	2,	/* modified CYS hydroxyethylcysteine	*/
		/*  56 */	5,	/* modified PHE d-phenylalanine	*/
		/*  57 */	5,	/* modified PHE l-phenylalaninol	*/
		/*  58 */	13,	/* modified PRO thioproline	*/
		/*  59 */	16,	/* modified SER covalently bound inhibitor	*/
		/*  60 */	16,	/* modified SER covalently bound inhibitor	*/
		/*  61 */	2,	/* modified CYS thio methylated cysteine	*/
		/*  62 */	16,	/* modified SER o-benzylsulfonyl-serine	*/
		/*  63 */	16,	/* modified SER phosphoserine	*/
		/*  64 */	2,	/* modified CYS post-translational modification	*/
		/*  65 */	2,	/* modified CYS s-nitroso-cysteine	*/
		/*  66 */	2,	/* modified CYS dioxyselenocysteine	*/
		/*  67 */	20,	/* modified TYR tyrosine-o-sulphonic acid	*/
		/*  68 */	16,	/* modified SER serine vanadate	*/
		/*  69 */	19,	/* modified TRP ne1-formal-tryptophan	*/
		/*  70 */	19,	/* modified TRP aza-tryptophan	*/
		/*  71 */	20,	/* modified TYR 3,5-diiodotryrosine	*/
		/*  72 */	20,	/* modified TYR sulfonated tyrosine	*/
		/*  73 */	20,	/* modified TYR 3-fluorotyrosine (see remark 600)	*/
		/*  74 */	0,	/* blank - can be filled if any	*/
		/*  75 */	0,	/* unexpected MODRES records are	*/
		/*  76 */	0,	/* encountered; up to 10 unknown	*/
		/*  77 */	0,	/* types of modified residue can	*/
		/*  78 */	0,	/* be handled, and added to those*/
		/*  79 */	0,	/* listed above; the number of	*/
		/*  80 */	0,	/* types of modified residue is	*/
		/*  81 */	0,	/* therefore dynamic, and held in*/
		/*  82 */	0,	/* variable:	  heterogen_modres	*/
		/*  83 */	0	/* 	*/
	};	

	/* 27-6-1
	A list of known 3-letter codes of solvent heterogens. These are ignored.
	Their only purpose is so that unknown heterogens can be signalled.
	Unknown non-solvent heterogens (i.e., present as part of the
	polypeptide chain) should be added to the local_heterogen3[] array,
	with HETEROGENS increased accordingly (see preproc.h).
XXX many aren't 'solvent' - change name of this array
	*/

	static char *local_solvent_heterogen3[NON_AA_HETEROGENS] =
		{
		/* these are very roughly in order of abundance in PDB files */
		"Hoh","Dod","Acn","Oar","Oph","Rph"," ca"," na","Hem","Nad",
		"Nag","So4","Fad","Adp","Glc","Gal","Man","Gtp","Bcl","Plp",
		"Nap","Atp","Hec","Po4","Fmn","Li1","Ndp","Gdp","Anp","Dpg",
		"Mes","Cyc","Lmz"," zn","C78","Gol","Mal","Fs4","Ptr","Bph",
		"2dp","2gp","Seo","Iph","Llp"," cl","U10","Dhe","Cob","Peb",
		"Coa","Imp","Nga","Ret","Bme","Mle","Xys","Gsp","Ap5","Srm",
		"Pte"," mg","Sah","Pmp","Myr","Amp","Pca","Hed","Flc","Bog",
		"Orn","5iu"," mn","Aba","Gnp","Act","Ump","Lda","Mpd","Cmp",
		"Bak","Sia","5gp","Arc","Inb","Ppg","Ptt","Hxc","Bar","Fok",
		"Ntp","  k","Pee","Oes","Ace","Cap","Cnc","Cb3","Fuc","Gtn",
		"Blv","Ben","Cit","2pe","Mae","Net","A3p","Xv6","Caa","Edo",
		"Bru","Dmp"," cu","Bai"," pc","Arg","Thp","Im1","Gtx","Pga",
		"Mtx","Pad","Pqq","Kcx","Tch","35g","Epy","Ddh","Bla","Bab",
		"Ata","Ham","Nah","Tpp","Bah","Pop","Oxl","Sph","U0e","Lpc",
		"Faa","Est","Fbp","Hea","Epe","Ppd","Bmt","U2g","Ctp","Thk",
		"Inn","Gtt","Nth","4ip","Tpq","Dip","3pg","Imo","Dct","Cro",
		"Esi","Cnd","3in","Aco","Cgs","Ap1","Psu","H2u","Naq","Ctc",
		"Pfg","Gcg","Esp","Tht","345","Ilg","Imu","Rs1","Nae","Tr1",
		"Asl","Tad","Mgd","Pub","Irp","Pps","Gar","Phb","1in","Nax",
		"Ium","Pcp","Sei","Acr","Dad","Sae","Cto","Plm"," pi","Tnb",
		"Bin","Squ","Pam","Pep","Flf","Gln","1pn","Ad2","Pth","G6p",
		"Ssb","Dcm","Upa","Gr3","Bcd","Pdp","Gai","P3m","Prh","Dcf",
		"Ama","Apa"," u1","Eoh","Pp7","Men","Nh2","Sty","Upg","Fes",
		"Dps","Mva","Ihp","Bba","Fka","Cgu","Nhr","Nhs","Fer","Lys",
		"Sig","Fs3","Stu","Pcg"," y3"," cd","Co3","Tm6","Atf","Dms",
		"Rea","Rap","Pal","Gdb","Csd","Cst","Fol","C8e","Pap","Atr",
		"Sap","146","Fra","Asf","Nle","Ipa","Trs","Dio","Nmx","Rbu",
		"2fp","Glu","Ths","Alf","Mln","E64","Cbs","Fk5","Ctr","Fmb",
		" co","Dal","Cxm","Its","Cho","Bsp","Sas","Ucp","Tfp","Hop",
		"Tm5","Sta","Fcb","Hci","Umg","T6p","Rbf","Sam","Sfg","Bam",
		"Fca","Sin","Pi7","Hax","5hp","Gr4","Sad","Tmq","Ccn","Ndc",
		"Ara","Tar","Pms","Tet","Udp","Hpa","Ipp","Su2","Pi1","Pi5",
		"B3i","Phc","Ply","Fmt","Etf","Acv","Rs2","Ndo","Bip","1c5",
		"Cme"," fe"," t3","Bpm","Gcp","Ctt"," hg","Gts","Nbn","Str",
		"Mk1","Sar","Eta","Aeb","Gua","Cha","Tsn","Bly","Mp1","Cod",
		"No3","Cpr","Hby","Nov","Dom","M7g","Hii","Zra","In5","Leo",
		"Cbd","Mcn","Cbz","Mhm","Pcm","Mmp","Pi6","Pp8","In6","Cac",
		"Cch","2an","Arq","Mo6","Sdk","Pme","Tes","T19","Ina","Md2",
		"Inp","Rtl","T42","Nmb","Ih3","Psi","Gum","Ura","  u","Npo",
		"Omt","Hda","5cm","Gmy","Ba1","Oxp","In3","Ah1","Gw3","Mdl",
		"Bho","4su","Gco","Pi4","Fkp","Eno","Azi","Ih2","Bis","Amy",
		"Mno","Pos","Igp","Adn","Fty","Cbl","Tpt","Mga","Alr","Eaa",
		"Cyg","Bm9","Pi2","I10","Dar","8ig","Pdc","Pi3","Phe","Fip",
		"Dpn","Zn3","Opa","13p","Gla","Glb","Ih1","Adr","Aho","Bez",
		"Mdc","Pyx","Mgt","Bm2","Ure","Css","Acy","T16","Dcb","1bh",
		"Dgl","Fct","Cmo","N25","Ct3","061","Tcl","Adc","T29","Ald",
		"Bio","C2f","Haz","Scn","Bzo","Fum","Int","Tmf","Iod","Fya",
		"Bo3","Npe","He1","Csw","Ope","Seb","Dsy","Ppb"," fa","F6p",
		"Dzf","Hom","Quo","Gto","Ipd","Ufp","Amd","Chi","Atg","Hag",
		"Flx","Gpx","R56","Mai","Zya","Ags","Iva","Pbb","Arb","Mpp",
		"Dmf","Dpm","Azl","Fle","Sep","Dtp","Sga","Fos","Pse","Roi",
		"Gt9","Mia","Gas","Sba","Erg","S80","Tpi","Cpb","Bsi","Avp",
		"Btb","Npp","Emr","3ga","G21","Gnh","Bef","Iop","Dao","Ap2",
		"W11"," xe","Cna","Tze","Sb2","Imd","Can","Tml","Cbx","Tnp",
		"A80","Ict","Gp8","Hmi","S27","In4","In7","Bcs","Cef","5pv",
		"Isa","Sb4","Rdc","Sbl","Pla","Idc","Ppe","P4p","R13","Prd",
		"Pxg","Al7","Al8","Trc","Lla","Edr","Bnz","Arm","Gep","Nox",
		"Nph","Opg","Bzb","1py","Put","Fup","G7m","I3p","Ipm","Tol",
		"Al2","Bhc","Hv8","Alz","Trz","Bm5","Ffc","Sb5","Hep","Fmp",
		"Pt1","Al1","Al6","Nem","Pba","Bma","Dhb","Cld","Clq","Pgy",
		"Sb6","Re9","Pls","P4c","R12","2ph","Acd","Ptp","Af1","Tk4",
		" kr","Cam","Xan","Gun","Al3","Al4","Al9","Cep","Ncr","Tpo",
		"4tb","Bnn","Boc","Tyi","C2o","Pp3","Esx","Ac1","Abi","Ova",
		"Ron","Lea","5mu","Daf","Al5","Bm1","Dfp","Aop","Apb","Urs",
		"Bo4","Dix","Diy","2bl","Mol","Pgc","Pho","Hbd","Pic","Mta",
		"Sbt","Flt","Tca","Tck","G16","Foe","2og","Gp6","Thc","Pto",
		"Ade","5it","Im2","Bcc","Img","Vg1","Be2","Gsh","Nar","Daa",
		"Dag","Bic","1da","761","Fe2","2dt","Hav","Nos","Cra","Enh",
		"Shh","Ets","Iic","Imh","Hpr","5nc","9di","Tps","Cfc","Neo",
		"Amu","Ngl","Of1","Dfx","Fba","Dhp","Cms","No2","Dme","Doc",
		"Msb","Aya","Mts","Ibz","Hex","Zno","Gma","Ota","Cyn","Sle",
		"Ads","Sme","Bdn","Cas","Ip4","Aip","Oaa","Oai","Ap3","Hc0",
		"Npl","Sbi","Hf1","Ctd","Cyq","Pts","Naf","Ipl","Stl","Cgn",
		"Oct","Ans","Dle","Mpi","Mpl","Hc1","Phs","Nop","Btn","Ias",
		"Tcr","Dur","Prl","Lac","Tia","Bct","Gse","Pxy","Ioh","Cer",
		"Pa5","Ddu","Dg2","Far","Brc","Lta","3ib","Bzs","2ma","1pa",
		"Dum","Tha","Pvl","Aib","Htr","Aln","Oba","Trn","Itm","Mic",
		"Dfo","2ap","Iyr","Ph2","3ap","4ap","C6c","Csb","Ica","Set",
		"Ept","Nvi","Mzm","R1p","Bcy","Htp","Cea","Cem","Cff","Dfi",
		"Nic","Mnp","C5c","Dnn","Dnp","Azm","Ici","Pph","Aan","3mp",
		"Gpe","Mbh","Spa","Tra","Ams","Api","Clt","Mo2","Feo","Moh",
		"Ehp","Gcu","Mpj","Haa","Dml","Omd","Buc","But","Mto"," au",
		"Po3","1nb","Hho","Ppt"," gd","3mt","Hmf","Hmr","Hpe","Hpg",
		"Ins","Ip3","Ahc","Mf3","Sva","Tsu","Aph","Efc","M1a","Gam",
		"San","Csp","Pmb","Eoa","Yrr","Cxp","Tem","Cyl","4nc","Be1",
		" ni","Hv7","Mfu","Alc"," sm","Vo4","Kph","Pbp","Dgn","Mo4",
		"Nlp","Pge","Ram","Mph","N4b","Msh","2hp","Tbu","Gld","Ppi",
		"Frd","Ach","Aea","G6d","Spd"," no","I4b","Ips","Tos","Hto",
		"Akg","Iso","Den","Mlt","Oho","Dmr","Sac","Ggl","Sch","Fla",
		"1mc","Hgb","Idm","Bzf","Thz","Slz","Soa","Ind","Aj3"," oh",
		"Met","Bhd","Ch2","Unx","Ni1","Wo4","Dmi","2ep","Bu1","Enc",
		" ar"," cp","Fms","Dtt","Abn","2pi","Af3","Oxe","Oxy","So3",
		"Snc","Snn","Mbr","Tmt","Tmz","Naw","Hsm","Dab","Das","Hv5",
		"Dce","Anl","Fco","Kth","2ez"," br","1mz","Nva","Cyh","24t",
		"Hma","25t","Oc6","Ch3","Amt","Ni2","Mmc","Ars","2fu","Emc",
		"Dpp","Seg","2mz","2no","Tfa","Cys","Ety","Rng","Oxm","So2",
		"Pyr","Cde"," ox","Oc1","Oc5","Alm","Sul","Ofo","Lom","Asp",
		"Moo","Bro","Hae","Npn","Cu1","Peo","C1o","Nme","Fgl","Sea",
		"Sfo"," gb","Owq","Cfo","Ohe","Mnc"," yb","Na2","Oc2","Oet",
		"Per","  f"," ce"," cr","Cya"," ir","Gte","Nao"," pb"," pt",
		" sr","Nh4","  o","  s"," ag","Val"," ho","H2s","4mo"," o ",
		" me"," mo"
		};


	/* Set the values of the arrays flagname, parname, amino_acid1,
	amino_acid3, heterogen3, heterogen_no_to_amino_acid and
	local_solvent_heterogen3, by
	setting each element (a pointer) to point to the same place as
	the equivalent pointer in the local_.... arrays. Its done like
	this because explicitly setting the contents of the non-local
	arrays upsets some compilers, if its done more than once. Its
	done more than once because the flagname, parname etc arrays
	must be defined (by way of the #include of global.h) for each
	source file which uses them.

	N.B. heterogen_no_to_amino_acid is an array of pointers to int,
	the others are arrays of pointers to char */

	for (i = 0; i < FLAGS; i++)
		/* flagname[i] & local_flagname[i] are pointers;
			set the correct value of the strings to which
			flagname[i] point, by resetting the pointer */

		flagname[i] = local_flagname[i];

	for (i = 0; i < PARS; i++)
		/* see above, re flagname */
		parname[i] = local_parname[i];


	for (i = 0; i < AMINO_ACIDS; i++)
		{
		/* see above, re flagname */
		amino_acid1[i] = local_amino_acid1[i];
		amino_acid3[i] = local_amino_acid3[i];
		}

	for (i = 0; i < HETEROGENS_MAX; i++)
		{
		strcpy(heterogen3[i],local_heterogen3[i]);
		map_heterogen_no_to_amino_acid[i] = local_map_heterogen_no_to_amino_acid[i];
		}

	for (i = 0; i < NON_AA_HETEROGENS; i++)
		{
		solvent_heterogen3[i] = local_solvent_heterogen3[i];
		}

	/* finally, initialize the global variable n_heterogens */

	n_heterogens = HETEROGENS;
	}
