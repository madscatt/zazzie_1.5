/* program: xtal2sas */

/*    written  24apr97 by Susan Krueger, NIST to run under unix */
/*    modified program agpr2h2o by Glenn Olah, LANL, originally  */
/*    written in 1992 for a PC using Borland C++ */

/*    Calculates the P(R) function given the .PDB file (crystal  */
/*    coordinates) using a MONTE CARLO method.  Spheres are built around */
/*    each amino acid centered about its alpha carbon. Calculates I(Q) */ 
/*    from the P(R) function.  */

#include <stdio.h>
#include <ctype.h>
#include <iostream>			    
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fcntl.h>
#include <time.h>     /* used to generate seed for random generator */

#define random(x) ( (random()/(float)((1<<31)-1)) * x ) 

#define pi 3.14159268
#define NA 6.0223e23

/* coordinates of data (x,y,z) and scattering at that coordinate (b) */
float x[1000000], y[1000000], z[1000000], b[1000000];
int FLAG_AA[1000000];

/* PDB input file information */
char MARKER[6], ATOM_NAME[7], RESIDUE_TYPE[5];
int RESIDUE_NUMBER, MOLECULE_NUMBER, unknown;
long ATOM_NUMBER;
float OCCUPANCY, Bfactor;

/* files */
FILE  *inf, *outf, *tmpf;
char outpfile[4][30], inpfile[2][50];
char fpart[30], tempfile[35];
int n2nd, nfile;

/* parameters pertaining to residues */
/* fraction of H-D exchange (fexch), number of exchangeables */
/* (nexch_AA), bD and bH added by SK */
char radiation_type[1], HorD[2];
int NUM_AA[21], NUM_ATOM[6], n, num, max_AA, ndat, hflag;
float vlarge_AA[21], tb=0., tV=0., tmw=0.; /* amino acid info. */
float mcV=0., vp, rho, fexch, bD=0.6671, bH=-0.3742; /* amino acid info. */
float rH2O, vH20, pD2O, Rg, delta, PR[1000], dmax, boxdiag, rprobe;
float xl, yl, zl, rl; /* box size around protein */
float xmin, ymin, zmin, xmax, ymax, zmax; /* used to determine box size */
float xmin0, ymin0, zmin0, xmax0, ymax0, zmax0; /* molecule dimensions */

/* parameters pertaining to I(Q) */
float Io_meas, qmax;
int qnum;

/* parameters pertaining to calculating P(R) and monte carlo */
long NTRIALS, SLAMITS;
long HITS, MLOOPS, DOUBLES, HDOUBLES;

/* 6 ATOM TYPES (1-6) */
char CHAR_ATOM[6]={'H','C','N','O','S','P'};
short e_ATOM[6]={1,6,7,8,16,15}; /* electron number of each type of atom */

/* van-der-Waals volumes for each ATOM TYPE */
float vdWs_volume_ATOM[6]={5.575,15.599,11.494,10.306, 24.838,27.833};
float vdWs_radius_ATOM[6]={1.1,1.55,1.40,1.35,1.81,1.88};

/* AMINO ACID information arrays. all information taken from */
/* Jacrot, B. "Determination of Molecular Weight by Neutron Scattering", */
/* Biopolymers Vol. 20,p. 2413-2426 (1981). */
/*char *AA[]={"GLY","ALA","VAL","LEU","ILE","PHE","TYR","TRP","ASP","GLU",*/
char *AA[]={"GLY","ALA","VAL","LEU","ILE","PHE","TYR","TRP","ASP","GLU", "SER","THR","ASN","GLN","LYS","ARG","HIS","MET","CYS","PRO","WAT"};

/* scattering length density of each amino acid */
/* x-ray scattering */
float bX_AA[21]={30.,38.,54.,62.,62.,78.,86.,98.,59.,67.,46.,54.,60.,68.,
       71.,85.,71.5,70.,54.,52.,10.};
/* protonated protein, neutron scattering in D2O */
/* these bD_AA values do not include exchange.  The right values are */
/* calculated later. SK */
/*float bD_AA[21]={2.769,2.686,2.520,2.437,2.437,5.180,6.802,8.118,4.886,4.803,
       4.308,4.224,6.580,6.497,5.752,9.714,6.521,2.805,4.013,2.227,1.914};*/
float bD_AA[21];
/*These valuse are updated on Dec 10, 2002 SKG and SK to reflect H/D exchange as follows:
  fexch = percent of exchangables. 
  b(D2O)=b(H2O)+fexch*nexch*(bD-bH).
  The old valuse are commented out, but left in the program. */
/* protonated  protein, neutron scattering in H2O */
float bH_AA[21]={1.728,1.645,1.479,1.396,1.396,4.139,4.719,6.035,3.845,3.762,
       2.225,2.142,3.456,3.373,1.586,3.466,4.959,1.764,1.930,2.227,-0.168};
/* deuterated protein, neutron scattering in D2O */
float bDD_AA[21]={4.85,6.852,10.854,12.850,12.850,13.51,14.09,16.45,8.010,
       10.01,7.432,9.432,9.704,11.70,15.12,17.00,11.73,11.14,7.137,9.516,1.914};
/* these bDH_AA values do not include exchange.  The right values are */
/* calculated later. */
/* deuterated protein, neutron scattering in H2O */
/*float bDH_AA[21]={3.81,5.811,9.813,11.809,11.809,12.217,12.007,14.367,6.969,
       8.969,5.349,7.349,6.580,8.58,10.95,10.95,10.17,10.10,5.05,9.516,-0.168};*/
float bDH_AA[21];

/* Volume and M.W. of each amino acid - taken from JACROT's paper */
float volume_AA[21]={66.4,91.5,141.7,167.9,168.8,203.4,203.6,237.6,113.6,140.6,
	   99.1,122.1,135.2,161.1,176.2,180.8,167.3,170.8,105.6,129.3,29.8};
float radius_AA[21]={2.512,2.7954,3.2342,3.4223,3.4284,3.6482,3.6495,3.8423,3.0044,
  3.22576,2.8707,3.0776,3.1839,3.3754,3.4778,3.5078,3.4182,3.4419,2.9322,3.1369,1.9276};
float MW_AA[21]={57.,71.,99.,113.,113.,147.,163.,186.,114.,128.,87.,101.,114.,
	   128.,129.,157.,136.5,131.,103.1,97.,18.};

/*skg: July 30 2001:  Flag to identify polar and nonpolar residues */
int Polar_AA[21]={0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,0,0,0,0,1};
/*SK: number of exchangeables for each residue */
float nexch_AA[21]= {1.,1.,1.,1.,1.,1.,2.,2.,1.,1.,2.,2.,3.,3.,
   4.,6.,1.5,1.,2.,0.,0.};
   

/* subroutine declaration */
void read_it();
void print_it();
int  line_read(float *xx, float *yy, float *zz);
void dumpout();
void protein_boxsize_AA(int max_AA);
void flag_it(char c, short *ix);
void calc_PR(int max_AA); /* calculate P(R) function */
void calc_Rg(); /* calculate Rg from P(R) function */
void calc_I(); /* calculate I(q) .vs. q */
float simpson(float *xx, float delta, int nmax); /* integration */

/*.................................................................*/
int main(void)
{

/* increment variables */
   int i, k, m, mk, max_ATOM, j_ATOM, j_AA, ic ;
   short iflag;
   float tmpx, tmpy, tmpz; /* temporary coordinates read in */

   printf("*** Determine P(R) function from *.PDB file ***\n");

/* read in input/output filenames and other info. */
   read_it();

/* increase volumes of each amino acid in order to guarantee that single */
/* HIT monte carlo volumes are correct... 
   SKG:  one doesn't really need to increase the volume, so vH20=0 */
   vH20=0.0;
   for (i=0;i<20;i++) { vlarge_AA[i] = volume_AA[i]*(1.+vH20); }
   for (i=0;i<20;i++) { radius_AA[i] = pow(3.*vlarge_AA[i]/4./pi,0.333333); }


/***************************/
/* calculate bD_AA and bDH_AA, taking into account H-D exchange, SK */   
   if (toupper(radiation_type[0])!='X') {
      for (k=0;k<20;k++)  {
         bD_AA[k] = bH_AA[k] + fexch*nexch_AA[k]*(bD-bH);
         bDH_AA[k] = bDD_AA[k] + fexch*nexch_AA[k]*(bH-bD); 
      }
      bD_AA[20]=1.914;  /*Water must not use exchangable equation!*/
      bDH_AA[20]=-0.168; /*Same as above*/
   }
   
/* initialization of loop counter */
   j_ATOM=0; j_AA= 0; max_AA= 0; max_ATOM= 0;
   tV= tmw= 0.;
   for (k=0;k<1000000;k++) { x[k]= y[k]= z[k]= b[k]= 0.; FLAG_AA[k]= 0; }
   for (k=0;k<21;k++) {NUM_AA[k]= 0; }
   for (k=0;k<6;k++) {NUM_ATOM[k]= 0; }

/* read from main file */
   for (k=0;k<50;k++) {
      if (inpfile[0][k]=='.') break;
      else fpart[k]=inpfile[0][k];
   }
   fpart[k]='\0';

   sprintf(tempfile,"%s.tmp",fpart);
   
   n2nd= 1;
   /*   if (inpfile[1][0]!=NULL) n2nd= 2; */
   for (mk=0;mk<n2nd;mk++) {

/* is the protein protonated or deuterated ? */
      if (toupper(radiation_type[0])!='X') {
         fflush(stdin); 
         printf("Is protein (%s) protonated or deuterated? <H/D>: ",
		inpfile[mk]);
         do {
	    HorD[mk]= getchar();
	 } while (toupper(HorD[mk])!='H' && toupper(HorD[mk])!='D');
      } else HorD[mk]= 'H';

      printf("Reading input file %s\n",inpfile[mk]);
      if ((inf=fopen(inpfile[mk],"r"))==NULL) { 
	 printf("\nFile open failure"); exit(2); }


   std::cout << "READING DATA FROM FILE" << std::endl ;   
   std::cout << "READING DATA FROM FILE" << std::endl ;   
   std::cout << "READING DATA FROM FILE" << std::endl ;   

/* main DO LOOP for reading in data */
      do {
         if (line_read(&tmpx,&tmpy,&tmpz) == 0) continue;
         if ((j_ATOM==0 && n==0) || ((j_ATOM-1)==max_ATOM && n==1))
            i= RESIDUE_NUMBER;
/*
	   printf("\nMARKER = %s\n", MARKER) ;
	   printf("ATOM NUM = %li\n",ATOM_NUMBER) ;
	   printf("ATOM NAME = %s\n", ATOM_NAME) ;
	   printf("RESIDUE TYPE = %s\n", RESIDUE_TYPE) ;
	   printf("MOL NUM = %i\n",MOLECULE_NUMBER) ;
	   printf("RES NUM = %i\n",RESIDUE_NUMBER) ;
	   printf("OCCUPANCY = %f\n",OCCUPANCY) ;
	   printf("Bfactor = %f\n",Bfactor) ;

	   printf("tmpx = %f\n",tmpx);
	   printf("tmpy = %f\n",tmpy);
	   printf("tmpz = %f\n",tmpz);
*/
       
/* sum up number of each kind of atom */
         flag_it(ATOM_NAME[0],&iflag);
         NUM_ATOM[iflag]++;

/* determine M.W. and volume from amino acid */
         if ((ATOM_NAME[0]=='C') && (ATOM_NAME[1]=='A')){
	    i++; j_AA++;
	    x[j_AA-1]= tmpx; y[j_AA-1]= tmpy; z[j_AA-1]= tmpz;
	    for (k=0;k<21;k++) {
	       if (strcmp(RESIDUE_TYPE,AA[k])==0) {
	          FLAG_AA[j_AA-1]= k;
	          if (toupper(radiation_type[0])=='X')
	             b[j_AA-1]= (bX_AA[k]/volume_AA[k] - bX_AA[20]/
			      volume_AA[20]);
		  else
		  if (toupper(HorD[mk])=='H')
		     b[j_AA-1]= ((1.-pD2O)*(bH_AA[k]/volume_AA[k] - bH_AA[20]/
			       volume_AA[20]) + pD2O*(bD_AA[k]/
			       volume_AA[k] - bD_AA[20]/volume_AA[20]));
		  else
		     b[j_AA-1]= ((1.-pD2O)*(bDH_AA[k]/volume_AA[k] - bDH_AA[20]/
			       volume_AA[20]) + pD2O*(bDD_AA[k]/
			       volume_AA[k] - bDD_AA[20]/volume_AA[20]));
		  NUM_AA[k]++; tV += vlarge_AA[k]; tmw += MW_AA[k];

/* pertinent only for neutrons */
/*************************************/
/* this next line isn't used, from what I can see.  SK */
/*	          tb += ((1.-pD2O)*bH_AA[k] + pD2O*bD_AA[k]);*/
		  break;
	       }
	    }
	 }
	 else{
	/*
		printf("ATOM NAME == %s\n", ATOM_NAME) ;*/
		int friendship = 1 ;
	 }
	 j_ATOM++;
      } while (!feof(inf));   /* end each file read DO LOOP */
      
      fclose(inf); /* close input file */
      max_ATOM= j_ATOM-1; max_AA= j_AA-1; 
      //max_AA = 1E6 ; 

   }  /* end of FOR LOOP k through input *.PDB files */



   std::cout << "DONE READING DATA FROM FILE" << std::endl ;   
   std::cout << "DONE READING DATA FROM FILE" << std::endl ;   
   std::cout << "DONE READING DATA FROM FILE" << std::endl ;   


   printf("# ATOMS READ       :");    printf("%5d",j_ATOM);
   printf("    # AMINO ACIDS READ :");   printf("%5d\n",j_AA);


/* vary the output file for statistics - files *.iq and/or *.pr should be */
/* averaged by the user after running program: xtal2sas.cpp */
   for (m=0;m<nfile;m++) {
      sprintf(outpfile[0],"%s_%d.pr",fpart,m+1);
      sprintf(outpfile[1],"%s_%d.inf",fpart,m+1);
      sprintf(outpfile[2],"%s_%d.iq",fpart,m+1);
      sprintf(outpfile[3],"%s_%d.crd",fpart,m+1);
      printf("\n***LOOP # : %d***\n",m+1);
      print_it(); /* write out info. to screen for observation and flags */

      
/* calculate size of BOX the protein fits in lased on the vector positions */
/* of all atoms */
      protein_boxsize_AA(max_AA);

/* calculate the correlation function */
      calc_PR(max_AA);

/* calculate Rg from P(R) function */
      calc_Rg();

/* calculate I(q) .vs. q */
      calc_I();

   printf("11. # data points in R-space = %3d\n",ndat);      
   printf("12. Rg = %f\n",Rg);
      
/* dump correlation function to output file */
      dumpout();
   }  /* end of FOR loop through output files */

}


/*.........................................................................*/
/* read in some required info */
void read_it(void)
{
   
/* open up main file containing names of the input files to be processed */
   fflush(stdin); 
   printf("PDB INPUT FILE : ");
   gets(inpfile[0]); 
   printf("PDB input file = %s\n",inpfile[0]);

/* scattering with x-rays or neutrons ? */
   fflush(stdin); 
   printf("X-RAYS or NEUTRONS (X/N): ");
   do { radiation_type[0]= getchar(); 
   }while (toupper(radiation_type[0])!='X' &&
   toupper(radiation_type[0])!='N');

   qnum= 400; 
   fflush(stdin); 
   printf("Number of I(Q) values: ");
   scanf("%d",&qnum);

   qmax= 0.4; 
   fflush(stdin); 
   printf("Maximum Q value: ");
   scanf("%f",&qmax);
   
   Io_meas= 1.0; 
   fflush(stdin); 
   printf("Intensity at I(0): ");
   scanf("%f",&Io_meas);

   hflag=0;
   fflush(stdin);
   printf("Flag = 0, water mass density that of bulk water.  \n");
   printf("flag = 1, surface polar residues hydrated with water mass density of bound water, SAAS caculation.  \n");
   printf("Flag = 2, all surface residues hydrated with water mass density of bound water.  \n");
   printf("Enter Hydration choice:  ");
   scanf("%d",&hflag);

   if (hflag !=0){
   rH2O= 0.;
   fflush(stdin);
   printf("hydration layer thickness (A): (0-5): ");
   scanf("%f",&rH2O);
   

   rprobe = 6.5;
   fflush(stdin);
   printf("Size of probe (A): (enter default as 6.5):  ");
   scanf("%f",&rprobe);
   }

/* input percent D2O in solvent */
   if (toupper(radiation_type[0])!='X') {
      fflush(stdin); 
      printf("Input percent D2O in solvent (0-100 percent): ");
      scanf("%f",&pD2O); pD2O= pD2O/100.;
   } else pD2O= 0.;

   
/*****************************/
/* added H-D exchange, SK */
/* input percent fraction of H-D exchange */
   if (toupper(radiation_type[0])!='X') {
      fexch=0.8;
      fflush(stdin); 
      printf("Input fraction of exchanged NH peptide (0.8-0.9): ");
      scanf("%f",&fexch); 
   } else fexch= 0.;

 

   nfile= 1;
   fflush(stdin); 
   printf("Number of output file iterations: (1-100): ");
   scanf("%d",&nfile);

   SLAMITS= 1000;
   fflush(stdin); 
   printf("Number of HITS: (100-10000): ");
   scanf("%ld",&SLAMITS);

}

/*.........................................................................*/
/* print info to screen for monitoring */
void print_it(void)
{
   fflush(stdin);
   if (toupper(radiation_type[0])=='X') printf("1. X-RAYS used\n"); else
      printf("1. Neutrons used\n");

/* percent solvent deuterated */
    printf("2. %%D2O= %7.3f , AA volume increase = %7.3f\n",100.*pD2O,vH20);

/* input file name */
   printf("3. PDB input file = %s -- %c\n",inpfile[0],toupper(HorD[0]));

/* output file which will contain correlation function */
   printf("4. Output: P(R) function = %s\n",outpfile[0]);
   
/* output file which will contain random point coords */
   printf("5. Output: Random point coordinates = %s\n",outpfile[3]);

/* output file which will contain information, M.W., volume ... */
   printf("6. Output: information file = %s\n",outpfile[1]);

/* output file which will contain Iq */
   printf("7. Output: I(q) .vs. q = %s\n",outpfile[2]);

/* input scaling factor for Io */
   printf("8. Io = %10.3f",Io_meas);

/* this is the parameter being looped through */
   printf(",  radius H2O(bound) = %5.2f\n",rH2O);

}

/*.........................................................................*/
/* read in one line (at a time) of the input file with a .PDB file format */
int line_read(float *xs, float *ys , float *zs)
{
   
   int ic, line_length;
   char str[10], LINE[100];

   ic= 0;

     do {
      fscanf(inf,"%c",&LINE[ic]);
   } while(LINE[ic++]!='\n');
   line_length= ic; 

   /*    fgets(LINE,80,stdin); line_length = strlen(LINE); */
   
    
   if (!((LINE[0] == 'A') && (LINE[1] == 'T') &&(LINE[2] == 'O')
	 &&(LINE[3] == 'M'))) return 0;

   for (ic=0;ic<4;ic++) MARKER[ic]= LINE[ic]; MARKER[5]='\0'; 
   for (ic=0;ic<6;ic++) str[ic]= LINE[ic+5]; str[6]= '\0'; ATOM_NUMBER= atol(str);
   for (ic=0;ic<4;ic++) ATOM_NAME[ic]=LINE[ic+12]; ATOM_NAME[4]='\0'; 
   for (ic=0;ic<3;ic++) RESIDUE_TYPE[ic]=LINE[ic+17]; RESIDUE_TYPE[4]='\0'; 
   for (ic=0;ic<2;ic++) str[ic]= LINE[ic+20]; str[2]='\0'; MOLECULE_NUMBER= atoi(str); 
   for (ic=0;ic<4;ic++) str[ic]= LINE[ic+22]; str[4]='\0'; RESIDUE_NUMBER= atoi(str); 
   for (ic=0;ic<8;ic++) str[ic]= LINE[ic+30]; str[8]='\0'; *xs= atof(str); 
   for (ic=0;ic<8;ic++) str[ic]= LINE[ic+38]; str[8]='\0'; *ys= atof(str); 
   for (ic=0;ic<8;ic++) str[ic]= LINE[ic+46]; str[8]='\0'; *zs= atof(str); 
   for (ic=0;ic<6;ic++) str[ic]= LINE[ic+54]; str[6]='\0'; OCCUPANCY= atof(str); 
   for (ic=0;ic<5;ic++) str[ic]= LINE[ic+61]; str[5]='\0'; Bfactor= atof(str); 
   if (line_length>=66) {
      for (ic=0;ic<3;ic++) str[ic]= LINE[ic+67];
      str[3]='\0'; unknown= atoi(str);
   }
   return 1;
}

/*.........................................................................*/
/* dumpout results - correlation function dumped and information file */
void dumpout(void)
{
   
   int ix, kx; float xd, yd;

   printf("Dumping out correlation function to file %s .\n",outpfile[0]);
   if ((outf=fopen(outpfile[0],"w"))==NULL)
     { printf("\nFile open failure"); exit(2); }

   for (ix=0;ix<ndat;ix++) {
      xd= float(ix)*delta;
      yd= PR[ix];
      fprintf(outf," %f        %f          0.0\n",xd,yd);
   }
   fclose(outf);

/* dump information to output file */
   printf("Dumping out information to file %s .\n",outpfile[1]);
   if ((outf=fopen(outpfile[1],"w"))==NULL)
     { printf("\nFile open failure"); exit(2); }
   if (toupper(radiation_type[0]=='X')) fprintf(outf,"1. X-rays used\n");
   else fprintf(outf,"1. Neutrons used\n");
   fprintf(outf,"2. PDB input file = %s -- %c\n",inpfile[0],toupper(HorD[0]));
   fprintf(outf,"3. I(q) .vs. q output file = %s\n",outpfile[2]);
   fprintf(outf,"4. Dimensions (A): x: %7.2f,%7.2f, y: %7.2f,%7.2f, z: %7.2f,%7.2f\n",xmin0,xmax0,ymin0,ymax0,zmin0,zmax0);   
   fprintf(outf,"5. Box diagonal length (A) = %7.2f\n",boxdiag);
   fprintf(outf,"6. Box Size (A): x: %7.2f , y: %7.2f , z: %7.2f\n",xl,yl,zl);
   fprintf(outf,"7. # data points in R-space = %3d\n",ndat);
   fprintf(outf,"8. Dmax = %7.3f , fraction D2O in solvent = %7.3f\n",dmax,pD2O);
   fprintf(outf,"9. M.W. = %7.2f  , volume (JACROT,increased) (A^3)= %7.2f\n",tmw,tV);
   fprintf(outf,"10. volume(s.v.=0.73) (A^3) = %7.2f\n",tmw*0.73*1.e24/NA);
   fprintf(outf,"11. Rg (A) = %6.2f , Io = %10.3f\n",Rg,Io_meas);
   fprintf(outf,"12. Level of hydration,%2d, rH2O (A) = %9.2f, probe radius = %9.2f\n",hflag,rH2O,rprobe);
   fprintf(outf,"13. # of atom types read:\n");
   for (kx=0;kx<6;kx++) fprintf(outf,"   %c : %d\n",
				CHAR_ATOM[kx],NUM_ATOM[kx]);
   fprintf(outf,"14. # of residues read:\n");
   for (kx=0;kx<21;kx++) fprintf(outf,"   %3d  %3s\n",NUM_AA[kx],AA[kx]);
   fclose(outf);

   printf("I(q) .vs. q was dumped out to %s .\n",outpfile[2]);
}

/*.........................................................................*/
/* determine minimum box size to fit protein in - calculate boxdiag, ndat */
void protein_boxsize_AA(int max_AA)
{
   
   int il;

   printf("Determine BOX SIZE to put protein in....\n");

/* determine maximum position */
   xmax= x[0] + radius_AA[FLAG_AA[0]];
   ymax= y[0] + radius_AA[FLAG_AA[0]];
   zmax= z[0] + radius_AA[FLAG_AA[0]];

/* determine minimum position */
   for(il=1;il<max_AA;il++) {
      if ((x[il] + radius_AA[FLAG_AA[il]])>xmax)
      	 xmax= x[il] + radius_AA[FLAG_AA[il]];
         xmax0=xmax;
      if ((y[il] + radius_AA[FLAG_AA[il]])>ymax)
         ymax= y[il] + radius_AA[FLAG_AA[il]];
         ymax0=ymax;
      if ((z[il] + radius_AA[FLAG_AA[il]])>zmax)
         zmax= z[il] + radius_AA[FLAG_AA[il]];
         zmax0=zmax;
   }

/* determine minimum position */
   xmin= x[0] - radius_AA[FLAG_AA[0]];
   ymin= y[0] - radius_AA[FLAG_AA[0]];
   zmin= z[0] - radius_AA[FLAG_AA[0]];

   printf(" max_AA = %i\n",max_AA) ;

   for(il=1;il<max_AA;il++) {
      if ((x[il] - radius_AA[FLAG_AA[il]])<xmin)
      	 xmin= x[il] - radius_AA[FLAG_AA[il]];
         xmin0=xmin;
      if ((y[il] - radius_AA[FLAG_AA[il]])<ymin)
	 ymin= y[il] - radius_AA[FLAG_AA[il]];
         ymin0=ymin;
      if ((z[il] - radius_AA[FLAG_AA[il]])<zmin)
	 zmin= z[il] - radius_AA[FLAG_AA[il]];
         zmin0=zmin;
   }

   /* detemine box size... */
/* slightly increase the box size so (WITH OUT A DOUBT) the molecule */
/* (including bound H2O) are completely in the box. */

   xl= fabs(xmax - xmin); yl= fabs(ymax - ymin); zl= fabs(zmax - zmin);
   xmin= xmin - xl/4. - rH2O*2.0; xmax= xmax + xl/4. + rH2O*2.0;
   ymin= ymin - yl/4. - rH2O*2.0; ymax= ymax + yl/4. + rH2O*2.0;
   zmin= zmin - zl/4. - rH2O*2.0; zmax= zmax + zl/4. + rH2O*2.0;
   xl= fabs(xmax - xmin); yl= fabs(ymax - ymin); zl= fabs(zmax - zmin);
   rl= 1./sqrt(xl*xl + yl*yl + zl*zl);
   boxdiag= 1./rl;
   if ((boxdiag*2)<1000.) ndat= boxdiag*2 +0.5; else ndat= 1000;
   delta= boxdiag/(float)ndat;
/*   printf("ndat=%ld, delta=%f\n",ndat,delta); */

/* print to screen */
   printf("9. Dimensions = x: %7.2f,%7.2f, y: %7.2f,%7.2f, z: %7.2f,%7.2f\n",
	  xmin0,xmax0,ymin0,ymax0,zmin0,zmax0);   
   printf("10. Box diagonal = %7.2f ; Box size = x: %7.2f , y: %7.2f , z: %7.2f\n",
	  boxdiag,xl,yl,zl);
}

/*.........................................................................*/
/* flag each atom type */
void flag_it(char A, short *iflag)
{
   
   switch (A) {
      case 'H' : *iflag=0; break;
      case 'C' : *iflag=1; break;
      case 'N' : *iflag=2; break;
      case 'O' : *iflag=3; break;
      case 'S' : *iflag=4; break;
      case 'P' : *iflag=5; break;
   }
}

/*.........................................................................*/
/* calculate the correlation function using a monte carlo method */
void calc_PR(int max_AA)
{
   
   int i, j, k, ib, ih, surface[1000000];
   float xo, yo, zo, x1, y1, z1, r1, r2, r, pr1, pr2, PRmax, rn;
   float xc[4];
   float xc1,yc1,zc1,wc1,xc2,yc2,zc2,wc2;
   long total_HITS;
   time_t t;

   for (k=0;k<1000;k++) { PR[k]= 0.; }

   MLOOPS= 1;
   for(k=0;k<MLOOPS;k++) {
      srand(1);srand((unsigned) time(&t)); /* seed the random generator */

/* open files */
      if ((tmpf=fopen(tempfile,"wb"))==NULL)
	{ printf("\nFile open failure"); exit(2); }      
      if ((outf=fopen(outpfile[3],"w"))==NULL)
        { printf("\nFile open failure"); exit(2); }

      printf("Generate random points ....\n");

      /*determine which residues are on surface (surface=1), only if requesting hydration*/
      
      for (j=0;j<max_AA;j++) {
	surface[j]=0;
	if (hflag != 0){
	for (i=0;i<j-4;i++) {
	  rn=sqrt((x[j]-x[i])*(x[j]-x[i])+(y[j]-y[i])*(y[j]-y[i])+(z[j]-z[i])*(z[j]-z[i]));
	  if (rn >= rprobe) surface[j]=1;
	  }
	for (i=j+4;i<max_AA;i++) {
	  rn=sqrt((x[j]-x[i])*(x[j]-x[i])+(y[j]-y[i])*(y[j]-y[i])+(z[j]-z[i])*(z[j]-z[i]));
	  if (rn >= rprobe) surface[j]=1;
	}
      }
      }

      printf("back from Generate random points ....\n");
	printf("SLAMITS = %li\n",SLAMITS);
/* loop for generating random points in box */
      NTRIALS= 0; HITS= 0; total_HITS= 0;
      do {
	 NTRIALS++;
         xo= (xl*(float) (random(50000))/50000.) + xmin;
         yo= (yl*(float) (random(50000))/50000.) + ymin;
         zo= (zl*(float) (random(50000))/50000.) + zmin;

/* see if random point is inside or outside protein */
         ib=ih= 0; pr1=pr2= 0.;
         for (j=0;j<max_AA;j++) {
	    x1= xo - x[j]; y1= yo- y[j]; z1= zo - z[j];
/* SKG: July 30, 2001:  Only want Polar residues to increase volume for H2O (working)*/
/* Phase II, take only surface residues to increase volume for H20 */
	    
            r1= sqrt(x1*x1 + y1*y1 + z1*z1);
	    if (surface[j] ==1) {
	      if (((hflag ==1 ) && (Polar_AA[FLAG_AA[j]] == 1 )) || (hflag==2)) {
            if (r1 <= (radius_AA[FLAG_AA[j]] + rH2O)) {
	       if (r1 <= radius_AA[FLAG_AA[j]]) {
		  ib++; total_HITS++; pr1 += b[j]; if (ib>=2) DOUBLES++;
/*                  printf("j=%ld, hits=%ld, ib=%ld, b[j]=%f, pr1 =%f\n", */
/*			 j,total_HITS,ib,b[j],pr1); */
	       
	       }else {
/* From S. J. Perkins paper (Eur. J. Biochem. 157, */
/* 169-180 (1986)), bound H2O has volume of 24.5 while bulk H2O has a */
/* volume of 29.8 */
		  ih++; if (ih>=2) HDOUBLES++;
		  if (toupper(radiation_type[0])=='X')
		     pr2 += (bX_AA[20]*(1./24.5 - 1./29.8));
		  else
		     pr2 += (((1.-pD2O)*bH_AA[20]+pD2O*bD_AA[20])*
			     (1./24.5 - 1./29.8));
	       }
	    }
	    }else {
	     if (r1 <= radius_AA[FLAG_AA[j]]) {
		  ib++; total_HITS++; pr1 += b[j]; if (ib>=2) DOUBLES++;
		  } 
	    } 
	    } else {
	          if (r1 <= radius_AA[FLAG_AA[j]]) {
		  ib++; total_HITS++; pr1 += b[j]; if (ib>=2) DOUBLES++;
		  }
	 }
	 }
/*         printf("ib = %ld, ih = %ld\n",ib,ih);   */
/*         printf("pr1 = %f, pr2 = %f\n",pr1,pr2); */
         if (ib!=0) {
	    ++HITS;
            xc[0]=pr1/float(ib); xc[1]=xo; xc[2]=yo; xc[3]=zo;
            fwrite(xc,sizeof(float),4,tmpf);
            fprintf(outf,"%8.4f  %8.4f  %8.4f  %8.4f\n",xo,yo,zo,xc[0]);
	 }
	 if (ih!=0 && ib==0) {
	    ++HITS;
            xc[0]=pr2/float(ih); xc[1]=xo; xc[2]=yo; xc[3]=zo;
            fwrite(xc,sizeof(float),4,tmpf);
            fprintf(outf,"%8.4f  %8.4f  %8.4f  %8.4f\n",xo,yo,zo,xc[0]);
	 }
//	  printf("HITS = %li\t",HITS);
      } while(HITS < SLAMITS);

      fclose(tmpf); /* temporary file closed to send end of file */
      
      printf("past loop for Ntrials....\n");
/* monte carlo integration for calculating volume.... */
      mcV = xl*yl*zl*float(HITS)/float(NTRIALS);

/* read temporary file */
      if ((tmpf=fopen(tempfile,"rb"))==NULL)
	{ printf("\nFile open failure"); exit(2); }

      printf("TOTAL HITS = %ld, TOTAL SHOOTS = %ld\n",HITS,NTRIALS);
      printf("Random point coordinates written to file %s .\n",outpfile[3]);
      printf("Calculating P(R) function ....\n");
      dmax= 0.;

      for(i=0;i<HITS;i++) {
	 fseek(tmpf,(long)(i*4*sizeof(float)),SEEK_SET);
         fread(xc,sizeof(float),4,tmpf);
         xc1= xc[1]; yc1= xc[2]; zc1= xc[3]; wc1= xc[0];

         for(j=i+1;j<HITS;j++) {
	    fread(xc,sizeof(float),4,tmpf);
            xc2= xc[1]; yc2= xc[2]; zc2= xc[3]; wc2= xc[0];
            r2= 0.;
            r2 += ((xc1-xc2)*(xc1-xc2));
            r2 += ((yc1-yc2)*(yc1-yc2));
            r2 += ((zc1-zc2)*(zc1-zc2));
            r= sqrt(r2);
            if (r>dmax) dmax= r;
            num= r/delta;
            PR[num] += float(wc1*wc2);
	 }
      }
      
      fclose(tmpf); fclose(outf);
   }

/* delete temporary file from hard disk */
   remove(tempfile);

   PRmax= 0.0;
   for (i=0;i<ndat;i++) {
      PR[i] /= (float)MLOOPS;
      if (PR[i]>PRmax) PRmax= PR[i]; 
   }
   i=ndat;
   do { i--; if (PR[i]!=0.) { ndat= i; break; } } while (i>0);
   PR[0]=0.; 
   for (i=0;i<ndat;i++)  {
      PR[i]= PR[i]/PRmax;
   }
   
}

/*.........................................................................*/
/* Calculation of Rg from the correlation function P(R) */
void calc_Rg()
{
   
   int i;
   float y[1000], normal;
   for (i=0;i<1000;i++) { y[i]= 0.; }
   normal= 0.;
   for (i=0;i<ndat;i++) y[i]= (float) i * delta * i * delta * PR[i];
   normal=simpson(PR,delta,ndat);
   Rg= sqrt(0.5*simpson(y,delta,ndat)/normal);
}

/*.........................................................................*/
/* Calculation of I(q) from the correlation function P(R) */
void calc_I()
{
   
   printf("%s\n", "iq data");
   int i, j;
   float y[1000], q, deltaq, Iq, Io_factor;
   if ((outf=fopen(outpfile[2],"w"))==NULL)
     { printf("\nFile open failure"); exit(2); }
   /*deltaq= qmax/qnum; q= (-deltaq);*/
   deltaq= qmax/(qnum-1); q= (-deltaq);
   for (j=0;j<qnum;j++) {
      q += deltaq;
      for (i=0;i<ndat;i++) {
	 if (((float)i*delta*q)!=0.)
	    y[i]= (float) PR[i]* sin(i*delta*q)/(i*delta*q) * delta;
	 else y[i]= (float) PR[i] * delta;
//	 printf("%f\t%f\n",q,y[i]);
      }
      Iq= 4.*pi*simpson(y,delta,ndat);
      if (j==0) Io_factor= Io_meas/Iq;
//      printf("%f\t%f\t%f\t%f\n",q,Io_factor,Io_meas,Iq);
      fprintf(outf,"%f   %f   0.\n",q,Io_factor*Iq);
   }
  
   fclose(outf);

}

/*.........................................................................*/
/* SIMPSON'S rule for integration */
float simpson(float y[], float dd, int N)
{
   
   float s1, s2, s3;
   int i;
   s1= 0.; s2= 0.;
   for(i=3;i<N;i+=2) {
      s1 += y[i];
      s2 += y[i-1];
   }
   s3= (y[0] + 4.* (s1 + y[1]) + 2.*s2 + y[N-1])/3.;
   return(s3);
}

