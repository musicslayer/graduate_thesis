//Matrix Builder
//Applies Newton's method to differential equations using a matrix.
//Corresponds to Vidal and Lovelace 2014 (Draft) and "New Treatment of the Pulsar Equation" (Master's Thesis).
//Created September 25, 2013 By Michael Joseph Vidal.

#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <io.h>
#include <math.h>
#include <string.h>
#include <time.h>

//Definitions for the types of simulations found in the program.
#define MONOPOLE 1
#define STANDARD 2
#define JETS 3
#define NULLSHEET 4

//This region of code serves as a control panel for the rest of the program. Comment or uncomment macros to get desired behavior.
typedef double T; //Choose your precision (float or double)

/*
Makes program read the matrix "J_INVERSE" from file instead of calculating it.
EVERY grid setting below must EXACTLY match the ones used in the file read.
READ may not work for Null Sheet cases.
*/
//#define READ 

#define SHORTCUT //Assumes constant Jacobian Matrix.
T V_FRAC=0.1; //Similar to relaxation parameter.

//These only affect the jets case. So far, only ckf_jets.c has this implemented.
#define JETS1 //Use boundary used by Takamori in jets case (this one is the typical choice).
//#define JETS2 //Use boundary used by Lovelace in jets case.

#define SMOOTH 1 //Amount of smoothing w.r.t. the light cylinder. SMOOTH 1 is the typical value.

#define NUM_ITERATIONS_MAX 100000 //Maximum number of iterations, regardless of convergence.

//Grid Parameters.
#define N_MAX_X 42
#define N_MAX_Y 42
#define N_LC 20
#define DX 0.05 //Make sure N_LC*DX=1.
#define DY 0.05
#define N_S_X 2 //Make sure the star is square i.e. DX*N_S_X=DY*N_S_Y.
#define N_S_Y 2

//Other simulation parameters.
#define RATIO 0.5 //Only used in TOTS case.
#define P_OP 1.0 //PSI EQUATORIAL.

int toggle=0; //Only used in double simulation cases.
T BETA=1.0; //Only used in jets case. (Defining it here allows us to change it mid simulation.

//Sizes of various objects used in the program.
const int SIZE_T=sizeof(T);
const int SIZE_T_X=(N_MAX_X-2)*sizeof(T);
const int SIZE_T_Y=(N_MAX_Y-2)*sizeof(T);
const int SIZE_T_X_Y=(N_MAX_X-2)*(N_MAX_Y-2)*sizeof(T);
const int SIZE_T_XY=(N_MAX_X-2)*(N_MAX_Y-2)*sizeof(T);
const int SIZE_T_XY_XY=(N_MAX_X-2)*(N_MAX_X-2)*(N_MAX_Y-2)*(N_MAX_Y-2)*sizeof(T);

//Files that are written to. Some are extra files that may be useful for debugging or for secondary data.
FILE *fp_r,*fp_z,*fp_a,*fp_b,*fp_j,*fp_j2,*fp_hhp,*fp_hhpPRIME,*fp_v,*fp_vxy,*fp_P_HHP,*fp_change;

T L_func[(N_MAX_X-2)];
T L_X[(N_MAX_X-2)*(N_MAX_Y-2)];
T L_Y[(N_MAX_X-2)*(N_MAX_Y-2)];
T CHANGE[(N_MAX_X-2)][(N_MAX_Y-2)];

T **A;
T *B;
T **J;
T **J_INVERSE;
T **HHP;
T **HHP_PRIME;
T *v0;
T *v1;
T *eq;
T **h;
T *p;

T **a1,**a2,**a3,**a4,**a5,**a6; 
T *bL1,*bL2,*bL3,*bL4,*bL5,*bL6;
T *bR1,*bR2,*bR3,*bR4,*bR5,*bR6;
T *bB1,*bB2,*bB3,*bB4,*bB5,*bB6;
T *bT1,*bT2,*bT3,*bT4,*bT5,*bT6;
T ***c;

void initialize(void);
void star(void);
__forceinline void create_matrix(void);
void writeToFiles(void);
__forceinline int RN(int,int);
__forceinline int CN_X(int);
__forceinline int CN_Y(int);
__forceinline T f(int,int);
__forceinline void hhpSet(T**,T**);
void equationBuilder(int);
void inverse(T**,int,T**);
T** Make2DTArray(int,int);
void resetCKF(void);
void resetTak(void);
void read(void);

/*
Case file choices. Include EXACLTY ONE of these to choose the case.
*/
//Single simulation cases.
#include "ckf_monopole.c"
//#include "ckf.c" //Pure CKF method from Contopoulos et al. 1999
//#include "ckf_jets.c" //Jets case
//#include "ckf_null.c" //Null Sheet case.
//#include "tak_monopole_test.c" //Monopole where we use the analytical expression for the answer.
//#include "tak_monopole.c" //Pure monopole case.
//#include "tak_monopole_jets.c" //Monopole case with jets.
//#include "tak.c" //Pure TOTS method from Takamori et al. 2012
//#include "tak_jets.c" //Jets case
//#include "tak_theory_jets.c" //Attempt to create the proposed theoretical jets case from Vidal and Lovelace 2014
//#include "tak_null.c" //Null Sheet case

/*
Double simulation cases. There are two general situations:
You start with TOTS method, and then go to CKF method,
or, You start with CKF method, and then go to TOTS method.
Note: Doing TOTS first will give CKF the grid of PSI values. Doing CKF first will give TOTS the grid of HHP values.
*/
//#include "ckf_monopole_jets.c" //CKF monopole, then Jets
//#include "ckf_tak.c" //CKF, then TOTS
//#include "ckf_null_tak.c" //CKF Null Sheet, then TOTS
//#include "tak_ckf.c" //TOTS, then CKF
//#include "tak_ckf_jets.c" //TOTS, then CKF Jets
//#include "tak_ckf_null.c" //TOTS, then CFK Null Sheet
//#include "tak_theory_jets_ckf.c" //TOTS theoretical jets case, then CKF

//The code that creates the matrix "A" and vector "b" for the linear system
#include "matrix_build.c"

//Code that calculates the inverse of a matrix.
#include "inverse.c"

//Code to solve the nonlinear equations.
#include "equation_solver.c" //Use with a "single simulation" case.
//#include "equation_solver_2case.c" //Use with a "double simulation" case.

//Additional code to reset part of the bottom boundary every iteration, which is neccessary only for Null Sheet cases.
#include "reset_ckf.c" //Only needed if ckf_null is an included case
#include "reset_tak.c" //Only needed if actual_null is an included case

//Code that reads the matrix "J_INVERSE" from file instead of calculating it.
#include "read_matrix.c"

int main(int argc, char *argv[])
{
	//Memory allocation and initialization.
	A=Make2DTArray((N_MAX_X-2)*(N_MAX_Y-2),(N_MAX_X-2)*(N_MAX_Y-2));
	B=malloc(SIZE_T_XY);
	J=Make2DTArray((N_MAX_X-2)*(N_MAX_Y-2),(N_MAX_X-2)*(N_MAX_Y-2));
	J_INVERSE=Make2DTArray((N_MAX_X-2)*(N_MAX_Y-2),(N_MAX_X-2)*(N_MAX_Y-2));
	v0=malloc(SIZE_T_XY);
	v1=malloc(SIZE_T_XY);
	eq=malloc(SIZE_T_XY);
	h=Make2DTArray((N_MAX_X-2)*(N_MAX_Y-2),(N_MAX_X-2)*(N_MAX_Y-2));
	p=malloc(SIZE_T_XY);

	HHP=Make2DTArray(N_MAX_X-2,N_MAX_Y-2);
	HHP_PRIME=Make2DTArray(N_MAX_X-2,N_MAX_Y-2);

	a1=Make2DTArray(N_MAX_X-2,N_MAX_Y-2);
	a2=Make2DTArray(N_MAX_X-2,N_MAX_Y-2);
	a3=Make2DTArray(N_MAX_X-2,N_MAX_Y-2);
	a4=Make2DTArray(N_MAX_X-2,N_MAX_Y-2);
	a5=Make2DTArray(N_MAX_X-2,N_MAX_Y-2);
	a6=Make2DTArray(N_MAX_X-2,N_MAX_Y-2);
	
	bL1=malloc(SIZE_T_Y);
	bL2=malloc(SIZE_T_Y);
	bL3=malloc(SIZE_T_Y);
	bL4=malloc(SIZE_T_Y);
	bL5=malloc(SIZE_T_Y);
	bL6=malloc(SIZE_T_Y);

	bR1=malloc(SIZE_T_Y);
	bR2=malloc(SIZE_T_Y);
	bR3=malloc(SIZE_T_Y);
	bR4=malloc(SIZE_T_Y);
	bR5=malloc(SIZE_T_Y);
	bR6=malloc(SIZE_T_Y);

	bB1=malloc(SIZE_T_X);
	bB2=malloc(SIZE_T_X);
	bB3=malloc(SIZE_T_X);
	bB4=malloc(SIZE_T_X);
	bB5=malloc(SIZE_T_X);
	bB6=malloc(SIZE_T_X);

	bT1=malloc(SIZE_T_X);
	bT2=malloc(SIZE_T_X);
	bT3=malloc(SIZE_T_X);
	bT4=malloc(SIZE_T_X);
	bT5=malloc(SIZE_T_X);
	bT6=malloc(SIZE_T_X);

	c=malloc(4*SIZE_T_X_Y);
	c[0]=Make2DTArray(N_MAX_X-2,N_MAX_Y-2);
	c[1]=Make2DTArray(N_MAX_X-2,N_MAX_Y-2);
	c[2]=Make2DTArray(N_MAX_X-2,N_MAX_Y-2);
	c[3]=Make2DTArray(N_MAX_X-2,N_MAX_Y-2);

	//Gives a rough idea of how much memory is used while the program is running.
	printf("SIZE T MBytes=%f\n",(T)SIZE_T/(1024*1024));
	//Big stuff
	printf("Big Stuff\n");
	printf("MBytes A=%f\n",(T)SIZE_T_XY_XY/(1024*1024));
	printf("MBytes J=%f\n",(T)SIZE_T_XY_XY/(1024*1024));
	printf("MBytes J-1=%f\n",(T)SIZE_T_XY_XY/(1024*1024));
	printf("MBytes h=%f\n",(T)SIZE_T_XY_XY/(1024*1024));

	//Small stuff
	printf("\nSmall Stuff\n");
	printf("MBytes B=%f\n",(T)SIZE_T_XY/(1024*1024));
	printf("MBytes v,eqBlock=%f\n",3.0*SIZE_T_XY/(1024*1024));
	printf("MBytes p=%f\n",(T)SIZE_T_XY/(1024*1024));
	printf("MBytes hhpBlock=%f\n",2.0*SIZE_T_X_Y/(1024*1024));
	printf("MBytes aBlock=%f\n",6.0*SIZE_T_X_Y/(1024*1024));
	printf("MBytes bBlock=%f\n",12.0*(SIZE_T_X+SIZE_T_Y)/(1024*1024));
	printf("MBytes cBlock=%f\n",4.0*SIZE_T_X_Y/(1024*1024));

#ifndef READ
	//Initialize everything to zero	
	for (int j=0;j<=(N_MAX_X-2)*(N_MAX_Y-2)-1;j++)
	{
		B[j]=0;
		for (int k=0;k<=(N_MAX_X-2)*(N_MAX_Y-2)-1;k++)
		{
			A[j][k]=0;
		}        
	}     	
#endif

	for(int x=0;x<N_MAX_X-2;x++)
	{		
		for(int y=0;y<N_MAX_Y-2;y++)
		{
			HHP[x][y]=0;
			HHP_PRIME[x][y]=0;
		}
	}

	//Prepare files. Note that some files are for secondary information that is largely unused.
	fp_r=fopen("r.dat","w+"); //R Values.
	fp_z=fopen("z.dat","w+"); //Z Values.
	fp_a=fopen("a.dat","w+"); //A Matrix.
	fp_b=fopen("b.dat","w+"); //B Vector.
	fp_j=fopen("j.dat","w+"); //Jacobian Matrix.
#ifdef READ
	fp_j2=fopen("j2.dat","r+"); //Inverse Jacobian Matrix.
#else
	fp_j2=fopen("j2.dat","w+"); //Inverse Jacobian Matrix.
#endif
	fp_hhp=fopen("hhp.dat","w+"); //HHP Matrix.
	fp_hhpPRIME=fopen("hppPRIME.dat","w+"); //D(HHP)/D(PHI) Matrix.
	fp_v=fopen("v.dat","w+"); //PSI Values (written in a line).
	fp_vxy=fopen("vxy.dat","w+"); //PSI Values (written in a grid).
	fp_P_HHP=fopen("P_HHP.dat","w+"); //PSI vs. HHP.
	fp_change=fopen("change.dat","w+"); //Matrix of change of PSI between iterations.

	//Exits the program immediately if a file cannot be opened, to prevent unnecessary time wasting.
	if(NULL==fp_hhp||NULL==fp_r||NULL==fp_z||NULL==fp_v)
	{ 		
   	 	printf("Error - Cannot open file.\n"); 
    		exit(0); 
	}	
	
	//Main program, separated into parts to show the user progress.
	printf("\n1\n");
	initialize(); //Case specifc discretization of the differential equation and boundary.
	printf("2\n");
	star(); //Case specific insertion of the star.
	printf("3\n");
	create_matrix(); //Creates the actual matrix "A" and vector "B".
	printf("4\n");
#ifdef READ
	read(); //Reads J_INVERSE.
#endif
	equationBuilder((N_MAX_X-2)*(N_MAX_Y-2)); //Creates equations from "A", "B", and "HHP" and applies Newton's method until convergence.
	printf("5\n");

	_fcloseall();
}

//There are two possible coordinate systems. (x,y) reflects a grid, and (s) is just all the grid points in a line.
//Converts between (x,y) and (s)
__forceinline int RN(int m, int n)
{
	return n*(N_MAX_X-2)+m;
}

//Converts between (s) and (x,y)
__forceinline int CN_X(int r)
{
	return r%(N_MAX_X-2);
}

__forceinline int CN_Y(int r)
{
	return r/(N_MAX_X-2);
}

void writeToFiles(void)
{		
	rewind(fp_r);
	rewind(fp_z);
	rewind(fp_a);
	rewind(fp_b);
	rewind(fp_j);
	rewind(fp_j2);
	rewind(fp_hhp);
	rewind(fp_hhpPRIME);
	rewind(fp_v);
	rewind(fp_vxy);
	rewind(fp_P_HHP);
	rewind(fp_change);

	//Write actual R values.
	for (int j=0;j<=N_MAX_X-2;j++)
	{
		fprintf(fp_r, "%f\n",j*DX);   
	}    
	
	//Write actual Z values.
	for (int k=N_MAX_Y-2;k>=0;k--)
	{
		fprintf(fp_z, "%f\n",k*DY);   
	}  
		
	//Write "A" matrix .
	for (int j=0;j<=(N_MAX_X-2)*(N_MAX_Y-2)-1;j++)
	{
		for (int k=0;k<=(N_MAX_X-2)*(N_MAX_Y-2)-1;k++)
		{
			fprintf(fp_a, "%40.20f ",A[j][k]);
		}
		fprintf(fp_a, "\n");
	}

	//Write "b" vector.
	for (int k=0;k<(N_MAX_X-2)*(N_MAX_Y-2);k++)
	{
		fprintf(fp_b, "%40.20f\n",B[k]);   
	}

	//Write "J" matrix.
	for (int j=0;j<=(N_MAX_X-2)*(N_MAX_Y-2)-1;j++)
	{
		for (int k=0;k<=(N_MAX_X-2)*(N_MAX_Y-2)-1;k++)
		{
 		           fprintf(fp_j, "%40.20f ",J[j][k]);
        		}
        		fprintf(fp_j, "\n");
    	}

#ifndef READ
	//Write "J-1" matrix .
	for (int j=0;j<=(N_MAX_X-2)*(N_MAX_Y-2)-1;j++)
	{
		for (int k=0;k<=(N_MAX_X-2)*(N_MAX_Y-2)-1;k++)
		{
 		           fprintf(fp_j2, "%40.20f ",J_INVERSE[j][k]);
        		}
        		fprintf(fp_j2, "\n");
    	}   	
#endif
	
	//Write HPP.
	for (int k=N_MAX_Y-3;k>=0;k--)
	{
		for (int j=0;j<N_MAX_X-2;j++)
		{
			fprintf(fp_hhp, "%40.20f ",HHP[j][k]);
        		}
        		fprintf(fp_hhp, "\n");
    	}
	
	//Write D(HHP)/D(PHI).
	for (int k=N_MAX_Y-3;k>=0;k--)
	{
		for (int j=0;j<N_MAX_X-2;j++)
		{			
			fprintf(fp_hhpPRIME, "%40.20f ",HHP_PRIME[j][k]);
        		}        		
		fprintf(fp_hhpPRIME, "\n");
    	}

	//Write "v" vector.
	fprintf(fp_v, "%40.20f ",0.0);	
	for (int j=0;j<N_MAX_X-2;j++)
	{
		//Deal with bottom boundary.
#if TYPE==MONOPOLE
		//Monopole P=P_OP
		fprintf(fp_v, "%40.20f ",P_OP);
#else
		//Others
		if(j<N_LC+1)
		{
			//All non-monopole cases Pz=0
			fprintf(fp_v, "%40.20f ",v0[RN(j,0)]);
		}
		else
		{
#if (TYPE==STANDARD||TYPE==JETS)
			//CKF, Takamori, Jets P=P_OP
			fprintf(fp_v, "%40.20f ",P_OP);
#elif TYPE==NULL
			//CKF_NULL, TAK_NULL H^2=(R^2-1)*(Pz)^2
			fprintf(fp_v, "%40.20f ",v0[RN(j,0)]-L_func[j]);
#endif
		}
		
#endif
	}	
	for (int k=0;k<N_MAX_Y-2;k++)
	{
		fprintf(fp_v, "%40.20f ",0.0);
		for (int j=0;j<N_MAX_X-2;j++)
		{
			fprintf(fp_v, "%40.20f ",v0[RN(j,k)]);
        		}	
    	}
	
	//Write "v" grid.
	for (int k=N_MAX_Y-3;k>=0;k--)
	{
		for (int j=0;j<N_MAX_X-2;j++)
		{
			fprintf(fp_vxy, "%40.20f ",v0[RN(j,k)]);			
        		}
        		fprintf(fp_vxy, "\n");		
    	}	
		
	//Write P vs HHP list.
	for(int k=0;k<=N_MAX_Y-3;k++)
	{		
		for(int j=0;j<=N_MAX_X-3;j++)
		{
			if(j>N_S_X||k>N_S_Y)
			{			
				fprintf(fp_P_HHP,"%40.20f %40.20f\n",v0[RN(j,k)],HHP[j][k]);
			}
		}
	}	
	
	//Write "change of PSI" grid.
	for (int k=N_MAX_Y-3;k>=0;k--)
	{
		for (int j=0;j<N_MAX_X-2;j++)
		{
			fprintf(fp_change, "%40.20f ",CHANGE[j][k]);
        		}
        		fprintf(fp_change, "\n");		
    	}	

	fflush(fp_r);
	fflush(fp_z);
	fflush(fp_a);
	fflush(fp_b);
	fflush(fp_j);
	fflush(fp_j2);
	fflush(fp_hhp);
	fflush(fp_hhpPRIME);
	fflush(fp_v);
	fflush(fp_vxy);
	fflush(fp_P_HHP);
	fflush(fp_change);
}

//Function to allocate memory for a 2D array of elements of type T.
T** Make2DTArray(int arraySizeX, int arraySizeY)
{  
	T** theArray;
	theArray = (T**) malloc(arraySizeX*sizeof(T*));  
	for (int i=0;i<arraySizeX;i++)
	{
		theArray[i] = (T*) malloc(arraySizeY*sizeof(T));
	}
	return theArray;  
}   
