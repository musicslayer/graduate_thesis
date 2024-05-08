/*
Refers to the type of simulation. 
If it is a double simulation, this must match whatever the end graph will be of.
If there is more than one type of simulation mixed together, this must match whatever the bottom boundary is of.
*/
#define TYPE MONOPOLE
#define i ((T)(s+1))
#define j ((T)(t+1))

void initialize(void)
{
	//Equation coefficients
	for(int s=0;s<=N_MAX_X-3;s++)
	{		
		for(int t=0;t<=N_MAX_Y-3;t++)
		{
			/*
			//These equations represent an alternate finite difference choice for the first derivative.
			//The upper right corner is impossible to solve for using this choice, so I don't recommend this option.
			//1st derivative is central difference
			a1[s][t]=1-i*i*DX*DX+0.5*(1/i+i*DX*DX);
			a2[s][t]=1-i*i*DX*DX-0.5*(1/i+i*DX*DX);
			a3[s][t]=1-i*i*DX*DX;
			a4[s][t]=1-i*i*DX*DX;
			a5[s][t]=-4+4*i*i*DX*DX;
			a6[s][t]=0;
			*/

			/*
			//1st derivative is backward difference
			a1[s][t]=1-i*i*DX*DX+1/i+i*DX*DX;
			a2[s][t]=1-i*i*DX*DX;
			a3[s][t]=1-i*i*DX*DX;
			a4[s][t]=1-i*i*DX*DX;
			a5[s][t]=-4+4*i*i*DX*DX-1/i-i*DX*DX;
			a6[s][t]=0;						
			*/
			
			//1st derivative is backward difference, DX and DY independent
			a1[s][t]=1-i*i*DX*DX+1/i+i*DX*DX;
			a2[s][t]=1-i*i*DX*DX;
			a3[s][t]=(DX*DX/(DY*DY))*(1-i*i*DX*DX);
			a4[s][t]=(DX*DX/(DY*DY))*(1-i*i*DX*DX);
			a5[s][t]=(-2-2*(DX*DX/(DY*DY)))*(1-i*i*DX*DX)-1/i-i*DX*DX;
			a6[s][t]=0;			
			
			//Value near LC is average of values to the left and right.
			if(s>=N_LC-1-SMOOTH&&s<=N_LC-1+SMOOTH)
			{		
				a1[s][t]=-0.5;
				a2[s][t]=-0.5;
				a3[s][t]=0;
				a4[s][t]=0;
				a5[s][t]=1;
				a6[s][t]=0;
			}	
			
			/*
			//Unused
			//Polynomial coefficients of HHP (constant term taken care of in function "f")
			//c[#] is the coefficient for (Pij)^(#+1)
			c[0][s][t]=4;
			c[1][s][t]=-6;
			c[2][s][t]=2;
			*/
		}
	}

	//Left Boundary P=0
	for(int t=1;t<=N_MAX_Y-4;t++)
	{
		int s=0;

		bL1[t]=1;
		bL2[t]=0;
		bL3[t]=0;
		bL4[t]=0;
		bL5[t]=0;
		bL6[t]=0;
	}

	//Right Boundary RPr+ZPz=0
	for(int t=1;t<=N_MAX_Y-4;t++)
	{
		int s=N_MAX_X-3;

		bR1[t]=-i;
		bR2[t]=i;
		bR3[t]=-j;
		bR4[t]=j;
		bR5[t]=0;
		bR6[t]=0;
	}

	//Bottom Boundary P=P_OP
	for(int s=1;s<=N_MAX_X-4;s++)
	{
		int t=0;

		bB1[s]=0;
		bB2[s]=0;
		bB3[s]=1;
		bB4[s]=0;
		bB5[s]=0;
		bB6[s]=P_OP;
	}

	//Top Boundary RPr+ZPz=0
	for(int s=1;s<=N_MAX_X-4;s++)
	{
		int t=N_MAX_Y-3;

		bT1[s]=-i;
		bT2[s]=i;
		bT3[s]=-j;
		bT4[s]=j;
		bT5[s]=0;
		bT6[s]=0;
	}

	//Bottom Left Corner
	{
		int s=0;
		int t=0;

		bL1[t]=1;
		bL2[t]=0;
		bL3[t]=0;
		bL4[t]=0;
		bL5[t]=0;
		bL6[t]=0;

		bB1[s]=0;
		bB2[s]=0;
		bB3[s]=1;
		bB4[s]=0;
		bB5[s]=0;
		bB6[s]=P_OP;
	}

	//Top Left Corner
	{
		int s=0;
		int t=N_MAX_Y-3;

		bL1[t]=1;
		bL2[t]=0;
		bL3[t]=0;
		bL4[t]=0;
		bL5[t]=0;
		bL6[t]=0;

		bT1[s]=0;
		bT2[s]=i;
		bT3[s]=-j;
		bT4[s]=j;
		bT5[s]=0;		
		bT6[s]=0;
	}

	//Bottom Right Corner
	{
		int s=N_MAX_X-3;
		int t=0;

		bR1[t]=-i;
		bR2[t]=i;
		bR3[t]=0;
		bR4[t]=j;
		bR5[t]=0;
		bR6[t]=j*P_OP;

		bB1[s]=0;
		bB2[s]=0;
		bB3[s]=1;
		bB4[s]=0;
		bB5[s]=0;
		bB6[s]=P_OP;
	}

	//Top Right Corner
	{
		int s=N_MAX_X-3;
		int t=N_MAX_Y-3;

		bR1[t]=-1;
		bR2[t]=1;
		bR3[t]=0;
		bR4[t]=0;
		bR5[t]=0;
		bR6[t]=0;

		bT1[s]=0;
		bT2[s]=0;
		bT3[s]=-1;
		bT4[s]=1;
		bT5[s]=0;
		bT6[s]=0;
	}	
}

//This function inserts the star.
void star(void){}

//This fills in HHP and d(HHP)/d(PSI), which is a function of "PSI" and may or may not include both linear and nonlinear terms.
__forceinline void hhpSet(T** HHP,T** HHP_PRIME)
{
	T HHP_L[N_MAX_Y-2];
	for(int KP=0;KP<N_MAX_Y-2;KP++)
	{		
		HHP_L[KP]=(v0[RN(N_LC-2,KP)]-v0[RN(N_LC-3,KP)]+v0[RN(N_LC+1,KP)]-v0[RN(N_LC,KP)])*(1/DX);
	}		
	
	//At r<rL, there are many possible cases (closed field lines, open field lines that do or do not cross the light cylinder).
	//For 0<r<=NR_S, we stay OUTSIDE star, because the star value is known and does not need to be altered.
	
	for(int jj=0;jj<=N_LC-2;jj++)
	{	
		for(int kk=0;kk<=N_MAX_Y-3;kk++)
		{
			int KP=0;
			{	
				const T PT=v0[RN(jj,kk)];			
							
				if(PT<v0[RN(N_LC-2,0)])
				{	
					if(PT<v0[RN(N_LC-2,N_MAX_Y-3)])
					{
						HHP[jj][kk]=HHP_L[N_MAX_Y-3]*PT/v0[RN(N_LC-2,N_MAX_Y-3)];
					}						
					else
					{	
						for(;KP<N_MAX_Y-4&&PT<v0[RN(N_LC-2,KP)];KP++);
						
						const T Q1=PT-v0[RN(N_LC-1,KP+1)];
						const T Q2=PT-v0[RN(N_LC-1,KP)];								
						HHP[jj][kk]=(Q1*HHP_L[KP]-Q2*HHP_L[KP+1])/(Q1-Q2);
					}
				}
				else
				{
					HHP[jj][kk]=0;
				}		
			}
			if(isnan(HHP[jj][kk])){HHP[jj][kk]=0;}
		}
	}

	for(int jj=N_LC-1;jj==N_LC-1;jj++)
	{				
		for(int kk=0;kk<=N_MAX_Y-3;kk++)
		{
			HHP[N_LC-1][kk]=HHP_L[kk];									
		}		
	}	

	for(int jj=N_LC;jj<=N_MAX_X-3;jj++)
	{				
		for(int kk=0;kk<=N_MAX_Y-3;kk++)
		{
			int KP=0;
			const T PT=v0[RN(jj,kk)];

			for(;KP<N_MAX_Y-4&&PT<v0[RN(N_LC-2,KP)];KP++);
			
			const T Q1=PT-v0[RN(N_LC-1,KP+1)];
			const T Q2=PT-v0[RN(N_LC-1,KP)];
			HHP[jj][kk]=(Q1*HHP_L[KP]-Q2*HHP_L[KP+1])/(Q1-Q2);			
			if(isnan(HHP[jj][kk])){HHP[jj][kk]=0;}
		}		
	}			

	for(int jj=N_LC-1;jj==N_LC-1;jj++)
	{				
		for(int kk=0;kk<=N_MAX_Y-3;kk++)
		{
			v0[RN(N_LC-1,kk)]=0.5*(v0[RN(N_LC-2,kk)]+v0[RN(N_LC,kk)]);		
		}		
	}	

	//TODO - HHP_PRIME is largely unnecessary, and this calculation may be wrong.
	for(int x=0;x<N_MAX_X-2;x++)
	{		
		for(int y=0;y<N_MAX_Y-2;y++)
		{
			HHP_PRIME[x][y]=0.1*sqrt(x*x+y*y);			
		}
	}	
}

//This represents the part of HHP that is a function of "R" and "Z" (NOT "PSI"). This includes any possible constant term.
__forceinline T f(int m, int n)
{
	return 0;
}

#undef i
#undef j
