__forceinline void create_matrix(void)
{
	//"Interior interior" points (points that do not touch a boundary)
	for(int s=1;s<=N_MAX_X-4;s++)
	{		
		for(int t=1;t<=N_MAX_Y-4;t++)
		{		
			A[RN(s,t)][RN(s-1,t)]=a1[s][t];
			A[RN(s,t)][RN(s+1,t)]=a2[s][t];
			A[RN(s,t)][RN(s,t-1)]=a3[s][t];
			A[RN(s,t)][RN(s,t+1)]=a4[s][t];
			A[RN(s,t)][RN(s,t)]=a5[s][t];				
			B[RN(s,t)]=a6[s][t]-DX*DX*f(s+1,t+1);		
		}		
	}

	//Left Edge
	for(int t=1;t<=N_MAX_Y-4;t++)
	{
		int s=0;

		//A[RN(s,t)][RN(s-1,t)]=0;
		A[RN(s,t)][RN(s+1,t)]=a2[s][t]-a1[s][t]*bL2[t]/bL1[t];
		A[RN(s,t)][RN(s,t-1)]=a3[s][t]-a1[s][t]*bL3[t]/bL1[t];
		A[RN(s,t)][RN(s,t+1)]=a4[s][t]-a1[s][t]*bL4[t]/bL1[t];
		A[RN(s,t)][RN(s,t)]=a5[s][t]-a1[s][t]*bL5[t]/bL1[t];
		B[RN(s,t)]=a6[s][t]-a1[s][t]*bL6[t]/bL1[t]-DX*DX*f(s+1,t+1);		
	}

	//Right Edge
	for(int t=1;t<=N_MAX_Y-4;t++)
	{
		int s=N_MAX_X-3;

		A[RN(s,t)][RN(s-1,t)]=a1[s][t]-a2[s][t]*bR1[t]/bR2[t];
		//A[RN(s,t)][RN(s+1,t)]=0;
		A[RN(s,t)][RN(s,t-1)]=a3[s][t]-a2[s][t]*bR3[t]/bR2[t];
		A[RN(s,t)][RN(s,t+1)]=a4[s][t]-a2[s][t]*bR4[t]/bR2[t];
		A[RN(s,t)][RN(s,t)]=a5[s][t]-a2[s][t]*bR5[t]/bR2[t];
		B[RN(s,t)]=a6[s][t]-a2[s][t]*bR6[t]/bR2[t]-DX*DX*f(s+1,t+1);		
	}

	//Bottom Edge
	for(int s=1;s<=N_MAX_X-4;s++)
	{
		int t=0;

		A[RN(s,t)][RN(s-1,t)]=a1[s][t]-a3[s][t]*bB1[s]/bB3[s];
		A[RN(s,t)][RN(s+1,t)]=a2[s][t]-a3[s][t]*bB2[s]/bB3[s];
		//A[RN(s,t)][RN(s,t-1)]=0;
		A[RN(s,t)][RN(s,t+1)]=a4[s][t]-a3[s][t]*bB4[s]/bB3[s];
		A[RN(s,t)][RN(s,t)]=a5[s][t]-a3[s][t]*bB5[s]/bB3[s];		
		B[RN(s,t)]=a6[s][t]-a3[s][t]*bB6[s]/bB3[s]-DX*DX*f(s+1,t+1);		
	}

	//Top Edge
	for(int s=1;s<=N_MAX_X-4;s++)
	{
		int t=N_MAX_Y-3;

		A[RN(s,t)][RN(s-1,t)]=a1[s][t]-a4[s][t]*bT1[s]/bT4[s];
		A[RN(s,t)][RN(s+1,t)]=a2[s][t]-a4[s][t]*bT2[s]/bT4[s];
		A[RN(s,t)][RN(s,t-1)]=a3[s][t]-a4[s][t]*bT3[s]/bT4[s];
		//A[RN(s,t)][RN(s,t+1)]=0;
		A[RN(s,t)][RN(s,t)]=a5[s][t]-a4[s][t]*bT5[s]/bT4[s];
		B[RN(s,t)]=a6[s][t]-a4[s][t]*bT6[s]/bT4[s]-DX*DX*f(s+1,t+1);
	}	

	//Bottom Left Corner
	{
		int s=0;
		int t=0;

		//A[RN(s,t)][RN(s-1,t)]=0;
		A[RN(s,t)][RN(s+1,t)]=a2[s][t]-a1[s][t]*bL2[t]/bL1[t]-a3[s][t]*bB2[s]/bB3[s];
		//A[RN(s,t)][RN(s,t-1)]=0;
		A[RN(s,t)][RN(s,t+1)]=a4[s][t]-a1[s][t]*bL4[t]/bL1[t]-a3[s][t]*bB4[s]/bB3[s];
		A[RN(s,t)][RN(s,t)]=a5[s][t]-a1[s][t]*bL5[t]/bL1[t]-a3[s][t]*bB5[s]/bB3[s];
		B[RN(s,t)]=a6[s][t]-a1[s][t]*bL6[t]/bL1[t]-a3[s][t]*bB6[s]/bB3[s]-DX*DX*f(s+1,t+1);
	}

	//Top Left Corner
	{
		int s=0;
		int t=N_MAX_Y-3;

		//A[RN(s,t)][RN(s-1,t)]=0;
		A[RN(s,t)][RN(s+1,t)]=a2[s][t]-a1[s][t]*bL2[t]/bL1[t]-a4[s][t]*bT2[s]/bT4[s];
		A[RN(s,t)][RN(s,t-1)]=a3[s][t]-a1[s][t]*bL3[t]/bL1[t]-a4[s][t]*bT3[s]/bT4[s];
		//A[RN(s,t)][RN(s,t+1)]=0;
		A[RN(s,t)][RN(s,t)]=a5[s][t]-a1[s][t]*bL5[t]/bL1[t]-a4[s][t]*bT5[s]/bT4[s];
		B[RN(s,t)]=a6[s][t]-a1[s][t]*bL6[t]/bL1[t]-a4[s][t]*bT6[s]/bT4[s]-DX*DX*f(s+1,t+1);
	}

	//Bottom Right Corner
	{
		int s=N_MAX_X-3;
		int t=0;

		A[RN(s,t)][RN(s-1,t)]=a1[s][t]-a2[s][t]*bR1[t]/bR2[t]-a3[s][t]*bB1[s]/bB3[s];
		//A[RN(s,t)][RN(s+1,t)]=0;
		//A[RN(s,t)][RN(s,t-1)]=0;
		A[RN(s,t)][RN(s,t+1)]=a4[s][t]-a2[s][t]*bR4[t]/bR2[t]-a3[s][t]*bB4[s]/bB3[s];
		A[RN(s,t)][RN(s,t)]=a5[s][t]-a2[s][t]*bR5[t]/bR2[t]-a3[s][t]*bB5[s]/bB3[s];
		B[RN(s,t)]=a6[s][t]-a2[s][t]*bR6[t]/bR2[t]-a3[s][t]*bB6[s]/bB3[s]-DX*DX*f(s+1,t+1);
	}	

	//Top Right Corner
	{
		int s=N_MAX_X-3;
		int t=N_MAX_Y-3;

		A[RN(s,t)][RN(s-1,t)]=a1[s][t]-a2[s][t]*bR1[t]/bR2[t]-a4[s][t]*bT1[s]/bT4[s];
		//A[RN(s,t)][RN(s+1,t)]=0;
		A[RN(s,t)][RN(s,t-1)]=a3[s][t]-a2[s][t]*bR3[t]/bR2[t]-a4[s][t]*bT3[s]/bT4[s];
		//A[RN(s,t)][RN(s,t+1)]=0;
		A[RN(s,t)][RN(s,t)]=a5[s][t]-a2[s][t]*bR5[t]/bR2[t]-a4[s][t]*bT5[s]/bT4[s];
		B[RN(s,t)]=a6[s][t]-a2[s][t]*bR6[t]/bR2[t]-a4[s][t]*bT6[s]/bT4[s]-DX*DX*f(s+1,t+1);
	}
}
