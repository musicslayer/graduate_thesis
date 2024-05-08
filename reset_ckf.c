void InsertionSort(void)
{
 	int i,j;
	T keyX,keyY;
	for(j=1;j<(N_MAX_X-2)*(N_MAX_Y-2);j++)// Start with 1 (not 0)
	{
		keyX=L_X[j];
		keyY=L_Y[j];
		for(i=j-1;(i >= 0)&&(L_X[i]<keyX);i--)// Smaller values move up
		{
			L_X[i+1] = L_X[i];
			L_Y[i+1] = L_Y[i];
		}
		L_X[i+1] = keyX;
		L_Y[i+1] = keyY;
	}
	return;
}

void resetCKF(void)
{
	int size=(N_MAX_X-2)*(N_MAX_Y-2);
	//We need to get PSI vs. HHP so we can take integral via trapezoidal method.
	for(int i=0;i<size;i++)
	{
		if(CN_Y(i)>0)
		{
			L_X[i]=v0[i];
			L_Y[i]=fabs(HHP[CN_X(i)][CN_Y(i)]); //Do we really need to take the absolute value.
		}
	}
	InsertionSort(); //This makes L_X sorted in DESCENDING order (L_X[0] is largest, L_X[size-1] is smallest).
	
	//Bottom Edge (past light cylinder).
	for(int s=N_LC;s<=N_MAX_X-4;s++)	
	{
		T sum=0;
		if(v0[s]<P_OP)
		{
			for(int i=size-1;i>0;i--)
			{
				if(L_X[i]>v0[s]){break;}			
				sum+=0.5*(L_X[i-1]-L_X[i])*(L_Y[i]+L_Y[i-1]);
			}
			sum=fabs(sum);//This ensures that the term under the square root is positive. We shouldn't really need this...
		}

		//Note: H is negative, so we want -sqrt(H^2).
		bB1[s]=0;
		bB2[s]=0;
		bB3[s]=-1;
		bB4[s]=0;
		bB5[s]=1;
		bB6[s]=-DY*sqrt(2*sum/((s+1)*(s+1)*DX*DX-1));
		
		L_func[s]=-DY*sqrt(2*sum/((s+1)*(s+1)*DX*DX-1)); //This is used just so that we can add the bottom boundary to plots.
	}

	//Bottom Right Corner.
	int s=N_MAX_X-3;
	T sum=0;
	if(v0[s]<P_OP)
	{
		for(int i=size-1;i>0;i--)
		{
			if(L_X[i]>v0[s]){break;}
			sum+=0.5*(L_X[i-1]-L_X[i])*(L_Y[i]+L_Y[i-1]);
		}
	}
	
	int t=0;
	int i=s+1;
	int j=t+1;
	
	//Note: H is negative, so we want -sqrt(H^2).
	bR1[t]=-i;
	bR2[t]=i;
	bR3[t]=0;
	bR4[t]=j;
	bR5[t]=-j;
	bR6[t]=DY*j*sqrt(2*sum/(i*i*DX*DX-1));

	bB1[s]=0;
	bB2[s]=0;
	bB3[s]=-1;
	bB4[s]=0;
	bB5[s]=1;
	bB6[s]=-DY*sqrt(2*sum/(i*i*DX*DX-1));
	
	L_func[s]=-DY*sqrt(2*sum/((s+1)*(s+1)*DX*DX-1));

	create_matrix();

	for(int x=0;x<size;x++)
	{
		for(int y=0;y<size;y++)
		{		
			J[x][y]=A[x][y];
			h[x][y]=A[x][y];
		}
	}
}
