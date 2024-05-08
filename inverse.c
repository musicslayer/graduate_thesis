//Returns the inverse of the matrix stored in 2D array "m"
void inverse(T** m,int size,T** g)
{	
	//Creates a copy of the input matrix so that the original is not destroyed
	T **f;
	f=Make2DTArray((N_MAX_X-2)*(N_MAX_Y-2),(N_MAX_X-2)*(N_MAX_Y-2));
	printf("MBytes f=%f\n",(T)SIZE_T_XY_XY/(1024*1024));

	//The memcpy calls won't work if "f" is dynamically allocated, so we have to manually copy "f" below.
	//memcpy(f,m,(N_MAX_X-2)*(N_MAX_Y-2)*(N_MAX_X-2)*(N_MAX_Y-2));

	//Initializes g as the identity matrix, which will turn into the inverse of the original matrix
	for(int x=0;x<size;x++)
	{   
		for(int y=0;y<size;y++)
		{
			g[x][y]=(x==y);
			f[x][y]=m[x][y];
		}
	}
	for(int j=0;j<size;j++)
	{    			
		printf("j=%i out of %i\n",j+1,size);
		/*
		//This code decides when to swap rows, and if the matrix is singular.
		//The memcpy calls won't work if "f" is dynamically allocated.
		//This needs to be rewritten
		if(fabs(f[j][j])<0.0000001)		
		{
			//const char key=_getch();
			//exit(0);
			int isSingular=1;
			//Swap rows to avoid a zero pivot
			for(int s=j;s<size;s++)
			{
				if(fabs(f[s][j])>0.0000001)
				{					
					T temp[(N_MAX_X-2)*(N_MAX_Y-2)];

					isSingular=0;
					memcpy(temp,f+j,SIZE_T_XY);
					memcpy(f+j,f+s,SIZE_T_XY);
					memcpy(f+s,temp,SIZE_T_XY);		
					memcpy(temp,g+j,SIZE_T_XY);
					memcpy(g+j,g+s,SIZE_T_XY);
					memcpy(g+s,temp,SIZE_T_XY);				
					break;
				}
			}
			//if(isSingular)
			{				
				//printf("YOUR JACOBIAN IS SINGULAR!!!!! CANNOT GO ON!!!!!");
				//printf("DET_J=%f\n",determinantAlt(J,(N_MAX-2)*(N_MAX-2)));
				//printf("DET_A=%f\n",DET_A);
				//exit(0);
				//return;
			}
			//multAccum*=-f[j][j];
		}
		else
		{
			//multAccum*=f[j][j];
		}
		*/
		T temp=1/f[j][j];		
		for(int l=size-1;l>j;l--)
		{
			f[j][l]*=temp;
		}
		for(int l=j;l>=0;l--)
		{	
			g[j][l]*=temp;
		}
		
		for(int k1=0;k1<j;k1++)
		{
			for(int k2=size-1;k2>j;k2--)
			{
				f[k1][k2]-=f[k1][j]*f[j][k2];		
			}	
			for(int k2=j;k2>=0;k2--)
			{
				g[k1][k2]-=f[k1][j]*g[j][k2];
			}			
		}
		for(int k1=j+1;k1<size;k1++)
		{
			for(int k2=size-1;k2>j;k2--)
			{
				f[k1][k2]-=f[k1][j]*f[j][k2];		
			}				
			for(int k2=j;k2>=0;k2--)
			{
				g[k1][k2]-=f[k1][j]*g[j][k2];			
			}				
		}
	}
}
