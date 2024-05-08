//Solves for PHI using Newton's method.
//Note: PHI is stored in v0 and v1.
void equationBuilder(int size)
{	
	T vINIT[(N_MAX_X-2)*(N_MAX_Y-2)];
	
#ifndef READ
	for(int x=0;x<size;x++)
	{
		for(int y=0;y<size;y++)
		{		
			J[x][y]=A[x][y];
			h[x][y]=A[x][y];
		}
	}
#endif

#ifdef SHORTCUT
#ifndef READ
	printf("pre inverse\n");
	inverse(J,size,J_INVERSE);
	printf("after inverse\n");
#endif
#endif

	//Our initial guess for every point is one.
	for(int i=0;i<size;i++)
	{		
		vINIT[i]=1.0;
		v0[i]=vINIT[i];
		v1[i]=vINIT[i];		
	}

	/*
	//An alternate way of setting initial values.
	//Our initial guess for every point is "random" between 0 and 1.
	srand(time(NULL)); 
	for(int i=0;i<size;i++)
	{				
		vINIT[i]=((T)rand())/RAND_MAX;
		v0[i]=vINIT[i];
		v1[i]=vINIT[i];		
	}
	*/

	T olddeltaALL=0;
	T deltaALL=0;
	for(int iiii=1;iiii<=NUM_ITERATIONS_MAX;iiii++)
	{		
		olddeltaALL=deltaALL;
		deltaALL=0;
		T deltaMAX=0;

		//This code controls user interaction. All user input is direct keyboard hits, NOT typing text into a command line.
		//To the best of my knowledge, this code to check for user input does NOT slow the program down significantly.
		if(_kbhit())
		{
			//This call to _getch() will not pause the program because the user already hit a key,
			//but in general the other calls will pause the program and wait for input.
			const char key=_getch();
			if (key=='b'||key=='B')
			{
				//Simulation is paused.
				printf("\nWriting data to files............\n");
				writeToFiles();
				printf("\nData has been written to files.\n\nEnter 'w' or 'W' to change the FRAC parameter.\nEnter 'e' or 'E' to change BETA.\nEnter 'r' or 'R' to reset the same simulation.\nEnter 'c' or 'C' to reset with another BETA value.\nEnter 'q' or 'Q' to quit.\nEnter anything else to continue.\n");				

				const char key=_getch();
				if (key=='w'||key=='W')
				{
					//Yes, the code needs to be structured this way!
					while(1)
					{
						printf("Enter new V_FRAC value:\n");
						scanf("%lf",&V_FRAC);
						fflush(stdin);//NOT STANDARD :(
						rewind(stdin);//NOT STANDARD :(
						printf("\nV_FRAC=%lf - Is this OK?\nEnter 'y' or 'Y' to continue\nEnter 'a' or 'A' to try again.\n\n",V_FRAC);

						char key;
						while(1)
						{
							key=_getch();
							if (key=='y'||key=='Y'||key=='a'||key=='A')
							{
								break;
							}							
						}	
						if (key=='y'||key=='Y')
						{
							break;
						}	
					}
				}
				else if (key=='e'||key=='E')
				{
					while(1)
					{
						printf("Enter new BETA value:\n");
						scanf("%lf",&BETA);
						fflush(stdin);//NOT STANDARD :(
						rewind(stdin);//NOT STANDARD :(
						printf("\nBETA=%lf - Is this OK?\nEnter 'y' or 'Y' to continue\nEnter 'a' or 'A' to try again.\n\n",BETA);

						char key;
						while(1)
						{
							key=_getch();
							if (key=='y'||key=='Y'||key=='a'||key=='A')
							{
								break;
							}							
						}	
						if (key=='y'||key=='Y')
						{
							break;
						}	
					}
				}
				else if (key=='r'||key=='R')
				{
					//Completely reset simulation.
					for(int i=0;i<size;i++)
					{		
						vINIT[i]=1;
						v0[i]=vINIT[i];
						v1[i]=vINIT[i];
						iiii=1;
					}
				}
				else if (key=='c'||key=='C')
				{
					//Completely reset simulation with different BETA value.
					for(int i=0;i<size;i++)
					{		
						vINIT[i]=1;
						v0[i]=vINIT[i];
						v1[i]=vINIT[i];
						iiii=1;
					}

					while(1)
					{
						printf("Enter new BETA value:\n");
						scanf("%lf",&BETA);
						fflush(stdin);//NOT STANDARD :(
						rewind(stdin);//NOT STANDARD :(
						printf("\nBETA=%lf - Is this OK?\nEnter 'y' or 'Y' to continue\nEnter 'a' or 'A' to try again.\n\n",BETA);

						char key;
						while(1)
						{
							key=_getch();
							if (key=='y'||key=='Y'||key=='a'||key=='A')
							{
								break;
							}							
						}	
						if (key=='y'||key=='Y')
						{
							break;
						}	
					}
				}		
				else if (key=='q'||key=='Q')
				{
					//Quits the program.
					printf("**************************************************\n");
					printf("%i out of %i iterations were completed.\n",iiii,NUM_ITERATIONS_MAX);
					_fcloseall();					
					exit(0);
				}
				else
				{
					printf("\n\nCONTINUE\n\n");
				}
			}
		}

		printf("\n         ITERATION=%i    V_FRAC=%f      BETA=%f\n",iiii,V_FRAC,BETA);

#ifndef SHORTCUT
		//Updates the Jacobian for each iteration of Newton's Method.
		for(int x=0;x<size;x++)
		{
			p[x]=DX*DX*HHP_PRIME[CN_X(x)][CN_Y(x)];
			J[x][x]=h[x][x]+p[x];
		}
			
		inverse(J,size,J_INVERSE);
#endif		
		//Calculates the nonlinear equations that we are trying to solve to be zero.
		for(int i=0;i<size;i++)
		{
			eq[i]=0;
			for(int k=0;k<size;k++)
			{
				eq[i]+=A[i][k]*v0[k];
			}
			eq[i]+=-B[i]+HHP[CN_X(i)][CN_Y(i)]*DX*DX;
		}		

		int s_max=0;
		//Newton's method.
		for(int i=0;i<size;i++)
		{		
			//We could do better by not including the star in our solver, since the values are already known.
			//if((CN_X(i)>(N_S_X-1))||(CN_Y(i)>(N_S_Y-1)))
			{
				//Update all variables.
				v1[i]=0;
				for(int j=0;j<size;j++)
				{
					v1[i]+=J_INVERSE[i][j]*eq[j];
				}				
				v1[i]=v0[i]-v1[i];

				//Keeps track of changes between iterations.
				deltaALL+=fabs(v0[i]-v1[i]);
				CHANGE[CN_X(i)][CN_Y(i)]=fabs(v0[i]-v1[i]);
				if((CN_X(i)>(N_S_X-1))||(CN_Y(i)>(N_S_Y-1)))
				{
					if (fabs(v0[i]-v1[i])>deltaMAX)
					{
						s_max=i;
						deltaMAX=fabs(v0[i]-v1[i]);
					}
				}
			}
		}
		printf("Max change at %i, or (%i,%i)\n",s_max,CN_X(s_max),CN_Y(s_max));
		printf("deltaALL=%f\n",deltaALL);
		printf("olddeltaALL=%f\n",olddeltaALL);
		printf("deltaMAX=%f\n",deltaMAX);
		
		//Similar effect to the relaxation parameter.
		for(int i=0;i<size;i++)
		{
			v0[i]=(V_FRAC)*(v1[i])+(1-V_FRAC)*(v0[i]);
		}		

		//Automatically end if a certain convergence level is reached.
		//User is provided with options as to what to do next.
		if(iiii>=5&&deltaMAX<0.000001)
		{
			printf("STABLE BREAK\n");
			printf("\nWriting data to files............\n");
			writeToFiles();
			
			printf("\nData has been written to files.\n\nEnter 'r' or 'R' to reset the same simulation.\nEnter 'c' or 'C' to reset with another BETA value.\nEnter 'e' or 'E' to continue with another BETA value.\nEnter anything else to quit.\n");
			const char key=_getch();

			if (key=='r'||key=='R')
			{
				//Completely reset simulation.
				for(int i=0;i<size;i++)
				{		
					vINIT[i]=1;
					v0[i]=vINIT[i];
					v1[i]=vINIT[i];					
				}
				iiii=0;
			}
			else if (key=='c'||key=='C')
			{
				//Completely reset simulation with different BETA value.
				for(int i=0;i<size;i++)
				{		
					vINIT[i]=1;
					v0[i]=vINIT[i];
					v1[i]=vINIT[i];					
				}
				iiii=0;

				while(1)
				{
					printf("Enter new BETA value:\n");
					scanf("%lf",&BETA);
					fflush(stdin);//NOT STANDARD :(
					rewind(stdin);//NOT STANDARD :(
					printf("\nBETA=%lf - Is this OK?\nEnter 'y' or 'Y' to continue\nEnter 'a' or 'A' to try again.\n\n",BETA);

					char key;
					while(1)
					{
						key=_getch();
						if (key=='y'||key=='Y'||key=='a'||key=='A')
						{
							break;
						}							
					}	
					if (key=='y'||key=='Y')
					{
						break;
					}	
				}
			}
			else if (key=='e'||key=='E')
			{
				//Continue simulation with different BETA value.
				iiii=0; //This is reset here so that the maximum iteration limit isn't reached prematurely if you switch BETA a lot.
				while(1)
				{
					printf("Enter new BETA value:\n");
					scanf("%lf",&BETA);
					fflush(stdin);//NOT STANDARD :(
					rewind(stdin);//NOT STANDARD :(
					printf("\nBETA=%lf - Is this OK?\nEnter 'y' or 'Y' to continue\nEnter 'a' or 'A' to try again.\n\n",BETA);

					char key;
					while(1)
					{
						key=_getch();
						if (key=='y'||key=='Y'||key=='a'||key=='A')
						{
							break;
						}							
					}	
					if (key=='y'||key=='Y')
					{
						break;
					}	
				}
			}
			else
			{	
				//Quits the program.
				printf("**************************************************\n");
				printf("%i out of %i iterations were completed.\n",iiii,NUM_ITERATIONS_MAX);
				_fcloseall();					
				exit(0);
			}
		}
		
		//Update HHP values.
		hhpSet(HHP,HHP_PRIME);		
	}	
}
