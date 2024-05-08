void read(void)
{
	rewind(fp_j2);

	int size=(N_MAX_X-2)*(N_MAX_Y-2);
	for (int j=0;j<=(N_MAX_X-2)*(N_MAX_Y-2)-1;j++)
	{
		printf("%i out of %i\n",j+1,size);
		for (int k=0;k<=(N_MAX_X-2)*(N_MAX_Y-2)-1;k++)
		{		
			fscanf(fp_j2, "%lf ",&J_INVERSE[j][k]);
		}
	}

	fclose(fp_j2);
}
