double MagneticLayer::countEpsilon()
{
	int N = itsCellsNumber;

	double M[3], heff[3], epsilon[3], eps;

	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
			M[0] = mx[i][j];
			M[1] = my[i][j];
			M[2] = mz[i][j];

			heff[0] = itsEffectiveField[i][j][0];
			heff[1] = itsEffectiveField[i][j][1];
			heff[2] = itsEffectiveField[i][j][2];

			cross(M,heff,epsilon);
			eps = sqrt(epsilon[0]*epsilon[0]+epsilon[1]*epsilon[1]+epsilon[2]*epsilon[2]);

			if(eps > 1E-4)
			{
				cout << "eps " << eps << "\n";				
				return eps;
			}

		}
	}
	
	return eps;
	


}
