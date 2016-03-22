void FerroelasticSubstrate::countDeformation(double T)
{
	int N = itsCellsNumber;
	int Ds = itsDomainSize;

	vector<double> A,B,C;
	A = itsDeformation12Curve;
	B = itsDeformation11Curve;
	C = itsDeformation22Curve;

	double u11, u12, u22;
	double u01 = 0;
	double u02 = 0;
	double u03 = 0;

	switch (itsCrystalID)
	{
		case 1: //GMO
			u12 = A[0] + A[1]*T + A[2]*pow(T,2) + A[3]*pow(T,3) + A[4]*pow(T,4);
			u11 = B[0] + B[1]*T + B[2]*pow(T,2) + B[3]*pow(T,3) + B[4]*pow(T,4);
			u22 = C[0] + C[1]*T + C[2]*pow(T,2) + C[3]*pow(T,3) + C[4]*pow(T,4);

			for(int i = 0; i < N; i++)
			{
				for(int j = 0; j < N; j++)
				{
					if(T > 159)
					{
						itsDeformation[i][j][0] = u01;
						itsDeformation[i][j][1] = u02;
						itsDeformation[i][j][2] = u03;
					} else {

						if(int(floor(i/Ds))%2 == 0)
						{
							itsDeformation[i][j][0] = u01-u11;
							itsDeformation[i][j][1] = u02-u22;
							itsDeformation[i][j][2] = u12-u03;
						} else {
							itsDeformation[i][j][0] = u01-u22;
							itsDeformation[i][j][1] = u02-u11;
							itsDeformation[i][j][2] = -u12-u03;
						}
					}
					
					

				}
			}
			break;

		case 2: //LCS

			u12 = A[0] + A[2]*exp((T-A[1])/A[3])+A[4]*exp((T-A[1])/A[5]);
			u11 = B[0] + B[1]*exp(-T/B[2]);
			u22 = C[0] + C[1]*exp(-T/C[2]);
			//u02 = -0.00183;

			for(int i = 0; i < N; i++)
			{
				for(int j = 0; j < N; j++)
				{
					if(T > 201)
					{
						itsDeformation[i][j][0] = u01;
						itsDeformation[i][j][1] = u02;
						itsDeformation[i][j][2] = u03;
					} else {
						itsDeformation[i][j][0] = u01-u11;
						itsDeformation[i][j][1] = u02-u22;

						if(int(floor(j/Ds))%2 == 0)
						{
							itsDeformation[i][j][2] = u03+u12;
						} else {
							itsDeformation[i][j][2] = u03-u12;
						}
					}

				}
			}
			break;
			
		case 3: //KDP
			u12 = A[0] + A[1]*exp(-T/A[2]);
			u11 = B[0] + B[1]*exp(-T/B[2]);

			for(int i = 0; i < N; i++)
			{
				for(int j = 0; j < N; j++)
				{
					if(T > 159)
					{
						itsDeformation[i][j][0] = u01;
						itsDeformation[i][j][1] = u02;
						itsDeformation[i][j][2] = u03;
					} else {

						if(int(floor(i/Ds))%2 == 0)
						{
							itsDeformation[i][j][0] = u01 + (u12+u11);
							itsDeformation[i][j][1] = u02 + (-u12+u11);
							itsDeformation[i][j][2] = u03 + u12;
						} else {
							itsDeformation[i][j][0] = u01 + (-u12+u11);
							itsDeformation[i][j][1] = u02 + (u12+u11);
							itsDeformation[i][j][2] = u03 - u12;
						}
					}
					
					

				}
			}
			break;
			
	}

}
				
