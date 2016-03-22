void MagneticLayer::countMagnetization(double (&M)[4], int option)
{
	int N = itsCellsNumber;

	double Mx = 0;
	double My = 0;
	double Mz = 0;

	if(option == 0)
	{

		for(int i = 0; i < N; i++)
		{
			for(int j = 0; j < N; j++)
			{
				Mx += mx[i][j];
				My += my[i][j];
				Mz += mz[i][j];
			}
		}

		itsMagnetization[0] = itsSatMagnetization*sqrt(Mx*Mx+My*My+Mz*Mz)/(N*N)/mu0;
		itsMagnetization[1] = itsSatMagnetization*Mx/(N*N)/mu0;
		itsMagnetization[2] = itsSatMagnetization*My/(N*N)/mu0;
		itsMagnetization[3] = itsSatMagnetization*Mz/(N*N)/mu0;

	} else {

		for(int i = N/4; i < 3*N/4; i++)
		{
			for(int j = N/4; j < 3*N/4; j++)
			{
				Mx += mx[i][j];
				My += my[i][j];
				Mz += mz[i][j];
			}
		}

		itsMagnetization[0] = 4*itsSatMagnetization*sqrt(Mx*Mx+My*My+Mz*Mz)/(N*N)/mu0;
		itsMagnetization[1] = 4*itsSatMagnetization*Mx/(N*N)/mu0;
		itsMagnetization[2] = 4*itsSatMagnetization*My/(N*N)/mu0;
		itsMagnetization[3] = 4*itsSatMagnetization*Mz/(N*N)/mu0;
	}

	M[0] = itsMagnetization[0];
	M[1] = itsMagnetization[1];
	M[2] = itsMagnetization[2];
	M[3] = itsMagnetization[3];

}	