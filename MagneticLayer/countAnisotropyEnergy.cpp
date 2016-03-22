void MagneticLayer::countAnisotropyEnergy()
{
	int N = itsCellsNumber;
	itsAnisotropyEnergy = 0;
	double L = itsLenght;
	double T = itsThickness;

	double Mx, My, Mz;
	double K1 = itsAnisotropyConstant1;
	double K2 = itsAnisotropyConstant2;


	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
			Mx = mx[i][j];
			My = my[i][j];
			Mz = mz[i][j];

			itsAnisotropyEnergy = itsAnisotropyEnergy + (K1*(Mx*Mx*My*My + Mx*Mx*Mz*Mz + My*My*Mz*Mz) + K2*Mx*Mx*My*My*Mz*Mz)*L*L*T;

			itsEffectiveField[i][j][0] = itsEffectiveField[i][j][0] - mu0*(2*K1*(Mx*My*My + Mx*Mz*Mz) + 2*K2*Mx*My*My*Mz*Mz)/pow(itsSatMagnetization,2);
			itsEffectiveField[i][j][1] = itsEffectiveField[i][j][1] - mu0*(2*K1*(My*Mz*Mz + My*Mx*Mx) + 2*K2*My*Mx*Mx*Mz*Mz)/pow(itsSatMagnetization,2);
			itsEffectiveField[i][j][2] = itsEffectiveField[i][j][2] - mu0*(2*K1*(Mz*My*My + Mz*Mx*Mx) + 2*K2*Mz*My*My*Mx*Mx)/pow(itsSatMagnetization,2);
		}
	}
}
