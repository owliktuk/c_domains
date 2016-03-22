void MagneticLayer::countZeemanEnergy()
{
	int N = itsCellsNumber;
	itsZeemanEnergy = 0;
	double L = itsLenght;
	double T = itsThickness;

	double Mx, My, Mz;


	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
			Mx = mx[i][j];
			My = my[i][j];
			Mz = mz[i][j];

			itsZeemanEnergy = itsZeemanEnergy - itsSatMagnetization*(itsExternalField[0]*Mx + itsExternalField[1]*My + itsExternalField[2]*Mz)*L*L*T;

			itsEffectiveField[i][j][0] = itsEffectiveField[i][j][0] + mu0*itsExternalField[0]/itsSatMagnetization;
			itsEffectiveField[i][j][1] = itsEffectiveField[i][j][1] + mu0*itsExternalField[1]/itsSatMagnetization;
			itsEffectiveField[i][j][2] = itsEffectiveField[i][j][2] + mu0*itsExternalField[2]/itsSatMagnetization;
		}
	}
}
