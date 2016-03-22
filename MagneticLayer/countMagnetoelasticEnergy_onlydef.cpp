void MagneticLayer::countMagnetoelasticEnergy(FerroelasticSubstrate &FS)
{
	int N = itsCellsNumber;
	itsMagnetoelasticEnergy = 0;
	double L = itsLenght;
	double T = itsThickness;

	double Mx, My, Mz;
	double itsMagnetoelasticField[3];

	double B1 = itsMagnetoelasticConstant1;
	double B2 = itsMagnetoelasticConstant2;
	double E1 = itsElasticConstant11;
	double E2 = itsElasticConstant12;
	double E4 = itsElasticConstant44;

	double Deformation[3];


	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
			Mx = mx[i][j];
			My = my[i][j];
			Mz = mz[i][j];

			FS.getDeformation(i,j,Deformation);

			itsMagnetoelasticEnergy = itsMagnetoelasticEnergy + B1*(Mx*Mx*Deformation[0] + My*My*Deformation[1]) - (E2/E1*B1*(Deformation[0]+Deformation[1]))*Mz*Mz + B2*Mx*My*Deformation[2];

			itsMagnetoelasticField[0] = 2*B1*Mx*Deformation[0] + B2*My*Deformation[2];
			itsMagnetoelasticField[1] = 2*B1*My*Deformation[1] + B2*Mx*Deformation[2];
			itsMagnetoelasticField[2] = -2*(E2/E1*B1*(Deformation[0]+Deformation[1]))*Mz;

			itsEffectiveField[i][j][0] = itsEffectiveField[i][j][0] - mu0*itsMagnetoelasticField[0]/pow(itsSatMagnetization,2);
			itsEffectiveField[i][j][1] = itsEffectiveField[i][j][1] - mu0*itsMagnetoelasticField[1]/pow(itsSatMagnetization,2);
			itsEffectiveField[i][j][2] = itsEffectiveField[i][j][2] - mu0*itsMagnetoelasticField[2]/pow(itsSatMagnetization,2);

		}
	}
	itsMagnetoelasticEnergy = itsMagnetoelasticEnergy*L*L*T;
}
