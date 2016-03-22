double MagneticLayer::countTotalEnergy(FerroelasticSubstrate &FS)
{
	int N = itsCellsNumber;
	itsTotalEnergy = 0;
	
	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
			itsEffectiveField[i][j][0] = 0;
			itsEffectiveField[i][j][1] = 0;
			itsEffectiveField[i][j][2] = 0;
		}
	}

	countExchangeEnergy();
	countStrayFieldEnergy();
	countZeemanEnergy();
	countMagnetoelasticEnergy(FS);
	//countAnisotropyEnergy();

	itsTotalEnergy = itsExchangeEnergy + itsZeemanEnergy + itsMagnetoelasticEnergy + itsStrayFieldEnergy;

	return itsTotalEnergy;	
}
