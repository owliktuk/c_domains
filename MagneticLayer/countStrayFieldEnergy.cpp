#include "../F110.cpp"

void MagneticLayer::countStrayFieldEnergy()
{
	int N = itsCellsNumber;
	itsStrayFieldEnergy = 0;
	double L = itsLenght;
	double T = itsThickness;

	double DemagnetizingField[3];

	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
			DemagnetizingField[0] = 0;
			DemagnetizingField[1] = 0;
			DemagnetizingField[2] = - itsSatMagnetization*mz[i][j]/mu0;

			itsStrayFieldEnergy += mx[i][j]*DemagnetizingField[0] + my[i][j]*DemagnetizingField[1] + mz[i][j]*DemagnetizingField[2];

			itsEffectiveField[i][j][0] = itsEffectiveField[i][j][0] + mu0*DemagnetizingField[0]/itsSatMagnetization;
			itsEffectiveField[i][j][1] = itsEffectiveField[i][j][1] + mu0*DemagnetizingField[1]/itsSatMagnetization;
			itsEffectiveField[i][j][2] = itsEffectiveField[i][j][2] + mu0*DemagnetizingField[2]/itsSatMagnetization;
		}
	}
	
	itsStrayFieldEnergy = -1/2*itsSatMagnetization*itsStrayFieldEnergy*L*L*T;

	
}
