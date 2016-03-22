#include "../F110.cpp"

void MagneticLayer::countStrayFieldEnergy(double **mx, double **my, double **mz)
{
	int N = itsCellsNumber;
	itsStrayFieldEnergy = 0;
	double L = itsLenght;
	double T = itsThickness;

	double x, y, z, x1, x2, y1, y2, z1, z2;
	double Wsx1, WsxN, Wsy1, WsyN, WszT, WszB;

	double DemagnetizingField[3];

	double **Lambda = itsVolumeCharge;
	double ***MagneticPotential = itsMagneticPotential;


	for(int i = 0; i < N+1; i++)
	{
		Lambda[i][0] = 0;
		Lambda[0][i] = 0;
		Lambda[N][i] = 0;
		Lambda[i][N] = 0;
		//obliczanie ładunków objętościowych
		for(int j = 0; j < N-1; j++)
		{			
			if(i < N-1)
				Lambda[i+1][j+1] = -(mx[i+1][j+1]+mx[i+1][j]-mx[i][j+1]-mx[i][j]+my[i+1][j+1]-my[i+1][j]+my[i][j+1]-my[i][j])/L/2;
		}
	}

	for(int i = 0; i < N+1; i++)
	{
		for(int j = 0; j < N+1; j++)
		{
			//zerowanie MagneticPotential
			MagneticPotential[i][j][0] = 0;
			MagneticPotential[i][j][1] = 0;
			
			for(int k = 0; k < N+1; k++)
			{
				for(int l = 0; l < N+1; l++)
				{
					MagneticPotential[i][j][0] = MagneticPotential[i][j][0] + Lambda[k][l]*itsVolumeChargeIntegral[abs(k-i)][abs(l-j)][0];
					MagneticPotential[i][j][1] = MagneticPotential[i][j][1] + Lambda[k][l]*itsVolumeChargeIntegral[abs(k-i)][abs(l-j)][1];

					if(k != N)
					{
						x = i*L;
						y=j*L;
						z=0;

						x2 = (k+1)*L;
						x1 = k*L;
						y2 = (l+1)*L;
						y1 = l*L;
						z2 = -T/2;
						z1 = T/2;

						//WszT = F110(x-x2,y-y2,z-z1) - F110(x-x1,y-y2,z-z1) - F110(x-x2,y-y1,z-z1) + F110(x-x1,y-y1,z-z1);
						//WszB = F110(x-x2,y-y2,z-z2) - F110(x-x1,y-y2,z-z2) - F110(x-x2,y-y1,z-z2) + F110(x-x1,y-y1,z-z2);

						//MagneticPotential[i][j][0] = MagneticPotential[i][j][0] + 2*mz[k][l]*WszT;

						z=T/2;
						WszT = F110(x-x2,y-y2,z-z1) - F110(x-x1,y-y2,z-z1) - F110(x-x2,y-y1,z-z1) + F110(x-x1,y-y1,z-z1);
						WszB = F110(x-x2,y-y2,z-z2) - F110(x-x1,y-y2,z-z2) - F110(x-x2,y-y1,z-z2) + F110(x-x1,y-y1,z-z2);

						MagneticPotential[i][j][1] = MagneticPotential[i][j][1] + mz[k][l]*(WszT-WszB);
					}
				}

			}

			MagneticPotential[i][j][0] = itsSatMagnetization*MagneticPotential[i][j][0]/4/PI/mu0;
			MagneticPotential[i][j][1] = itsSatMagnetization*MagneticPotential[i][j][1]/4/PI/mu0;

		}
	}

	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
			DemagnetizingField[0] = -(MagneticPotential[i+1][j+1][0]+MagneticPotential[i+1][j][0]-MagneticPotential[i][j+1][0]-MagneticPotential[i][j][0])/2/L;
			DemagnetizingField[1] = -(MagneticPotential[i+1][j+1][0]-MagneticPotential[i+1][j][0]+MagneticPotential[i][j+1][0]-MagneticPotential[i][j][0])/2/L;
			DemagnetizingField[2] = -(MagneticPotential[i][j][1] + MagneticPotential[i+1][j][1] + MagneticPotential[i+1][j+1][1] + MagneticPotential[i][j+1][1] - (MagneticPotential[i][j][0] + MagneticPotential[i+1][j][0] + MagneticPotential[i+1][j+1][0] + MagneticPotential[i][j+1][0]))/2/T; //- itsSatMagnetization*mz[i][j]/mu0;

			itsStrayFieldEnergy += mx[i][j]*DemagnetizingField[0] + my[i][j]*DemagnetizingField[1] + mz[i][j]*DemagnetizingField[2];

			itsEffectiveField[i][j][0] = itsEffectiveField[i][j][0] + mu0*DemagnetizingField[0]/itsSatMagnetization;
			itsEffectiveField[i][j][1] = itsEffectiveField[i][j][1] + mu0*DemagnetizingField[1]/itsSatMagnetization;
			itsEffectiveField[i][j][2] = itsEffectiveField[i][j][2] + mu0*DemagnetizingField[2]/itsSatMagnetization;
		}
	}
	
	itsStrayFieldEnergy = -1/2*itsSatMagnetization*itsStrayFieldEnergy*L*L*T;

}
