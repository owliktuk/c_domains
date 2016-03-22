#include "../F110.cpp"

void MagneticLayer::countStrayFieldEnergy()
{
	int N = itsCellsNumber;
	itsStrayFieldEnergy = 0;
	double L = itsLenght;
	double T = itsThickness;

	double x, y, z, x1, x2, y1, y2, z1, z2;
	double Wsx1, WsxN, Wsy1, WsyN;

	double DemagnetizingField[3];

	double **Lambda = itsVolumeCharge;
	double **MagneticPotential = itsMagneticPotential;

	double *SigmaX1 = new double [N];
	double *SigmaXN = new double [N];
	double *SigmaY1 = new double [N];
	double *SigmaYN = new double [N];
	
	for (int j=0; j < N; j++){
		SigmaX1[j] = -mx[0][j];
		SigmaXN[j] = mx[N-1][j];
		SigmaY1[j] = -my[j][0];
		SigmaYN[j] = my[j][N-1];
	}

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
			MagneticPotential[i][j] = 0;
			
			for(int k = 0; k < N+1; k++)
			{
				//potancjal od ladunkow objetosciowych
				for(int l = 0; l < N+1; l++)
				{
					MagneticPotential[i][j] = MagneticPotential[i][j] + Lambda[k][l]*itsVolumeChargeIntegral[abs(k-i)][abs(l-j)];

				}

				//potencjal od scianek bocznych
				if(k != N)
				{
					x = i*L;
					y=j*L;
					z=0;

					y1 = k*L;
					y2 = (k+1)*L;
					x1 = 1/2*L;
					z1 = T/2;
					z2 = -T/2;

					Wsx1 = F110(y-y2,z-z2,x-x1) - F110(y-y1,z-z2,x-x1) - F110(y-y2,z-z1,x-x1) + F110(y-y1,z-z1,x-x1);

					x1 = (N+1/2)*L;
					WsxN = F110(y-y2,z-z2,x-x1) - F110(y-y1,z-z2,x-x1) - F110(y-y2,z-z1,x-x1) + F110(y-y1,z-z1,x-x1);

					x1 = (k)*L;
					x2 = (k+1)*L;
					y1 = 1/2*L;
					Wsy1 = F110(x-x2,z-z2,y-y1) - F110(x-x1,z-z2,y-y1) - F110(x-x2,z-z1,y-y1) + F110(x-x1,z-z1,y-y1);

					y1 = (N+1/2)*L;
					WsyN = F110(x-x2,z-z2,y-y1) - F110(x-x1,z-z2,y-y1) - F110(x-x2,z-z1,y-y1) + F110(x-x1,z-z1,y-y1);					

					MagneticPotential[i][j] = MagneticPotential[i][j] + SigmaX1[k]*Wsx1 + SigmaXN[k]*WsxN + SigmaY1[k]*Wsy1 + SigmaYN[k]*WsyN;

				}

			}

			MagneticPotential[i][j] = itsSatMagnetization*MagneticPotential[i][j]/4/PI/mu0;

		}
	}

	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
			DemagnetizingField[0] = -(MagneticPotential[i+1][j+1]+MagneticPotential[i+1][j]-MagneticPotential[i][j+1]-MagneticPotential[i][j])/2/L;
			DemagnetizingField[1] = -(MagneticPotential[i+1][j+1]-MagneticPotential[i+1][j]+MagneticPotential[i][j+1]-MagneticPotential[i][j])/2/L;
			DemagnetizingField[2] = - itsSatMagnetization*mz[i][j]/mu0; //przyblizenie pola demagnetyzacji w kierunku x3

			itsStrayFieldEnergy += mx[i][j]*DemagnetizingField[0] + my[i][j]*DemagnetizingField[1] + mz[i][j]*DemagnetizingField[2];

			itsEffectiveField[i][j][0] = itsEffectiveField[i][j][0] + mu0*DemagnetizingField[0]/itsSatMagnetization;
			itsEffectiveField[i][j][1] = itsEffectiveField[i][j][1] + mu0*DemagnetizingField[1]/itsSatMagnetization;
			itsEffectiveField[i][j][2] = itsEffectiveField[i][j][2] + mu0*DemagnetizingField[2]/itsSatMagnetization;
		}
	}
	
	itsStrayFieldEnergy = -1/2*itsSatMagnetization*itsStrayFieldEnergy*L*L*T;

	delete [] SigmaX1;
	delete [] SigmaXN;
	delete [] SigmaY1;
	delete [] SigmaYN;
	
}
