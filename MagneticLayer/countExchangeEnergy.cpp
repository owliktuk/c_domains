void MagneticLayer::countExchangeEnergy()
{
	int N = itsCellsNumber;
	itsExchangeEnergy = 0;
	double L = itsLenght;
	double T = itsThickness;
	
	double Dmxdx, Dmxdy, Dmydx, Dmydy, Dmzdx, Dmzdy;


	for(int i = 2; i < N-2; i++)
	{
		for(int j = 2; j < N-2; j++)
		{
			Dmxdx = (-mx[i-2][j]+16*mx[i-1][j]-30*mx[i][j]+16*mx[i+1][j]-mx[i+2][j])/12/pow(itsLenght,2);
			Dmxdy = (-mx[i][j-2]+16*mx[i][j-1]-30*mx[i][j]+16*mx[i][j+1]-mx[i][j+2])/12/pow(itsLenght,2);
			Dmydx = (-my[i-2][j]+16*my[i-1][j]-30*my[i][j]+16*my[i+1][j]-my[i+2][j])/12/pow(itsLenght,2);
			Dmydy = (-my[i][j-2]+16*my[i][j-1]-30*my[i][j]+16*my[i][j+1]-my[i][j+2])/12/pow(itsLenght,2);
			Dmzdx = (-mz[i-2][j]+16*mz[i-1][j]-30*mz[i][j]+16*mz[i+1][j]-mz[i+2][j])/12/pow(itsLenght,2);
			Dmzdy = (-mz[i][j-2]+16*mz[i][j-1]-30*mz[i][j]+16*mz[i][j+1]-mz[i][j+2])/12/pow(itsLenght,2);

			itsExchangeEnergy = itsExchangeEnergy - itsExchangeConstant*(mx[i][j]*(Dmxdx+Dmxdy) + my[i][j]*(Dmydy+Dmydx) + mz[i][j]*(Dmzdx+Dmzdy))*L*L*T;

			itsEffectiveField[i][j][0] = itsEffectiveField[i][j][0] + 2*mu0*itsExchangeConstant*(Dmxdx + Dmxdy)/pow(itsSatMagnetization,2);
			itsEffectiveField[i][j][1] = itsEffectiveField[i][j][1] + 2*mu0*itsExchangeConstant*(Dmydy + Dmydx)/pow(itsSatMagnetization,2);
			itsEffectiveField[i][j][2] = itsEffectiveField[i][j][2] + 2*mu0*itsExchangeConstant*(Dmzdx + Dmzdy)/pow(itsSatMagnetization,2);
		}
	}
}
