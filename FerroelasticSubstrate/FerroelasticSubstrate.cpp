#include "FerroelasticSubstrate.hpp"

FerroelasticSubstrate::FerroelasticSubstrate(int CellsNumber, int DomainSize, double Temperature, int crystalID)
{
	itsCellsNumber = CellsNumber;
	int N = itsCellsNumber;
	itsDomainSize = DomainSize;
	itsCrystalID = crystalID;

	itsDeformation = new double **[N];
	for (int j=0; j < N; j++)
	{
		itsDeformation[j] = new double * [N];
		for (int k=0; k < N; k++){
			itsDeformation[j][k] = new double[3];
		}
	}

	switch (itsCrystalID)
	{
		case 1: //GMO

			itsDeformation11Curve.resize(5);
			itsDeformation11Curve[0] = 0.00505;
			itsDeformation11Curve[1] = -8.35E-6;
			itsDeformation11Curve[2] = 8.969E-8;
			itsDeformation11Curve[3] = -1.653E-9;
			itsDeformation11Curve[4] = 3.158E-12;

			itsDeformation22Curve.resize(5);
			itsDeformation22Curve[0] = 0.00248;
			itsDeformation22Curve[1] = -4.34E-6;
			itsDeformation22Curve[2] = 1.94E-9;
			itsDeformation22Curve[3] = -3.996E-10;
			itsDeformation22Curve[4] = -1.587E-14;

			itsDeformation12Curve.resize(5);
			itsDeformation12Curve[0] = 0.00133;
			itsDeformation12Curve[1] = 1.937E-6;
			itsDeformation12Curve[2] = -7.857E-8;
			itsDeformation12Curve[3] = 6.656E-10;
			itsDeformation12Curve[4] = -2.82E-12;
			break;

		case 2: //LCS

			itsDeformation12Curve.resize(6);
			itsDeformation12Curve[0] = 0.00336;
			itsDeformation12Curve[1] = 154.71791;
			itsDeformation12Curve[2] = -2.67E-7;
			itsDeformation12Curve[3] = 5.54232;
			itsDeformation12Curve[4] = -3.242E-4;
			itsDeformation12Curve[5] = 34.83432;

			itsDeformation11Curve.resize(3);
			itsDeformation11Curve[0] = 0.00552;
			itsDeformation11Curve[1] = -2.74E-6;
			itsDeformation11Curve[2] = -26.69951;

			itsDeformation22Curve.resize(3);
			itsDeformation22Curve[0] = 0.00202;
			itsDeformation22Curve[1] = -4.88E-6;
			itsDeformation22Curve[2] = -34.235;
			break;
			
		case 3: //KDP-

			itsDeformation12Curve.resize(3);
			itsDeformation12Curve[0] = 0.0039;
			itsDeformation12Curve[1] = -3.76E-14;
			itsDeformation12Curve[2] = -5;


			itsDeformation11Curve.resize(3);
			itsDeformation11Curve[0] = 0.00168;
			itsDeformation11Curve[1] = -8.427E-11;
			itsDeformation11Curve[2] = -7.456;
			break;
	}
	
	countDeformation(Temperature);
}

void FerroelasticSubstrate::getDeformation(int i, int j, double (&Def)[3]) const
{
	Def[0] = itsDeformation[i][j][0];
	Def[1] = itsDeformation[i][j][1];
	Def[2] = itsDeformation[i][j][2];
}

#include "countDeformation.cpp"
  
