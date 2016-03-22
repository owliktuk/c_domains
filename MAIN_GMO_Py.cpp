#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <cstdio>

#include "MagneticLayer/MagneticLayer.hpp"


using std::cout;
using std::cin;

int main()
{
	double theta = 0;
	double phi = 0;
	double T = 165;
	double alpha = 1E-5;
	double koniec, epsilon, epsilon0, koniec_eps;
	double Hx, Magnetization[4];
	int failed = 0;
	double Energy, new_Energy;
	
	std::ofstream plik;	    
	
	MagneticLayer ML = MagneticLayer(100, 1E-8, 5E-9, 13E-12, -3694, -0.23E4, 1.24E6, 1.24E6, 2.495E11, 1.1209E11, 0.68705E11, 1.005);
	MagneticLayer new_ML = MagneticLayer(100, 1E-8, 5E-9, 13E-12, -3694, -0.23E4, 1.24E6, 1.24E6, 2.495E11, 1.1209E11, 0.68705E11, 1.005);
	FerroelasticSubstrate FS = FerroelasticSubstrate(100, 100, T, 1);

	int relaxed = 0;
	int ilosc_drog = 4;
	for(int droga = 1; droga < ilosc_drog; droga++)
	{

		
		if(T < 50) T = 50;
		if(T > 165) T = 165;

		while(T >= 50 && T <= 165)
		{
			plik.open("GMOPy_wyniki_phi0_theta0_N100_Ds100_L1-8_D100.txt", std::ofstream::app);
			
			koniec = 0;
			koniec_eps = 0;
			epsilon0 = -5;
			failed = 0;

			if(T < 165 && relaxed == 0)
			{
				T = 165;
				relaxed = 1;
			}

			if(T == 165 && relaxed == 0) Hx = 0; else Hx = 8000;


			ML.setExternalField(Hx*cos(PI*theta/180)*cos(PI*phi/180), Hx*cos(PI*theta/180)*sin(PI*phi/180), Hx*sin(PI*theta/180));
			new_ML.setExternalField(Hx*cos(PI*theta/180)*cos(PI*phi/180), Hx*cos(PI*theta/180)*sin(PI*phi/180), Hx*sin(PI*theta/180));

			FS.countDeformation(T);

			Energy = ML.countTotalEnergy(FS);

			alpha = 1E-5;
			while(!koniec && failed < 40 && koniec_eps < 50)
			{
				ML.countNewConfiguration(new_ML, alpha);
				new_Energy = new_ML.countTotalEnergy(FS);

				if(new_Energy < Energy)
				{
					ML.copyConfiguration(new_ML);
					Energy = new_Energy;
					alpha = 2*alpha;
					epsilon = ML.countEpsilon();
					failed = 0;

					if(epsilon < 1E-4) koniec = 1;

					if(epsilon == epsilon0) //zabezpieczenie przed zawieszeniem sie obliczen
					{
					   koniec_eps++;
					} else {
					  koniec_eps = 0;
					  epsilon0 = epsilon;
					}
					
				} else {
					failed++;
					alpha = alpha/2;
				}

			}

			ML.countMagnetization(Magnetization);
			plik << T << "\t" << Magnetization[0] << "\t" << Magnetization[1] << "\t" << Magnetization[2] << "\t" << Magnetization[3] << "\n";
			cout << "T: " << T << " M: " << Magnetization[0] << "\t" << Magnetization[1] << "\t" << Magnetization[2] << "\t" << Magnetization[3] << "\n";

			if(T == 50 && droga > ilosc_drog-2)
			{
				ML.exportConfiguration(1);
				cout << "eksportuje 50... \n";
			} else if(T == 165 && droga > ilosc_drog-2)
			{
				ML.exportConfiguration(0);
				cout << "eksportuje 165... \n";
			}

			if(T > 120) T += pow(-1,droga)*5; else T += pow(-1,droga)*10;

			plik.close();
			
		}
	}

				

	return 0;
}
