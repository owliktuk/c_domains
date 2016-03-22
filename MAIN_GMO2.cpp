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
	double T = 50;
	double alpha = 1E-5;
	double koniec, epsilon, epsilon0, koniec_eps;
	double Hx, Magnetization[4];
	int failed = 0;
	double Energy, new_Energy;
	
	std::ofstream plik;	    
	
	MagneticLayer ML = MagneticLayer(100, 1E-8, 50E-9, 8E-12, -3694, -0.23E4, 10.2E6, 10.2E6, 3.115E11, 1.257E11, 0.929E11, 0.62);
	MagneticLayer new_ML = MagneticLayer(100, 1E-8, 50E-9, 8E-12, -3694, -0.23E4, 10.2E6, 10.2E6, 3.115E11, 1.257E11, 0.929E11, 0.62);
	FerroelasticSubstrate FS = FerroelasticSubstrate(100, 50, T, 1);

	int relaxed = 0;
	int ilosc_drog = 5;
	for(int droga = 1; droga < ilosc_drog; droga++)
	{

		
		if(T < 50) T = 50;
		if(T > 165) T = 165;

		while(T >= 50 && T <= 165)
		{
			plik.open("GMO_wyniki_phi0_theta0_N100_Ds50_L1-8_D1_u12.txt", std::ofstream::app);
			
			koniec = 0;
			koniec_eps = 0;
			epsilon0 = -5;
			failed = 0;

			if(T > 50 && relaxed == 0)
			{
				T = 50;
				relaxed = 1;
			}

			if(T == 50 && relaxed == 0) Hx = 0; else Hx = 8000;


			ML.setExternalField(Hx*cos(PI*theta/180)*cos(PI*phi/180), Hx*cos(PI*theta/180)*sin(PI*phi/180), Hx*sin(PI*theta/180));
			new_ML.setExternalField(Hx*cos(PI*theta/180)*cos(PI*phi/180), Hx*cos(PI*theta/180)*sin(PI*phi/180), Hx*sin(PI*theta/180));

			if(droga == 1)
			{
				FS.setDomainSize(100); //pojedyncza domena na poczatku
			} else if(droga != 1 && T < 159)
			{
				FS.setDomainSize(55);
			}
			
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

			ML.countMagnetization(Magnetization, 0);
			plik << T << "\t" << Magnetization[0] << "\t" << Magnetization[1] << "\t" << Magnetization[2] << "\t" << Magnetization[3] << "\t";
			ML.countMagnetization(Magnetization, 1);
			plik << Magnetization[0] << "\t" << Magnetization[1] << "\t" << Magnetization[2] << "\t" << Magnetization[3] << "\n";
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

			if(T > 120) T += pow(-1,droga+1)*5; else T += pow(-1,droga+1)*10;

			plik.close();
			
		}
	}

				

	return 0;
}
