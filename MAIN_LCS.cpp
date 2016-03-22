/* Application counts tempearture dependence of longitudinal magnetization and magnetic configuration of thin ferromagnetic film on 
 * LiCsSO4 ferroelastic substrate.
*/

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
	double theta = 0;	//out-of-plane direction of H
	double phi = 180;	//in-plane direction of H
	double T = 210;  	//T in Kelvins
	double alpha = 1E-5;	
	double koniec, epsilon, epsilon0, koniec_eps;
	double Hx, Magnetization[4];
	int failed = 0;
	double Energy, new_Energy;
	
	std::ofstream plik;	
	
	//MagneticLayer(CellsNumber, Lenght, Thickness, ExchangeConstant, AnisotropyConstant1, 
	//AnisotropyConstant2, MagnetoelasticConstant1, MagnetoelasticConstant2, 
	//Elastic11, Elastic12, Elastic44, SatMagnetization);
	MagneticLayer ML = MagneticLayer(100, 1E-8, 5E-9, 13E-12, -3694, -0.23E4, 1.24E6, 1.24E6, 2.495E11, 1.1209E11, 0.68705E11, 1.005);
	MagneticLayer new_ML = MagneticLayer(100, 1E-8, 5E-9, 13E-12, -3694, -0.23E4, 1.24E6, 1.24E6, 2.495E11, 1.1209E11, 0.68705E11, 1.005);

	//FerroelasticSubstrate(CellsNumber, DomainSize, Temperature, crystalID);
	FerroelasticSubstrate FS = FerroelasticSubstrate(100, 50, T, 2);

	int relaxed = 0;
	int ilosc_drog = 4;
	for(int droga = 1; droga < ilosc_drog; droga++)
	{

		
		if(T == 100) T = 110;
		if(T == 215) T = 210;

		while(T > 100 && T < 215)
		{

			plik.open("LCS_wyniki.txt", std::ofstream::app);
			
			koniec = 0;
			koniec_eps = 0;
			epsilon0 = -5;
			failed = 0;

			if(T == 205 && relaxed == 0)
			{
				T = 210;
				relaxed = 1;
			}

			if(T == 210 && relaxed == 0) Hx = 0; else Hx = 8000; //Magnetic field strenght in A/m


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

					if(epsilon == epsilon0)
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

			//longitudinal magnetization

			ML.countMagnetization(Magnetization, 0);
			plik << T << "\t" << Magnetization[0] << "\t" << Magnetization[1] << "\t" << Magnetization[2] << "\t" << Magnetization[3] << "\t";
			ML.countMagnetization(Magnetization, 1);
			plik << Magnetization[0] << "\t" << Magnetization[1] << "\t" << Magnetization[2] << "\t" << Magnetization[3] << "\n";
			cout << "T: " << T << " M: " << Magnetization[0] << "\t" << Magnetization[1] << "\t" << Magnetization[2] << "\t" << Magnetization[3] << "\n";

			//magnetic configuration at T=210K and T=110K

			if(T == 210 && droga > ilosc_drog-2)
			{
				ML.exportConfiguration(0);
				cout << "eksportuje 210... \n";
			} else if(T == 110 && droga > ilosc_drog-2)
			{
				ML.exportConfiguration(1);
				cout << "eksportuje 110... \n";
			}

			if(T > 170) T += pow(-1,droga)*5; else T += pow(-1,droga)*10;

			plik.close();

		}
	}
		

	return 0;
}
