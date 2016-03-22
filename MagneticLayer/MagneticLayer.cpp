#include "MagneticLayer.hpp"
//#include "../FerroelasticSubstrate/FerroelasticSubstrate.hpp"


MagneticLayer::MagneticLayer(int CellsNumber, double Lenght, double Thickness, double ExchangeConstant, double AnisotropyConstant1, double AnisotropyConstant2, double MagnetoelasticConstant1, double MagnetoelasticConstant2, double Elastic11, double Elastic12, double Elastic44, double SatMagnetization)
{
	itsCellsNumber = CellsNumber;
	itsLenght = Lenght;
	itsThickness = Thickness;
	itsExchangeConstant = ExchangeConstant;
	itsAnisotropyConstant1 = AnisotropyConstant1;
	itsAnisotropyConstant2 = AnisotropyConstant2;
	itsMagnetoelasticConstant1 = MagnetoelasticConstant1;
	itsMagnetoelasticConstant2 = MagnetoelasticConstant2;
	itsElasticConstant11 = Elastic11;
	itsElasticConstant12 = Elastic12;
	itsElasticConstant44 = Elastic44;
	itsSatMagnetization = SatMagnetization;
	itsExternalField[0] = 0;
	itsExternalField[1] = 0;
	itsExternalField[2] = 0;
	
	
	mx = new double * [CellsNumber];
	my = new double * [CellsNumber];
	mz = new double * [CellsNumber];
	for (int j=0; j < CellsNumber; j++){
		mx[j] = new double[CellsNumber];
		my[j] = new double[CellsNumber];
		mz[j] = new double[CellsNumber];
	}

	double Phi, Theta;
	srand (time(NULL));
	for (int i=0; i < CellsNumber; i++){
		for (int j=0; j < CellsNumber; j++){
			Phi =  rand() % 360;
			Theta = rand() % 11 -5;
			mx[i][j] = cos(PI*Theta/180)*cos(PI*Phi/180);
			my[i][j] = cos(PI*Theta/180)*sin(PI*Phi/180);
			mz[i][j] = sin(PI*Theta/180);
		}
	}

	itsEffectiveField = new double **[CellsNumber];
	for (int j=0; j < CellsNumber; j++)
	{
		itsEffectiveField[j] = new double * [CellsNumber];
		for (int k=0; k < CellsNumber; k++){
			itsEffectiveField[j][k] = new double[3];
		}
	}

	
	itsVolumeCharge = new double * [CellsNumber+1];
	itsMagneticPotential = new double * [CellsNumber+1];
	itsVolumeChargeIntegral = new double * [CellsNumber+1];
	for (int j=0; j < CellsNumber+1; j++){
		itsVolumeCharge[j] = new double[CellsNumber+1];
		itsMagneticPotential[j] = new double [CellsNumber+1];
		itsVolumeChargeIntegral[j] = new double [CellsNumber+1];
	}

	countVolumeChargeIntegral();

}

void MagneticLayer::copyConfiguration(const MagneticLayer &rhs)
{

	for (int i=0; i < itsCellsNumber; i++){
		for (int j=0; j < itsCellsNumber; j++){
			mx[i][j] = rhs.getConfiguration(i,j,0);
			my[i][j] = rhs.getConfiguration(i,j,1);
			mz[i][j] = rhs.getConfiguration(i,j,2);
			itsEffectiveField[i][j][0] = rhs.getEffectiveField(i,j,0);
			itsEffectiveField[i][j][1] = rhs.getEffectiveField(i,j,1);
			itsEffectiveField[i][j][2] = rhs.getEffectiveField(i,j,2);
		}
	}
}

MagneticLayer::~MagneticLayer()
{
	delete itsVolumeChargeIntegral;
	delete itsEffectiveField;
	delete itsVolumeCharge;
	delete itsMagneticPotential;
	delete mx;
	delete my;
	delete mz;
}

void MagneticLayer::setExternalField(double Hx, double Hy, double Hz)
{
	itsExternalField[0] = Hx;
	itsExternalField[1] = Hy;
	itsExternalField[2] = Hz;
}

void MagneticLayer::setConfiguration(int i, int j, double new_mx, double new_my, double new_mz)
{
	mx[i][j] = new_mx;
	my[i][j] = new_my;
	mz[i][j] = new_mz;
}

double MagneticLayer::getConfiguration(int i, int j, int component) const
{
	switch(component)
	{
		case 0:
			return mx[i][j];
		case 1:
			return my[i][j];
		case 2:
			return mz[i][j];
	}
}

void MagneticLayer::exportConfiguration(int w) const
{
	ofstream plikx;
	ofstream pliky;
	ofstream plikz;
	if(w == 0)
	{
		plikx.open("configurationX0.txt");
		pliky.open("configurationY0.txt");
		plikz.open("configurationZ0.txt");
		
	} else {

		plikx.open("configurationX.txt");
		pliky.open("configurationY.txt");
		plikz.open("configurationZ.txt");

	}

	for (int i=0; i < itsCellsNumber; i++){
		for (int j=0; j < itsCellsNumber; j++){
			plikx << mx[i][j] << " ";
			pliky << my[i][j] << " ";
			plikz << mz[i][j] << " ";
		}
		plikx << "\n";
		pliky << "\n";
		plikz << "\n";
	}
	plikx.close();
	pliky.close();
	plikz.close();
}

#include "countExchangeEnergy.cpp"
#include "countAnisotropyEnergy.cpp"
#include "countStrayFieldEnergy.cpp"
#include "countZeemanEnergy.cpp"
#include "countMagnetoelasticEnergy.cpp"
#include "countNewConfiguration.cpp"
#include "countVolumeChargeIntegral.cpp"
#include "countTotalEnergy.cpp"
#include "countEpsilon.cpp"
#include "countMagnetization.cpp"
	
  
