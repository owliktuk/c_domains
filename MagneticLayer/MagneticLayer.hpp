#include "../FerroelasticSubstrate/FerroelasticSubstrate.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <cstdio>

const double mu0 = 12.566E-7;
const double PI = 3.14159265359;

using namespace std;
class MagneticLayer
{
      public:
             MagneticLayer(int CellsNumber, double Lenght, double Thickness,
double ExchangeConstant, double AnisotropyConstant1, double AnisotropyConstant2, double MagnetoelasticConstant1, double MagnetoelasticConstant2, double Elastic11, double Elastic12, double Elastic44, double SatMagnetization);	     
	     ~MagneticLayer();

             void setConfiguration(int i, int j, double mx, double my, double mz);
             void setExternalField(double Hx, double Hy, double Hz);
	     void exportConfiguration(int w) const;
             
             double getConfiguration(int i, int j, int component) const;
             double getEffectiveField(int i, int j, int k) const { return itsEffectiveField[i][j][k]; }
             
             double countTotalEnergy(FerroelasticSubstrate &FS);
	     void countNewConfiguration(MagneticLayer &Temp_ML, double alpha);
	     void copyConfiguration(const MagneticLayer &);
	     double countEpsilon();
	     void countMagnetization(double (&M)[4], int opt);

      private:
	      double **mx;
	      double **my;
	      double **mz;
	      int itsCellsNumber;
	      double itsLenght;
	      double itsThickness;
	      double itsExchangeConstant;
	      double itsAnisotropyConstant1;
	      double itsAnisotropyConstant2;
	      double itsMagnetoelasticConstant1;
	      double itsMagnetoelasticConstant2;
	      double itsElasticConstant11;
	      double itsElasticConstant12;
	      double itsElasticConstant44;
	      double itsSatMagnetization;
	      double itsMagnetization[4];
	      
	      double itsExchangeEnergy;
	      double itsAnisotropyEnergy;
	      double itsZeemanEnergy;
	      double itsStrayFieldEnergy;
	      double **itsVolumeChargeIntegral;
	      double itsMagnetoelasticEnergy;
	      double itsTotalEnergy;
	      double ***itsEffectiveField;
	      double itsExternalField[3];

	      double **itsMagneticPotential;
	      double **itsVolumeCharge;
	      
	      void countExchangeEnergy();
	      void countAnisotropyEnergy();
	      void countZeemanEnergy();
	      void countStrayFieldEnergy();
	      void countVolumeChargeIntegral();
	      void countMagnetoelasticEnergy(FerroelasticSubstrate &FS);
              
};
