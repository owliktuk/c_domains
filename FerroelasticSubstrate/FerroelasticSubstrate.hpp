#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

class FerroelasticSubstrate
{
      public:
             FerroelasticSubstrate(int CellsNumber, int DomainSize, double Temperature, int crystalID);
	     ~FerroelasticSubstrate() { delete itsDeformation; };

	     void countDeformation(double Temperature);
	     void getDeformation(int i, int j, double (&Def)[3]) const;
	     void setDomainSize(double DomainSize) { itsDomainSize = DomainSize; }

      private:
	      int itsCrystalID;
              double ***itsDeformation;
	      int itsCellsNumber;
	      double itsDomainSize;
	      vector<double> itsDeformation11Curve;
	      vector<double> itsDeformation12Curve;
	      vector<double> itsDeformation22Curve;
              
};
