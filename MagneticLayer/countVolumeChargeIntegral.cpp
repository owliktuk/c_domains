#include "../F111.cpp"

void MagneticLayer::countVolumeChargeIntegral()
{
	int N = itsCellsNumber;
	double L = itsLenght;
	double T = itsThickness;

	double x, y, z, x1, x2, y1, y2, z1, z2;	

	for(int k=0; k < N+1; k++)
	{
		for(int l=0; l < N+1; l++)
		{

			//volume charges
			x2=(k+1/2)*L;
			x1 = (k-1/2)*L;
			x = 0;

			y2=(l+1/2)*L;
			y1 = (l-1/2)*L;
			y = 0;

			z2 = -T/2;
			z1 = T/2;
			z = 0;

			itsVolumeChargeIntegral[k][l] = F111(x-x2,y-y2,z-z2) - F111(x-x1,y-y2,z-z2) - F111(x-x2,y-y1,z-z2) + F111(x-x1,y-y1,z-z2) - F111(x-x2,y-y2,z-z1) + F111(x-x1,y-y2,z-z1) + F111(x-x2,y-y1,z-z1) - F111(x-x1,y-y1,z-z1);


		}
	}

}