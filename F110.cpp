#include <cmath>

double F110(double x, double y, double z)
{
	double u = x*x;
	double v = y*y;
	double w = z*z;
	double r = sqrt(u+v+w);

	double Lx = atanh(x/r);
	double Ly = atanh(y/r);
	double Pz = z*atan(y*x/(z*r));

	if(isnan(Lx)) Lx = 0;
	if(isnan(Ly)) Ly = 0;
	if(isnan(Pz)) Pz = 0;

	if(isinf(Lx)) Lx = 0;
	if(isinf(Ly)) Ly = 0;
	if(isinf(Pz)) Pz = 0;


	double F = y*Lx + x*Ly - Pz;

	return F;
}