#include <cmath>

double F111(double x, double y, double z)
{
	double u = x*x;
	double v = y*y;
	double w = z*z;
	double r = sqrt(u+v+w);

	double Lx = atanh(x/r);
	double Ly = atanh(y/r);
	double Lz = atanh(z/r);
	double Px = x*atan(y*z/(x*r));
	double Py = y*atan(x*z/(y*r));
	double Pz = z*atan(y*x/(z*r));

	if(isnan(Lx)) Lx = 0;
	if(isnan(Ly)) Ly = 0;
	if(isnan(Lz)) Lz = 0;
	if(isnan(Px)) Px = 0;
	if(isnan(Py)) Py = 0;
	if(isnan(Pz)) Pz = 0;

	if(isinf(Lx)) Lx = 0;
	if(isinf(Ly)) Ly = 0;
	if(isinf(Lz)) Lz = 0;
	if(isinf(Px)) Px = 0;
	if(isinf(Py)) Py = 0;
	if(isinf(Pz)) Pz = 0;


	double F = x*y*Lz + x*z*Ly + y*z*Lx - 1/2*(x*Px + y*Py + z*Pz);

	return F;
}