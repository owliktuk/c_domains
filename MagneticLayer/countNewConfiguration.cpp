#include "../cross.cpp"

void MagneticLayer::countNewConfiguration(MagneticLayer &new_ML, double alpha)
{
	int N = itsCellsNumber;

	double M[3], heff[3], result[3], result_p[3];

	double new_mx, new_my, new_mz;
	double norma;
	double Phi, Theta;
	
	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
			M[0] = mx[i][j];
			M[1] = my[i][j];
			M[2] = mz[i][j];
			
			heff[0] = itsEffectiveField[i][j][0];
			heff[1] = itsEffectiveField[i][j][1];
			heff[2] = itsEffectiveField[i][j][2];

			cross(M,heff,result_p);
			cross(M,result_p,result);

			new_mx = mx[i][j] - alpha*result[0];
			new_my = my[i][j] - alpha*result[1];
			new_mz = mz[i][j] - alpha*result[2];

			norma = sqrt(new_mx*new_mx+new_my*new_my+new_mz*new_mz);
			new_mx = new_mx/norma;
			new_my = new_my/norma;
			new_mz = new_mz/norma;	

			new_ML.setConfiguration(i,j,new_mx,new_my,new_mz);
		}
	}
}
