void cross(double (&W1)[3], double (&W2)[3], double (&result)[3])
{
	result[0] = W1[1]*W2[2] - W1[2]*W2[1];
	result[1] = W1[2]*W2[0] - W1[0]*W2[2];
	result[2] = W1[0]*W2[1] - W1[1]*W2[0];
}