/*
 * godunov.cpp
 *
 *  Created on: Mar 10, 2017
 *      Author: siddhant
 */


//Godunov 1st order accurate scheme to numerically simulate 1D Dam Break Problem

#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <algorithm>

//Defining number of grid points and constants
const int N = 101;
const int L = 100;
const double g = 9.81;
const double dt = 0.001;
const double t = 4;

using namespace std;

//Inline function to compute the maximum value of each row
inline double max_vec(vector<vector<double>> a)
{
	for (int i = 0; i < N; ++i)
	{
		if (a[0][i] > a[1][i])
		{
			return a[0][i];
		}
		else
		{
			return a[1][i];
		}
	}
};

int main(int argc, char* argv[])
{

double dx;
dx = L/(N - 1);

vector<vector<double>> Qold(2, vector<double>(N, 0.0));
vector<vector<double>> Qnew(2, vector<double>(N, 0.0));
vector<double> c(N, 0.0);
vector<vector<double>> eigen(2, vector<double>(N, 0.0));
vector<vector<double>> Eold(2, vector<double>(N, 0.0));
vector<vector<double>> F1(2, vector<double>(N, 0.0));
vector<vector<double>> F2(2, vector<double>(N, 0.0));

vector<double> Unew(N, 0.0);
vector<double> alpha(N, 0.0);

//Defining initial value of h and v

for (int i = 0; i < N; ++i)
{
	if (i*dx < 50)
	{
		Qold[0][i] = 10;
		Qold[1][i] = 0;
	}
	else
	{
		Qold[0][i] = 1;
		Qold[1][i] = 0;
	}
}

double tinitial = 0.0;

while (tinitial <= t)
{
	tinitial = tinitial + dt;
	for (int i = 0; i < N; ++i)
	{
		c[i] = sqrt(g*Qold[0][i]);
	}

	//Calculating the Eigen values
	for (int i = 0; i < N; ++i)
	{
		eigen[0][i] = (Qold[1][i]/Qold[0][i]) + c[i];
		eigen[1][i] = (Qold[1][i]/Qold[0][i]) - c[i];
	}

	//Computing the absolute values
	for (unsigned int i = 0; i < 2; ++i)
	{
		for (unsigned int j = 0; j < N; ++j)
		{
			if (eigen[i][j] < 0){
				eigen[i][j] *= -1;}
		}
	}

	for (int i = 0; i < N; ++i)
	{
		alpha[i] = max_vec(eigen);
	}

	for (int i = 0; i < N; ++i)
	{
		Eold[0][i] = Qold[1][i];
		Eold[1][i] = (((Qold[1][i] * Qold[1][i])) / Qold[0][i]) + (0.5 * g * Qold[0][i] * Qold[0][i]);
	}

	for (int i = 0; i < N - 1; ++i)
	{
		F1[0][i] = 0.5*(Eold[0][i] + Eold[0][i+1]) - 0.5*max(alpha[i], alpha[i+1])*(Qold[0][i+1] - Qold[0][i]);
		F1[1][i] = 0.5*(Eold[1][i] + Eold[1][i+1]) - 0.5*max(alpha[i], alpha[i+1])*(Qold[1][i+1] - Qold[1][i]);
	}

	for (int i = 1; i < N; ++i)
	{
		F2[0][i] = 0.5*(Eold[0][i] + Eold[0][i-1]) - 0.5*max(alpha[i], alpha[i-1])*(Qold[0][i] - Qold[0][i-1]);
		F2[1][i] = 0.5*(Eold[1][i] + Eold[1][i-1]) - 0.5*max(alpha[i], alpha[i-1])*(Qold[1][i] - Qold[1][i-1]);
	}

	for (int i = 1; i < N - 1; ++i)
	{
		Qnew[0][i] = Qold[0][i] - (dt/dx)*(F1[0][i] - F2[0][i]);
		Qnew[1][i] = Qold[1][i] - (dt/dx)*(F1[1][i] - F2[1][i]);
	}

	for (int l = 0; l < 2; ++l)
	{
		for (int i = 0; i < 1; ++i)
		{
			Qnew[l][i] = Qold[l][i] - (dt/dx)*(F1[l][i + 1] - F1[l][i]);
		}
	}

	for (int l = 0; l < 2; ++l)
	{
		for (int i = N - 1; i < N; ++i)
		{
			Qnew[l][i] = Qold[l][i] - (dt/dx)*(F2[l][i] - F2[l][i - 1]);
		}
	}

	for (int l = 0; l < 2; ++l)
	{
		for (int i = 0; i < N; ++i)
		{
			Qold[l][i] = Qnew[l][i];
		}
	}

	for (int i = 0; i < N; ++i)
	{
		Unew[i] = Qnew[1][i]/Qold[0][i];
	}

}

//Writing the results for final water flow profile and velocity profile of the water
ofstream write_Qnew("Qnew.dat");
for (int i = 0; i < N; ++i)
{
	write_Qnew << Qnew[0][i] << "\n";
}

ofstream write_Unew("Unew.dat");
for (int i = 0; i < N; ++i)
{
	write_Unew << Unew[i] << "\n";
}


return 0;
}




