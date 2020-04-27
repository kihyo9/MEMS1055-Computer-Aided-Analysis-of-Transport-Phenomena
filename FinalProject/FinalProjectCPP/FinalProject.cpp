// FinalProject.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <cmath>
#include "homework6.h"
#include "twoDimFlow.h"
#include <algorithm>
using namespace std;
double* initialCondition(double*, double*, double);

void testWriting() {
	ofstream myfile;
	cout << "Writing numbers file...\n";
	myfile.open("numbers.txt");

	int xMin = -10, xMax = 10;
	double dx = 0.01;
	int tMin = 0, tMax = 10;
	double dt = 1.2;
	myfile << xMin << " " << xMax << " " << dx << "\n";
	myfile << tMin << " " << tMax << " " << dt << "\n";

	int timestepCount = 0;
	for (double i = tMin; i <= tMax; i += dt) {
		myfile << "Timestep_#" << timestepCount++ << " " << i << "s" << "\n";
		for (double j = xMin; j <= xMax; j += dx) {
			double amp = (0.2 * i + 1);
			double freq = (3. / (1 + 0.4 * i));
			double bias = (i / 4.);
			double num = amp * sin(j * freq / (2 * 3.1416)) + bias;

			if (j == 5) {
				myfile << num;
			}
			else {
				myfile << num << " ";
			}
		}

		myfile << "\n";
	}
	cout << "Finished writing numbers file...\n\n";
	myfile.close();
}

int main()
{
	cout << "___Starting simulation___\n\n";

	//testWriting();
	//homework6 hw6;
	//hw6.solve();

	twoDimFlow tdf;

	//geometry
	const double aLength = 1.5;
	const double bLength = 7.5;
	double totalLength = aLength + bLength;
	const double aHeight = 0.5;
	const double bHeight = 0.5;
	double totalHeight = aHeight + bHeight;



	//solve
	const double nu = 0.15;
	const double u_m = 0.03;
	double Re = u_m*(aLength + bLength)/nu;
	cout << "Re is " << Re << "\n";

	const double dx = 0.1;
	const double dy = 0.1;
	double maxu = u_m;
	double maxv = 0;

	double dt_crit = tdf.maxdt(Re, dx, dy, maxu, maxv);
	double dt = 0.8 * dt_crit;
	cout << "dt is " << dt << "\n";

	double CFL = nu * dt / min(dx, dy);
	cout << "CFL is " << CFL << "\n";

		//write initial conditions to file
		tdf.initialConditions(u_m, aLength, bLength, aHeight, bHeight, dx, dy);
		double dummy = 0;

	//number of gridponts
	int dx_steps = int(round(totalLength / dx)) + 1;
	int dy_steps = int(round(totalHeight / dy)) + 1;
	int block_xsteps = int(round(aLength / dx)) + 1;
	int block_ysteps = int(round(bHeight / dy)) + 1;

	tdf.mainSolver(dx, dy, dt, dx_steps, dy_steps, block_xsteps, block_ysteps);


	// exit
	cout << "\nPress any key to exit... ";

	cin.get();

	return 0;
}


// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
