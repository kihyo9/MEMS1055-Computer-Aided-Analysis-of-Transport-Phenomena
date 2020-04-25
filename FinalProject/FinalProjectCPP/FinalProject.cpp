// FinalProject.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <cmath>
#include "homework6.h"
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

void solve() {
	homework6 hw6;

	double u = 1;
	double gamma = 1;
	double L = 1;
	double dx = 0.01;
	double t_last = L / u;
	const int dx_steps = 101;
	double dx_data[dx_steps];
	double t_0[dx_steps];
	string fileName = "explicit-data.txt";

	//make dx_data (x-axis)
	for (int i = 0; i < dx_steps; i++) {
		dx_data[i] = i * dx;
	}

	//make initial condition timestep
	initialCondition(dx_data, t_0, L);

	//calculate dt
	double r = 0.25;
	double dt;
	if (gamma == 0.) {
		dt = 2.5e-4;
	}
	else {
		dt = r * pow(dx,2) / gamma;
	}
	const int dt_steps = int(ceil(t_last / dt)) + 1;

	//print simulation conditions
	hw6.printStability(dx, dt, u, gamma, r);

	//solve
	hw6.solveExplicit(dt, dx, gamma, u, L, t_0, dx_steps, dt_steps);

}

double* initialCondition(double* dx_data, double* t_0, double L) {
	for (int i = 0; i < 101; i++) {
		if (dx_data[i] <= L / 4 || dx_data[i] >= 3 * L / 4) {
			t_0[i] = 0;
		}
		else if (dx_data[i] > L / 4 && dx_data[i] < L / 2) {
			t_0[i] = dx_data[i] * 4. / L + -1;
		}
		else {
			t_0[i] = -dx_data[i] * 4. / L + 3;
		}
	}

	return t_0;
}

int main()
{
	cout << "___Starting simulation___\n\n";

	//testWriting();

	solve();


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
