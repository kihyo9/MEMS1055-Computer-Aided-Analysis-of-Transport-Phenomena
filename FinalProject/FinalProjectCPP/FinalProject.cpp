// FinalProject.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <cmath>
#include "twoDimFlow.h"
#include <algorithm>
#include <ctime> 
#include <chrono>
#include <sstream>
#include <iomanip>
using namespace std;

int main()
{
	auto t = std::time(nullptr);
	auto tm = *std::localtime(&t);

	std::ostringstream oss;
	oss << std::put_time(&tm, "%d-%m-%Y %H-%M-%S");
	auto str = oss.str();

	
	string histfile = "run-" + str + ".txt";
	ofstream myfile;
	myfile.open(histfile);
	cout << "___Starting simulation___\n\n";
	myfile << "___Starting simulation___\n\n";
	

	twoDimFlow tdf;

	//geometry
	const double aLength = 1.5;
	const double bLength = 7.5;
	double totalLength = aLength + bLength;
	const double aHeight = 0.5;
	const double bHeight = 0.5;
	double totalHeight = aHeight + bHeight;



	//solve
	auto start = std::chrono::system_clock::now();
	const double nu = 0.15;
	const double u_m = 0.3;
	double Re = u_m*(aLength + bLength)/nu;
	cout << "Re is " << Re << "\n";
	myfile << "Re is " << Re << "\n";

	const double dx = 0.1;
	const double dy = 0.1;
	double maxu = u_m;
	double maxv = 0;

	double dt_crit = tdf.maxdt(Re, dx, dy, maxu, maxv);
	double dt = 0.8 * dt_crit;
	cout << "dt is " << dt << "\n";
	myfile << "dt is " << dt << "\n";

	double CFL = nu * dt / min(dx, dy);
	cout << "CFL is " << CFL << "\n";
	myfile << "CFL is " << CFL << "\n";
	myfile.close();

		//write initial conditions to file
		tdf.initialConditions(u_m, aLength, bLength, aHeight, bHeight, dx, dy);
		double dummy = 0;

	//number of gridponts
	int dx_steps = int(round(totalLength / dx)) + 1;
	int dy_steps = int(round(totalHeight / dy)) + 1;
	int block_xsteps = int(round(aLength / dx)) + 1;
	int block_ysteps = int(round(bHeight / dy)) + 1;

	//large array
	int size = dx_steps * dy_steps;
	double** coeff = new double* [size];
	for (int i = 0; i < size; i++)
		coeff[i] = new double[size] {};

	int timesteps = tdf.mainSolver(dx, dy, dt, dx_steps, dy_steps, block_xsteps, block_ysteps, Re, histfile, coeff);

	//performance
	ofstream perf;
	perf.open("perf.txt", std::ios_base::app);
	auto finish = std::chrono::system_clock::now();
	std::time_t end_time = chrono::system_clock::to_time_t(finish);
	std::time_t start_time = chrono::system_clock::to_time_t(start);
	chrono::duration<double> elapsed = finish - start;

	perf << "Start time: " << std::ctime(&start_time) << "\n";
	perf << "dx = " << dx << ", " << "dy = " << dy << ", " << "dt = " << dt << ", " << "timesteps = " << timesteps << ", " << "\n";
	perf << elapsed.count() << "s\n\n";
	perf.close();


	// exit
	cout << "\nPress any key to exit... ";

	cin.get();

	return 0;
}