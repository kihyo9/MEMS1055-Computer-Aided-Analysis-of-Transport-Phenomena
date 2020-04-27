#pragma once
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <stdlib.h>
#define MAX_SIZE 7000000
using namespace std;
class twoDimFlow
{
public:
	void mainSolver(double dx, double dy, double dt, int dx_steps, int dy_steps, int block_xsteps, int block_ysteps) {

		//calculate initial vorticies
		solverA(dx, dy, dx_steps, dy_steps, block_xsteps, block_ysteps);

		//hard-coded
		double answer[1100] = {};
		solverC(dx, dy, dx_steps, dy_steps, block_xsteps, block_ysteps, answer); //get vorticities for timestep 0
		//solverB(u, v);
		int x = 0;

		/*
		Start of timestep iterations
		*/		
		//while (!isConverged(w,phi,wm1,phim1)) {
		//	//retained data
		//	double wm1 = 0;
		//	double phim1 = 0;

		//	//new vorticies
		//	solverB(u, v);

		//	//new stream function
		//	solverC(dx_steps, dy_steps);

		//	//new velocities
		//	solverDu(phi);
		//	solverDv(phi);
		//}
		/*
		End of timestep iterations
		*/
	}

	/*
	Purpose: Initialize vorticity via initial velocities
	Complexity: O(xy), simple first-order scheme for each gridpoint
	Equation: w = dv/dx - du/dx
	Notes: Used only once at timestep zero
	*/
	void solverA(double dx, double dy, int dx_steps, int dy_steps, int block_xsteps, int block_ysteps) {
		cout << "Writing initial vorticities...\n";
		ifstream ufile;
		ifstream vfile;
		ofstream vort;
		ufile.open("u.txt");
		vfile.open("v.txt");
		vort.open("vort.txt");

		//hard-coded
		double u1[11] = {};
		double u2[11] = {};
		double u3[11] = {};

		//hard-coded
		double v1[11] = {};
		double v2[11] = {};
		double v3[11] = {};
		

		//walls: i = 0, i = block_xsteps, i = dx_steps, j = 0, j = block_ysteps, j = dy_steps
		for (int i = 0; i < dx_steps; i++) {
			int stop = i < block_xsteps - 1 ? block_ysteps : dy_steps;

			//inlet
			if (i == 0) {
				for (int a = 0; a < block_ysteps; a++) {
					ufile >> u2[a];
					vfile >> v2[a];
				}

				for (int a = 0; a < block_ysteps; a++) {
					ufile >> u3[a];
					vfile >> v3[a];
				}
				//calculate and write the vorticities
				for (int j = 0; j < stop; j++) {
					if (j == stop - 1) {
						vort << -(3 * u2[j] - 4 * u2[j - 1] + u2[j - 2]) / (2 * dy);
					}
					else if (j == 0) {
						vort << -(-3 * u2[j] + 4 * u2[j + 1] - u2[j + 2]) / (2 * dy);
					}
					else {
						vort << (v3[j] - v2[j]) / (dx) - (u2[j + 1] - u2[j - 1]) / (2 * dy);
					}

					if (j < stop - 1) {
						vort << " ";
					}
					else {
						vort << "\n";
					}
				}
			}
			//outlet
			else if (i == dx_steps - 1) {
				//calculate and write the vorticities
				for (int j = 0; j < stop; j++) {
					if (j == stop - 1) {
						vort << -(u3[j] - u3[j-1]) / (dy);
					}
					else if (j == 0) {
						vort << -(u3[j + 1] - u3[j]) / (dy);
					}
					else {
						vort << (v3[j] - v2[j]) / (dx) - (u3[j + 1] - u3[j - 1]) / (2 * dy);
					}

					if (j < stop - 1) {
						vort << " ";
					}
					else {
						vort << "\n";
					}
				}
			}
			//cliff
			else if (i == block_xsteps - 1 || i == block_xsteps - 2) {
				for (int a = 0; a < stop; a++) {
					u1[a] = u2[a];
					u2[a] = u3[a];
					v1[a] = v2[a];
					v2[a] = v3[a];
					ufile >> u3[a];
					vfile >> v3[a];
				}

				//line before cliff needs the u3 and v3 to read extra elements
				if (i == block_xsteps - 2) {
					for (int a = stop; a < dy_steps; a++) {
						ufile >> u3[a];
						vfile >> v3[a];
					}
				}


				//calculate and write the vorticities
				for (int j = 0; j < stop; j++) {
					if (j > dy_steps - block_ysteps + 1) {
						vort << (v3[j] - v2[j]) / (dx);
					}
					//corner
					else if (j == dy_steps - block_ysteps + 1) {
						vort << (v3[j] - v2[j]) / (dx) - (u2[j] - u2[j-1]) / (dy);
					}
					else if (j == 0) {
						vort << -(-3 * u2[j] + 4 * u2[j + 1] - u2[j + 2]) / (2 * dy);
					}
					else {
						vort << (v3[j] - v1[j]) / (2 * dx) - (u2[j + 1] - u2[j - 1]) / (2 * dy);
					}

					if (j < stop - 1) {
						vort << " ";
					}
					else {
						vort << "\n";
					}
				}
			}
			//other lines
			else {
				for (int a = 0; a < stop; a++) {
					u1[a] = u2[a];
					u2[a] = u3[a];
					v1[a] = v2[a];
					v2[a] = v3[a];

					ufile >> u3[a];
					vfile >> v3[a];
				}

				//calculate and write the vorticities
				for (int j = 0; j < stop; j++) {
					if (j == stop - 1) {
						vort << -(3 * u2[j] - 4 * u2[j - 1] + u2[j - 2]) / (2 * dy);
					}
					else if (j == 0) {
						vort << -(-3 * u2[j] + 4 * u2[j + 1] - u2[j + 2]) / (2 * dy);
					}
					else {
						vort << (v3[j] - v1[j]) / (2 * dx) - (u2[j + 1] - u2[j - 1]) / (2 * dy);
					}

					if (j < stop - 1) {
						vort << " ";
					}
					else {
						vort << "\n";
					}
				}
			}
		}
		ufile.close();
		vfile.close();
		vort.close();


	}


	/*
	Purpose: Uses current vorticies and velocities to step forward in time to new vorticities
	Complexity: O(xy), second-order central scheme for diffusion, and first-order upwind for advective terms
	Equation: dw/dt + duw/dx + dvw/dy = (nu)(d2w/dx2 + d2w/dy2)
	*/
	void solverB(double u, double v) {

	}


	/*
	Purpose: Find stream function values via vorticity
	Complexity: O(nxy), where n is the number of iterations until convergence
	Equation: Successive over-relaxtion iterative method (SOR)
	Notes: convergence criteria is (new+old)/(new) < 0.001; not necessary for timestep 0
	*/
	double* solverC(double dx, double dy, int dx_steps, int dy_steps, int block_xsteps, int block_ysteps, double * answer) {
		double w = 0.5;
		int count = 0;
		double threshold = 0.0001;

		//hard-coded
		const int size = 1100;
		double b[size] = {};
		//double coeff[900][900] = {};
		//double* coeff[900];
		//for (int i = 0; i < 900; i++)
		//	coeff[i] = (double*)malloc(900 * sizeof(double));

		double** coeff = new double* [size];
		for (int i = 0; i < size; i++)
			coeff[i] = new double[size];
		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++)
				coeff[i][j] = 0.;

		//vorticity data
		ifstream vort;
		vort.open("vort.txt");

		int gridpoints = dx_steps * dy_steps;
		for (int i = 0; i < gridpoints; i++) {
			int x = i/ dy_steps;
			int y = i % dy_steps;

			//inlet wall
			if (x == 0 && y <= dy_steps - block_ysteps) {				
				vort >> b[i];
				//hard-coded
				b[i] = -(9 / 25) * ((2 / 3) * pow(y, 3) - (1 / 2) * pow(y, 2));
				coeff[i][i] = 1;
			}
			//black wall
			else if(y == 0){
				vort >> b[i];
				b[i] = 0;
				coeff[i][i] = 1;
			}
			//red wall A
			else if (y == dy_steps - block_ysteps && x < block_xsteps) {
				vort >> b[i];
				//hard-coded
				b[i] = 0.015;
				coeff[i][i] = 1;
			}
			//red wall B
			else if (y == dy_steps && x >= block_xsteps - 1) {
				vort >> b[i];
				//hard-coded
				b[i] = 0.015;
				coeff[i][i] = 1;
			}
			//red wall cliff
			else if (x == block_xsteps - 1 && y >= dy_steps - block_ysteps) {
				vort >> b[i];
				//hard-coded
				b[i] = 0.015;
				coeff[i][i] = 1;
			}
			//outlet
			else if (x == dy_steps - 1) {
				vort >> b[i];
				//hard-coded
				b[i] = -(9/100)*((1/3)*pow(y,3) - (1/2)*pow(y,2));
				coeff[i][i] = 1;
			}
			//interior points
			else if(x >= block_xsteps - 1 || y <= dy_steps - block_ysteps){
				vort >> b[i];
				coeff[i][i] = -(1 / dx + 1 / dy);
				coeff[i][i + dy_steps] = -(1 / (2 * dx));
				coeff[i][i - dy_steps] = -(1 / (2 * dx));
				coeff[i][i + 1] = -(1 / (2 * dy));
				coeff[i][i - 1] = -(1 / (2 * dy));
			}
			//wall
			else {
				coeff[i][i] = 1;
			}
		}


		while (!converged(coeff, answer, b, gridpoints, threshold, count)) {
			count++;
			for (int i = 0; i < gridpoints; i++) {
				double sigma = 0;
				for (int j = 0; j < gridpoints; j++) {
					if (i != j) {
						sigma += coeff[i][j] * answer[j];
					}
				}
				answer[i] = (1 - w) * answer[i] + (w / coeff[i][i]) * (b[i] - sigma);
				//cout << "answer[" << i << "]: " << answer[i] << "\n";
			}
		}

		return answer;
	}

	//hard-coded
	bool converged(double** coeff, double answer[], double b[], int size, double threshold, int count) {
		double bigsum = 0;
		for (int i = 0; i < size; i++) {
			double sum = 0;
			for (int j = 0; j < size; j++) {
				sum += answer[j] * coeff[i][j];
				if (j == 999) {
					int fun = 0;
				}
			}
			bigsum += sum - b[i];
		}
		cout << "Iteration " << count << ": " << bigsum << "\n";
		return bigsum < threshold && bigsum > -1*(threshold);
	}


	/*
	Purpose: Find velocities with the current stream function
	Complexity: O(xy), simple first-order scheme for each gridpoint
	Equation: u = dPhi/dy, v = -dPhi/dx
	*/
	void solverDu(double phi) {

	}

	void solverDv(double phi) {

	}


	/*
	Purpose: Evaulate if the timesteps have converged to steady-state
	Complexity: O(xy), differencing every gridpoint
	Equation: (new+old)/(new) < 0.001
	*/
	bool isConverged(double w, double phi, double wm1, double phim1) {

		//compare between current iteration and last
		bool wConverged = true;
		bool phiConverged = true;
		
		if (wConverged && phiConverged) {
			return true;
		}
		else {
			return false;
		}
	}


	/*
	Purpose: Create initial velocities with some parameters
	Complexity: O(xy), for every gridpoint
	Equation: (new+old)/(new) < 0.001
	*/
	void initialConditions(double u_m, double aLength, double bLength, double aHeight, double bHeight, double dx, double dy) {
		double totalLength = aLength + bLength;
		double totalHeight = aHeight + bHeight;
		ofstream ufile;
		ofstream vfile;
		cout << "Writing initial conditions...\n";
		ufile.open("u.txt");
		vfile.open("v.txt");

		int dx_steps = int(round(totalLength / dx)) + 1;
		int dy_steps = int(round(totalHeight / dy)) + 1;
		int block_xsteps = int(round(aLength / dx)) + 1;
		int block_ysteps = int(round(aHeight / dy)) + 1;

		for (int i = 0; i < dx_steps; i++) {
			double x = dx * i;
			for (int j = 0; j < dy_steps; j++) {
				double y = dy * j;
				
				if (i < block_xsteps - 1) {
					if (j < block_ysteps) {
						ufile << init1(u_m, aHeight, y);
						vfile << 0;
					}
					else {
						ufile << "\n";
						vfile << "\n";
						break; //wall
					}					
				}
				else {
					vfile << 0;
					double yLimit = aHeight + (x - aLength) * (bHeight / bLength);
					if (y > yLimit) {
						ufile << 0;
					}
					else {
						ufile << init2(u_m, aHeight, y, yLimit);
					}
				}
				if (j < dy_steps - 1) {
					ufile << " ";
					vfile << " ";
				}
				else {
					ufile << "\n";
					vfile << "\n";
				}
			}
		}
		ufile.close();
		vfile.close();
	}

	double init1(double u_m, double aHeight, double y){
		return -(6 * u_m) * (y - aHeight) * (y) / (aHeight);
	}

	double init2(double u_m, double aHeight, double y, double yLimit) {
		return -(6 * u_m * aHeight/yLimit) * (y - yLimit) * (y) / (yLimit);
	}


	/*
	Purpose: Find dt_crit
	Complexity: O(1)
	Equation: see below
	*/
	double maxdt(double Re, double  dx, double  dy, double  maxu, double  maxv) {
		return 1 / ((maxu / dx) + (maxv / dy) + (2*Re)*(pow(dx,-2) + pow(dy, -2)));
	}
};

