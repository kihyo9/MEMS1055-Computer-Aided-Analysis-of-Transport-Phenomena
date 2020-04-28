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
	int mainSolver(double dx, double dy, double dt, int dx_steps, int dy_steps, int block_xsteps, int block_ysteps, double Re, string histfile) {
		int timestep = 0;
		ofstream hist;
		hist.open(histfile, std::ios_base::app);
		cout << "\nTimestep #" << timestep << ": " << dt * timestep << "s\n";
		hist << "\nTimestep #" << timestep << ": " << dt * timestep << "s\n";

		//calculate initial vorticies
		solverA(dx, dy, dx_steps, dy_steps, block_xsteps, block_ysteps);

		//calculate the initial stream function
		cout << "Writing initial stream function..." << "\n";
		hist << "Writing initial stream function..." << "\n";
		hist.close();
		int size = dx_steps * dy_steps;
		double* answer = new double[size]{};
		solverC(dx, dy, dx_steps, dy_steps, block_xsteps, block_ysteps, answer, histfile); //get vorticities for timestep 0

		/*
		Start of timestep iterations
		*/
		while (!FileAndConvergenceCheck(dt, timestep, dx_steps, dy_steps, block_xsteps, block_ysteps, histfile)) {
			timestep++;
			hist.open(histfile, std::ios_base::app);
			cout << "\nTimestep #" << timestep << ": " << dt * timestep << "s\n";
			hist << "\nTimestep #" << timestep << ": " << dt * timestep << "s\n";
			hist.close();

			//new vorticies
			solverB(dt, dx, dy, dx_steps, dy_steps, block_xsteps, block_ysteps, Re);

			//new stream function
			answer = new double[size] {};
			solverC(dx, dy, dx_steps, dy_steps, block_xsteps, block_ysteps, answer, histfile);

			//new velocities
			solverD(dx, dy, dx_steps, dy_steps, block_xsteps, block_ysteps);
		}
		/*
		End of timestep iterations
		*/

		return timestep;
	}

	/*
	Purpose: Compare files, if they exist, for convergence and do some file clean-up
	Notes: Throws if files do not exist; returns true if converged, false otherwise
	Criterion: (new+old)/(new) < threshold
	*/
	bool FileAndConvergenceCheck(double dt, int& timestep, int dx_steps, int dy_steps, int block_xsteps, int block_ysteps, string histfile, double threshold = 0.0001) {
		if (timestep != 0) {
			//check if files exist
			if (!(exists_test0("u2.txt") && exists_test0("v2.txt") && exists_test0("stream2.txt") && exists_test0("vort2.txt")))
				throw std::exception("Missing new file.");

			if (!(exists_test0("u.txt") && exists_test0("v.txt") && exists_test0("stream.txt") && exists_test0("vort.txt")))
				throw std::exception("Missing base file.");

			//compare
			double uchange = compare("u.txt", "u2.txt", dx_steps, dy_steps, block_xsteps, block_ysteps);
			double vchange = compare("v.txt", "v2.txt", dx_steps, dy_steps, block_xsteps, block_ysteps);
			double streamchange = compare("stream.txt", "stream2.txt", dx_steps, dy_steps, block_xsteps, block_ysteps);
			double vortchange = compare("vort.txt", "vort2.txt", dx_steps, dy_steps, block_xsteps, block_ysteps);

			if (uchange < threshold && vchange < threshold && streamchange < threshold && vortchange < threshold)
				return true;

			//delete old
			if (exists_test0("u.txt")) {
				if (remove("u.txt") != 0)
					throw std::exception("Base file removal failure.");
			}
			if (exists_test0("v.txt")) {
				if (remove("v.txt") != 0)
					throw std::exception("Base file removal failure.");
			}
			if (exists_test0("stream.txt")) {
				if (remove("stream.txt") != 0)
					throw std::exception("Base file removal failure.");
			}
			if (exists_test0("vort.txt")) {
				if (remove("vort.txt") != 0)
					throw std::exception("Base file removal failure.");
			}

			//rename
			if (rename("u2.txt", "u.txt") != 0)
				throw std::exception("Renaming failure.");
			if (rename("v2.txt", "v.txt") != 0)
				throw std::exception("Renaming failure.");
			if (rename("stream2.txt", "stream.txt") != 0)
				throw std::exception("Renaming failure.");
			if (rename("vort2.txt", "vort.txt") != 0)
				throw std::exception("Renaming failure.");

			ofstream hist;
			hist.open(histfile, std::ios_base::app);
			ofstream info;
			info.open("info-" + histfile, std::ios_base::app);

			cout << "u-change: " << uchange << ", v-change: " << vchange << ", streamfunc-change: " << streamchange << ", vorticity-change: " << vortchange << "\n";
			hist << "u-change: " << uchange << ", v-change: " << vchange << ", streamfunc-change: " << streamchange << ", vorticity-change: " << vortchange << "\n";
			info << uchange << " " << vchange << " " << streamchange << " " << vortchange << "\n";

			hist.close();
			info.close();

			return false;
		}
		else {
			//check if files exist
			if (!(exists_test0("u2.txt") && exists_test0("v2.txt") && exists_test0("stream2.txt") && exists_test0("vort2.txt")))
				throw std::exception("Missing new file.");

			//delete old
			if (exists_test0("u.txt")) {
				if (remove("u.txt") != 0)
					throw std::exception("Base file removal failure.");
			}
			if (exists_test0("v.txt")) {
				if (remove("v.txt") != 0)
					throw std::exception("Base file removal failure.");
			}
			if (exists_test0("stream.txt")) {
				if (remove("stream.txt") != 0)
					throw std::exception("Base file removal failure.");
			}
			if (exists_test0("vort.txt")) {
				if (remove("vort.txt") != 0)
					throw std::exception("Base file removal failure.");
			}

			//rename
			if (rename("u2.txt", "u.txt") != 0)
				throw std::exception("Renaming failure.");
			if (rename("v2.txt", "v.txt") != 0)
				throw std::exception("Renaming failure.");
			if (rename("stream2.txt", "stream.txt") != 0)
				throw std::exception("Renaming failure.");
			if (rename("vort2.txt", "vort.txt") != 0)
				throw std::exception("Renaming failure.");

			return false;
		}
	}

	inline bool exists_test0(const std::string& name) {
		ifstream f(name.c_str());
		return f.good();
	}

	//hard-coded
	double compare(string old_file, string new_file, int dx_steps = 91, int dy_steps = 11, int block_xsteps = 16, int block_ysteps = 6) {
		int size = dx_steps * dy_steps;

		ifstream oldf;
		oldf.open(old_file);

		ifstream newf;
		newf.open(new_file);

		double diffsum = 0;
		for (int i = 0; i < size; i++) {
			int x = i / dy_steps;
			int y = i % dy_steps;

			if ((x >= block_xsteps - 1 || y <= dy_steps - block_ysteps) && (x < dx_steps)) {
				double old_value;
				oldf >> old_value;

				double new_value;
				newf >> new_value;

				if (new_value != 0.) {
					diffsum += abs((old_value - new_value)/new_value);
				}
				else {
					diffsum += abs(old_value - new_value);
				}
				
			}
		}
		oldf.close();
		newf.close();

		return diffsum / size;
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
		vort.open("vort2.txt");

		double* u1 = new double[dy_steps]{};
		double* u2 = new double[dy_steps] {};
		double* u3 = new double[dy_steps] {};

		double* v1 = new double[dy_steps] {};
		double* v2 = new double[dy_steps] {};
		double* v3 = new double[dy_steps] {};
		

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
	Purpose: Uses current vorticies, strea, function and velocities to step forward in time to new vorticities
	Complexity: O(xy), second-order central scheme for diffusion, and first-order upwind for advective terms
	Equation: dw/dt + duw/dx + dvw/dy = (nu)(d2w/dx2 + d2w/dy2)
	*/
	void solverB(double dt, double dx, double dy, int dx_steps, int dy_steps, int block_xsteps, int block_ysteps, double Re) {
		int size = dx_steps * dy_steps;

		ifstream ufile;
		ufile.open("u.txt");
		double** udata = new double* [dx_steps];
		for (int i = 0; i < dx_steps; i++)
			udata[i] = new double[dy_steps] {};		

		ifstream vfile;
		vfile.open("v.txt");
		double** vdata = new double* [dx_steps];
		for (int i = 0; i < dx_steps; i++)
			vdata[i] = new double[dy_steps] {};

		ifstream vortfile;
		vortfile.open("vort.txt");
		double** vortdata = new double* [dx_steps];
		for (int i = 0; i < dx_steps; i++)
			vortdata[i] = new double[dy_steps] {};

		ifstream streamfile;
		streamfile.open("stream.txt");
		double** streamdata = new double* [dx_steps];
		for (int i = 0; i < dx_steps; i++)
			streamdata[i] = new double[dy_steps] {};

		for (int i = 0; i < size; i++) {
			int x = i / dy_steps;
			int y = i % dy_steps;
			int stop = x < block_xsteps - 1 ? block_ysteps : dy_steps;
			if ((x >= block_xsteps - 1 || y <= dy_steps - block_ysteps) && (x < dx_steps)) {
				ufile >> udata[x][y];
				vfile >> vdata[x][y];
				vortfile >> vortdata[x][y];
				streamfile >> streamdata[x][y];
			}
		}
		ufile.close();
		vfile.close();
		vortfile.close();
		streamfile.close();
		
		ofstream new_w;
		new_w.open("vort2.txt");
		//write out coefficients
		int gridpoints = dx_steps * dy_steps;
		for (int i = 0; i < gridpoints; i++) {
			int x = i / dy_steps;
			int y = i % dy_steps;

			//inlet wall
			if (x == 0 && y < dy_steps - block_ysteps && y > 0) {
				new_w << -2 * (streamdata[x+1][y] - streamdata[x][y]) / pow(dx, 2);
			}
			//black wall
			else if (y == 0) {
				new_w << -2 * (streamdata[x][y+1] - streamdata[x][y]) / pow(dx, 2);
			}
			//red wall A
			else if (y == dy_steps - block_ysteps && x < block_xsteps) {
				new_w << -2 * (streamdata[x][y-1] - streamdata[x][y]) / pow(dx, 2);
			}
			//red wall B
			else if (y == dy_steps - 1 && x >= block_xsteps - 1) {
				new_w << -2 * (streamdata[x][y-1] - streamdata[x][y]) / pow(dx, 2);
			}
			//red wall cliff
			else if (x == block_xsteps - 1 && y >= dy_steps - block_ysteps) {
				new_w << -2 * (streamdata[x + 1][y] - streamdata[x][y]) / pow(dx, 2);
			}
			//outlet
			else if (x == dx_steps - 1) {
				new_w << -2 * (streamdata[x-1][y] - streamdata[x][y]) / pow(dx, 2);
			}
			//interior points
			else if (x >= block_xsteps - 1 || y <= dy_steps - block_ysteps) {
				double term1 = udata[x][y] >= 0 ? (udata[x+1][y] * vortdata[x+1][y] - udata[x][y] * vortdata[x][y]) / dx \
					: (udata[x][y] * vortdata[x][y] - udata[x-1][y] * vortdata[x-1][y]/dx);

				double term2 = vdata[x][y] >= 0 ? (vdata[x][y+1] * vortdata[x][y+1] - vdata[x][y] * vortdata[x][y]) / dy \
					: (vdata[x][y] * vortdata[x][y] - vdata[x][y-1] * vortdata[x][y-1] / dy);

				double term3 = (vortdata[x + 1][y] - 2 * vortdata[x][y] + vortdata[x - 1][y]) / pow(dx,2);

				double term4 = (vortdata[x][y + 1] - 2 * vortdata[x][y] + vortdata[x][y - 1]) / pow(dy,2);

				new_w << vortdata[x][y] + dt * (-term1 - term2 + term3 + term4);
			}
			//wall
			else {
			}

			int stop = x < block_xsteps - 1 ? block_ysteps : dy_steps;
			if ((x >= block_xsteps - 1 || y <= dy_steps - block_ysteps) && (x < dx_steps)) {
				if (y == stop - 1) {
					new_w << "\n";
				}
				else {
					new_w << " ";
				}
			}
		}
		new_w.close();
	}


	/*
	Purpose: Find stream function values via vorticity
	Complexity: O(nxy), where n is the number of iterations until convergence
	Equation: Successive over-relaxtion iterative method (SOR)
	Notes: convergence criteria is (new+old)/(new) < 0.001; not necessary for timestep 0
	*/
	void solverC(double dx, double dy, int dx_steps, int dy_steps, int block_xsteps, int block_ysteps, double * answer, string histfile) {
		int size = dx_steps * dy_steps;
		//const int size = 1100;
		double* b = new double[size] {};
		//double coeff[900][900] = {};
		//double* coeff[900];
		//for (int i = 0; i < 900; i++)
		//	coeff[i] = (double*)malloc(900 * sizeof(double));

		double** coeff = new double* [size];
		for (int i = 0; i < size; i++)
			coeff[i] = new double[size] {};
		//for (int i = 0; i < size; i++)
		//	for (int j = 0; j < size; j++)
		//		coeff[i][j] = 0.;

		//vorticity data
		ifstream vort;
		vort.open("vort.txt");

		//write out coefficients
		int gridpoints = dx_steps * dy_steps;
		for (int i = 0; i < gridpoints; i++) {
			int x = i/ dy_steps;
			int y = i % dy_steps;

			//inlet wall
			if (x == 0 && y <= dy_steps - block_ysteps) {				
				vort >> b[i];
				//hard-coded - affected by pipe geometry and mean intial velocity
				b[i] = -(9. / 25.) * ((2. / 3.) * pow(y*dy, 3) - (1. / 2.) * pow(y*dy, 2));
				coeff[i][i] = 1.;
			}
			//black wall
			else if(y == 0){

				vort >> b[i];
				b[i] = 0;
				coeff[i][i] = 1.;
			}
			//red wall A
			else if (y == dy_steps - block_ysteps && x < block_xsteps) {
				vort >> b[i];
				//hard-coded - affected by pipe geometry and mean intial velocity
				b[i] = 0.015;
				coeff[i][i] = 1.;
			}
			//red wall B
			else if (y == dy_steps-1 && x >= block_xsteps - 1) {
				vort >> b[i];
				//hard-coded - affected by pipe geometry and mean intial velocity
				b[i] = 0.015;
				coeff[i][i] = 1.;
			}
			//red wall cliff
			else if (x == block_xsteps - 1 && y >= dy_steps - block_ysteps) {
				vort >> b[i];
				//hard-coded - affected by pipe geometry and mean intial velocity
				b[i] = 0.015;
				coeff[i][i] = 1.;
			}
			//outlet
			else if (x == dx_steps - 1) {
				vort >> b[i];
				//hard-coded - affected by pipe geometry and mean intial velocity
				b[i] = -(9./100.)*((1./3.)*pow(y * dy,3) - (1./2.)*pow(y * dy,2));
				coeff[i][i] = 1;
			}
			//interior points
			else if(x >= block_xsteps - 1 || y <= dy_steps - block_ysteps){
				vort >> b[i];
				double dx2 = pow(dx, 2);
				double dy2 = pow(dy, 2);
				coeff[i][i] = (2. / dx2 + 2. / dy2);
				coeff[i][i + dy_steps] = -(1. / (dx2));
				coeff[i][i - dy_steps] = -(1. / (dx2));
				coeff[i][i + 1] = -(1. / (dy2));
				coeff[i][i - 1] = -(1. / (dy2));

			}
			//wall
			else {
				coeff[i][i] = 1.;
			}
		}
		vort.close();

		//SOR
		double w = 1.3;
		int count = 0;
		double threshold = 0.0001 / size;
		
		while (!converged(coeff, answer, b, gridpoints, threshold)) {
			count++;
			for (int i = 0; i < gridpoints; i++) {
				double sigma = 0;
				for (int j = 0; j < gridpoints; j++) {
					if (i != j) {
						sigma += coeff[i][j] * answer[j];
					}
				}
				answer[i] = (1 - w) * answer[i] + (w / coeff[i][i]) * (b[i] - sigma);
			}
		}
		ofstream hist;
		hist.open(histfile, std::ios_base::app);
		hist << "SOR iterations: " << count << "\n";
		cout << "SOR iterations: " << count << "\n";
		hist.close();

		//write to file
		ofstream stream;
		stream.open("stream2.txt");

		for (int i = 0; i < size; i++) {
			int x = i / dy_steps;
			int y = i % dy_steps;
			int stop = x < block_xsteps - 1 ? block_ysteps : dy_steps;
			if ((x >= block_xsteps - 1 || y <= dy_steps - block_ysteps) && (x < dx_steps)) {
				stream << answer[i];
				if (y == stop - 1) {
					stream << "\n";
				}
				else {
					stream << " ";
				}
			}
		}
		stream.close();
	}

	bool converged(double** coeff, double answer[], double b[], int size, double threshold) {
		double bigsum = 0;
		for (int i = 0; i < size; i++) {
			double sum = 0;
			for (int j = 0; j < size; j++) {
				sum += answer[j] * coeff[i][j];
			}
			bigsum += sum - b[i];
		}
		return bigsum < threshold && bigsum > -1*(threshold);
	}


	/*
	Purpose: Find velocities with the current stream function
	Complexity: O(xy), simple first-order scheme for each gridpoint
	Equation: u = dPhi/dy, v = -dPhi/dx
	*/
	void solverD(double dx, double dy, int dx_steps, int dy_steps, int block_xsteps, int block_ysteps) {
		int size = dx_steps * dy_steps;
		ifstream streamfile;
		streamfile.open("stream.txt");
		double** streamdata = new double* [dx_steps];
		for (int i = 0; i < dx_steps; i++)
			streamdata[i] = new double[dy_steps] {};
		for (int i = 0; i < size; i++) {
			int x = i / dy_steps;
			int y = i % dy_steps;
			int stop = x < block_xsteps - 1 ? block_ysteps : dy_steps;
			if ((x >= block_xsteps - 1 || y <= dy_steps - block_ysteps) && (x < dx_steps)) {
				streamfile >> streamdata[x][y];
			}
		}
		streamfile.close();

		ofstream ufile;
		ufile.open("u2.txt");
		ofstream vfile;
		vfile.open("v2.txt");

		int gridpoints = dx_steps * dy_steps;
		for (int i = 0; i < gridpoints; i++) {
			int x = i / dy_steps;
			int y = i % dy_steps;

			//inlet
			if (x == 0 && y < dy_steps - block_ysteps && y > 0) {
				//hard-coded - maybe not anymore?
				ufile << init1(y*dy);
				vfile << 0;
			}
			//black wall
			else if (y == 0) {
				ufile << 0;
				vfile << 0;
			}
			//red wall A
			else if (y == dy_steps - block_ysteps && x < block_xsteps) {
				ufile << 0;
				vfile << 0;
			}
			//red wall B
			else if (y == dy_steps - 1 && x >= block_xsteps - 1) {
				ufile << 0;
				vfile << 0;
			}
			//red wall cliff
			else if (x == block_xsteps - 1 && y >= dy_steps - block_ysteps) {
				ufile << 0;
				vfile << 0;
			}
			//outlet
			else if (x == dx_steps - 1) {
				//hard-coded - maybe not anymore?
				ufile << init2(y*dy);
				vfile << 0;
			}
			//interior points
			else if (x >= block_xsteps - 1 || y <= dy_steps - block_ysteps) {
				ufile << (streamdata[x][y+1] - streamdata[x][y-1]) / (2*dy);
				vfile << (streamdata[x+1][y] - streamdata[x-1][y]) / (2*dx);
			}
			//wall
			else {
			}

			int stop = x < block_xsteps - 1 ? block_ysteps : dy_steps;
			if ((x >= block_xsteps - 1 || y <= dy_steps - block_ysteps) && (x < dx_steps)) {
				if (y == stop - 1) {
					ufile << "\n";
					vfile << "\n";
				}
				else {
					ufile << " ";
					vfile << " ";
				}
			}
		}

		ufile.close();
		vfile.close();
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
		ufile.open("u2.txt");
		vfile.open("v2.txt");

		int dx_steps = int(round(totalLength / dx)) + 1;
		int dy_steps = int(round(totalHeight / dy)) + 1;
		int block_xsteps = int(round(aLength / dx)) + 1;
		int block_ysteps = int(round(aHeight / dy)) + 1;

		for (int i = 0; i < dx_steps; i++) {
			double x = dx * i;
			for (int j = 0; j < dy_steps; j++) {
				double y = dy * j;
				
				if (i < block_xsteps - 1) {
					if (j < block_ysteps) { //region A
						ufile << init1(y, u_m, aHeight);
						vfile << 0;
					}
					else {
						ufile << "\n";
						vfile << "\n";
						break; //wall
					}					
				}
				else { //region B
					vfile << 0;
					double yLimit = aHeight + (x - aLength) * (bHeight / bLength);
					if (y > yLimit) { //red wall B
						ufile << 0;
					}
					else {
						ufile << init2(y, u_m, aHeight, yLimit);
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

	double init1(double y, double u_m = 0.03, double aHeight = 0.5){
		return -(6 * u_m) * (y - aHeight) * (y) / pow(aHeight, 2);
	}

	double init2(double y, double u_m = 0.03, double aHeight = 0.5, double yLimit = 1) {
		return -(6 * u_m * aHeight/yLimit) * (y - yLimit) * (y) / pow(yLimit,2);
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

