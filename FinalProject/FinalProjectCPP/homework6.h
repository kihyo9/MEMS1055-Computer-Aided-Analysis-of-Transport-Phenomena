#pragma once
#include <fstream>
#include <cmath>
#include <iostream>
#include <chrono>
#include <ctime>    
using namespace std;
class homework6
{
public:

	void solve() {		
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
			dt = r * pow(dx, 2) / gamma;
		}
		const int dt_steps = int(ceil(t_last / dt)) + 1;

		//print simulation conditions
		printStability(dx, dt, u, gamma, r);

		//solve
		solveExplicit(dt, dx, gamma, u, L, t_0, dx_steps, dt_steps);

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

	void solveExplicit(double dt, double dx, double gamma, double u, double L, double t_0[],
				const int dx_steps, const int dt_steps) {

		string solverType = "Explicit";
		string fileName = "explicit-data.txt";

		ofstream myfile;
		cout << "Writing numbers file...\n";
		myfile.open(fileName);
		double xMax = dx_steps * dx;
		double tMax = dt_steps * dt;
		myfile << 0 << " " << xMax << " " << dx << "\n";
		myfile << 0 << " " << tMax << " " << dt << "\n";
		myfile.close();

		//Data
		double last_step[101];
		copy(t_0, last_step, dx_steps);
		int linesWritten = 0;

		//start clock
		auto start = std::chrono::system_clock::now();

		//start timestep iterations
		int timestepCount = 0;
		for (int i = 0; i < dt_steps; i++) {
			//output progress
			if (timestepCount % 1000 == 0) cout << "Timestep " << timestepCount << ", " << i*dt << "s" << "\n";
			timestepCount++;
				
			double t_step[101];
			if (i == 0) {
				copy(t_0, t_step, dx_steps);
			}
			else {
				for (int j = 0; j < dx_steps; j++) {
					t_step[j] = nextGridpointExplicit(last_step[(j - 1 + dx_steps) % dx_steps], last_step[j], last_step[(j + 1 + dx_steps) % dx_steps], dx, dt, u, gamma);
				}
			}

			bool newLine = writeToFile(fileName, t_step, i, linesWritten, dt_steps - 1, dt, dx_steps);
			if (newLine) {
				linesWritten++;
			}

			copy(t_step, last_step, dx_steps);

		}

		//end clock
		double timeElasped = performanceTime(start, solverType, gamma, dt_steps, u, L, dt);

		//result info
		cout << "That took " << timeElasped << "s\n";
		cout << "Timestep size: " << dt << "\n";
		cout << "Number of timesteps: " << dt_steps << "\n";
	}

	void solveUpwind(double dt, double dx, double gamma, double u, double L, double t_0[],
		const int dx_steps, const int dt_steps) {

		string solverType = "Upwind";
		string fileName = "upwind-data.txt";

		ofstream myfile;
		cout << "Writing numbers file...\n";
		myfile.open(fileName);
		double xMax = dx_steps * dx;
		double tMax = dt_steps * dt;
		myfile << 0 << " " << xMax << " " << dx << "\n";
		myfile << 0 << " " << tMax << " " << dt << "\n";
		myfile.close();

		//Data
		double last_step[101];
		copy(t_0, last_step, dx_steps);
		int linesWritten = 0;

		//start clock
		auto start = std::chrono::system_clock::now();

		//start timestep iterations
		int timestepCount = 0;
		for (int i = 0; i < dt_steps; i++) {
			//output progress
			if (timestepCount % 1000 == 0) cout << "Timestep " << timestepCount << ", " << i * dt << "s" << "\n";
			timestepCount++;

			double t_step[101];
			if (i == 0) {
				copy(t_0, t_step, dx_steps);
			}
			else {
				for (int j = 0; j < dx_steps; j++) {
					t_step[j] = nextGridpointUwpind(last_step[(j - 1 + dx_steps) % dx_steps], last_step[j], last_step[(j + 1 + dx_steps) % dx_steps], dx, dt, u, gamma);
				}
			}

			bool newLine = writeToFile(fileName, t_step, i, linesWritten, dt_steps - 1, dt, dx_steps);
			if (newLine) {
				linesWritten++;
			}

			copy(t_step, last_step, dx_steps);

		}

		//end clock
		double timeElasped = performanceTime(start, solverType, gamma, dt_steps, u, L, dt);

		//result info
		cout << "That took " << timeElasped << "s\n";
		cout << "Timestep size: " << dt << "\n";
		cout << "Number of timesteps: " << dt_steps << "\n";
	}

	double nextGridpointExplicit(double Tm1, double T0, double T1, double dx, double dt, double u, double gamma) {
		double coeffTm1 = u * dt / (2 * dx) + gamma * dt / (pow(dx, 2));
		double coeffT0 = 1 - 2 * gamma * dt / pow(dx,2);
		double coeffT1 = gamma * dt / pow(dx,2) - u * dt / (2 * dx);
		return coeffTm1 * Tm1 + coeffT0 * T0 + coeffT1 * T1;
	}

	double nextGridpointUwpind(double Tm1, double T0, double T1, double dx, double dt, double u, double gamma) {
		double coeffTm1 = u * dt / (dx) + gamma * dt / (pow(dx, 2));
		double coeffT0 = 1 - u * dt / dx - 2 * gamma * dt / pow(dx, 2);
		double coeffT1 = gamma * dt / pow(dx, 2);
		return coeffTm1 * Tm1 + coeffT0 * T0 + coeffT1 * T1;
	}

	int pyMod(int num, int mod) {
		if (num >= 0) {
			return num % mod;
		}
		else {
			return (num+mod) % mod;
		}
	}

	void printStability(double dx, double dt, double u, double gamma, double r) {
		cout << "Gamma = " << gamma << ", u = " << u << ", timestep size = " << dt << "\n";
		cout << "CFL is " << dt * u / dx << " (must be < " << pow((1 - 4 * pow(r, 2)), 0.5) << ")\n";
		cout << "r is " << r << " (must be < 0.5)\n";
		cout << "G is " << 2 * gamma / dt - pow(u, 2) << " (must be > 0)\n";
	}

	bool writeToFile(string fileName, double t_step[], int n, int linesWritten, int dt_steps, double dt, int dx_steps) {
		if (sequenceGenerator(linesWritten, dt_steps) == n) {
			ofstream myfile;
			myfile.open(fileName, std::ios_base::app);
			myfile << "Timestep_#" << n << " " << n * dt << "s" << "\n";
			for (int i = 0; i < dx_steps; i++) {
				if (i == dx_steps - 1) {
					myfile << t_step[i] << "\n";
				}
				else {
					myfile << t_step[i] << " ";
				}
			}
			myfile.close();
			return true;
		}
		else {
			return false;
		}
	}

	int sequenceGenerator(int i, int max) {
		if (i == 0) {
			return 0;
		}
		else if (i == 1) {
			return 1;
		}
		else {
			int ans = int(round(pow(2.16, i - 1)));
			if (ans < max) {
				return ans;
			}
			else {
				return max;
			}
		}
	}

	template<class time>
	double performanceTime(const time start, string solverType, double gamma, int dt_steps, double u, double L, double dt) {
		ofstream myfile;
		myfile.open("performance.txt", std::ios_base::app);
		auto finish = std::chrono::system_clock::now();
		std::time_t end_time = chrono::system_clock::to_time_t(finish);
		chrono::duration<double> elapsed = finish - start;

		myfile << "Start time: " << std::ctime(&end_time);
		myfile << solverType << " gamma = " << gamma << ", " << "timesteps = " << dt << ", " << "u = " << u << ", " << "L = " << L << "\n";
		myfile << elapsed.count() << "s\n\n";
		myfile.close();
		return elapsed.count();
	}

	template<class arrayType, class arrayType2>
	arrayType2* copy(arrayType copying[], arrayType2 target[], int dx_steps) {
		for (int i = 0; i < dx_steps; i++) {
			target[i] = copying[i];
		}

		return target;
	}
};

