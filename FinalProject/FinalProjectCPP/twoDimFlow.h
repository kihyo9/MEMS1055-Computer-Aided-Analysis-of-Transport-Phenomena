#pragma once
class twoDimFlow
{
public:
	void mainSolver(double uInit, double vInit, double dx, double dy, double dt) {
		//quantities of interest
		double u;
		double v;
		double w;
		double phi;

		//calculate initial vorticies
		solverA(w);

		//retained data
		double wm1;
		double phim1;

		/*
		Start of timestep iterations
		*/		
		while (!isConverged(w,phi,wm1,phim1)) {
			//new vorticies
			solverB(u, v);

			//new stream function
			solverC(w);

			//new velocities
			solverDu(phi);
			solverDv(phi);
		}
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
	void solverA(double w) {

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
	void solverC(double w) {

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
		bool wConverged;
		bool phiConverged;
		
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
	double initialConditions(double someParameters) {
		return 0.0;
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

