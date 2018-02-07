#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include <numeric>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;

  Dvector vars(n_vars);
  initializeVars(vars, state);
	for (size_t i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  Dvector vars_lower(n_vars);
  Dvector vars_upper(n_vars);

  setupVarBounds(vars_lower, vars_upper);

	Dvector constraints_lower(n_constraints);
  Dvector constraints_upper(n_constraints);

  setupConstraintBounds(constraints_lower, constraints_upper, state);
    
  FG_eval fg_eval(coeffs);

  std::string options = getOptions();
  
  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lower, vars_upper, constraints_lower,
      constraints_upper, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

	vector<double> result;
	result.push_back(solution.x[delta_start]);
	result.push_back(solution.x[a_start]);

	for (size_t i=0; i<N-1; i++) {
		result.push_back(solution.x[x_start+i+1]);
		result.push_back(solution.x[y_start+i+1]);
	}

	return result;
}

void MPC::initializeVars(Dvector& vars, const Eigen::VectorXd state) {
  for (size_t i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }
}

void MPC::setupVarBounds(Dvector& lower, Dvector& upper) {
  for (size_t i=0; i<delta_start; i++) {
    lower[i] = -1.0e19;// numeric_limits<double>::min();
    upper[i] = 1.0e19;// numeric_limits<double>::max();
  }
  for (size_t i=delta_start; i<a_start; i++) {
    lower[i] = -0.436332;
    upper[i] = 0.436332;
  }
  for (size_t i=a_start; i<n_vars; i++) {
    lower[i] = -1.0;
    upper[i] = 1.0;
  }
}

void MPC::setupConstraintBounds(Dvector& lower, Dvector& upper, const Eigen::VectorXd state) {
  for (size_t i=0; i<n_constraints; i++) {
    lower[i] = 0;
    upper[i] = 0;
  }
  addState(lower, state);
  addState(upper, state);
}

void MPC::addState(Dvector& v, const Eigen::VectorXd state) {
  v[x_start] = state[0];
  v[y_start] = state[1];
  v[psi_start] = state[2];
  v[v_start] = state[3];
  v[cte_start] = state[4];
  v[epsi_start] = state[5];

}

std::string MPC::getOptions() {
  std::string options;

  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  return options;
}

// FG_eval implementation

void FG_eval::updateCost(ADvector& fg, const ADvector& vars) {
  fg[0] = 0;
  addStateErrorToCost(fg, vars);
  addActuatorUseToCost(fg, vars);
  addSequentialActuatorToCost(fg, vars);
}

void FG_eval::addStateErrorToCost(ADvector& fg, const ADvector& vars) {
  for (size_t t=0; t<N; t++) {
    fg[0] += 2000 * CppAD::pow(vars[cte_start+t], 2);
    fg[0] += 2000 * CppAD::pow(vars[epsi_start+t], 2);
    fg[0] += CppAD::pow(vars[v_start+t] - ref_v, 2);
  }
}

void FG_eval::addActuatorUseToCost(ADvector& fg, const ADvector& vars) {
  for (size_t t=0; t<N-1; t++) {
    fg[0] += 5 * CppAD::pow(vars[delta_start+t], 2); 
    fg[0] += 5 * CppAD::pow(vars[a_start+t], 2);
  }
}

void FG_eval::addSequentialActuatorToCost(ADvector& fg, const ADvector& vars) {
  for (size_t t=0; t<N-2; t++) {
    fg[0] += 2000 * CppAD::pow(vars[delta_start+t+1] - vars[delta_start+t], 2);
    fg[0] += 10 * CppAD::pow(vars[a_start+t+1] - vars[a_start+t], 2);
  }
}

void FG_eval::updateConstraints(ADvector& fg, const ADvector& vars) {
	fg[1+x_start] = vars[x_start];
  fg[1+y_start] = vars[y_start];
  fg[1+psi_start] = vars[psi_start];
  fg[1+v_start] = vars[v_start];
  fg[1+cte_start] = vars[cte_start];
  fg[1+epsi_start] = vars[epsi_start];

	for (size_t t=1; t<N; t++) {
  	updatePrediction(fg, vars, t);
  }
}

void FG_eval::updatePrediction(ADvector& fg, const ADvector& vars, int t) {
	CppAD::AD<double> x1 = vars[x_start + t];
  CppAD::AD<double> y1 = vars[y_start + t];
  CppAD::AD<double> psi1 = vars[psi_start + t];
  CppAD::AD<double> v1 = vars[v_start + t];
  CppAD::AD<double> cte1 = vars[cte_start + t];
  CppAD::AD<double> epsi1 = vars[epsi_start + t];

  // The state at time t.
  CppAD::AD<double> x0 = vars[x_start + t - 1];
  CppAD::AD<double> y0 = vars[y_start + t - 1];
  CppAD::AD<double> psi0 = vars[psi_start + t - 1];
  CppAD::AD<double> v0 = vars[v_start + t - 1];
  CppAD::AD<double> cte0 = vars[cte_start + t - 1];
  CppAD::AD<double> epsi0 = vars[epsi_start + t - 1];

  // Only consider the actuation at time t.
  CppAD::AD<double> delta0 = vars[delta_start + t - 1];
  CppAD::AD<double> a0 = vars[a_start + t - 1];

  CppAD::AD<double> f0 = coeffs[0] + coeffs[1]*x0 + coeffs[2]*x0*x0 + coeffs[3]*x0*x0*x0;
  CppAD::AD<double> psides0 = CppAD::atan(coeffs[1] + 2*coeffs[2]*x0 + 3*coeffs[3]*x0*x0);

	fg[1+x_start+t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
  fg[1+y_start+t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
  fg[1+psi_start+t] = psi1 - (psi0 - v0 * delta0 / Lf * dt);
  fg[1+v_start+t] = v1 - (v0 + a0 * dt);
  fg[1+cte_start+t] =
          cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
  fg[1+epsi_start+t] =
          epsi1 - ((psi0 - psides0) - v0 * delta0 / Lf * dt);
}

