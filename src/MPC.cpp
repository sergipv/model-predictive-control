#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include <numeric>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  size_t i;

  Dvector vars(n_vars);
  initializeVars(vars, state);

  Dvector vars_lower(n_vars);
  Dvector vars_upper(n_vars);

  setupVarBounds(vars_lower, vars_upper);

  Dvector constraints_lower(n_constraints);
  Dvector constraints_upper(n_constraints);

  setupConstraintBounds(constraints_lower, constraints_upper, state);

  FG_eval fg_eval(coeffs);

  std::string options;
  // Uncomment this if you'd like more print information
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
  addState(vars, state);
}

void MPC::addState(Dvector& v, const Eigen::VectorXd state) {
  v[x_start] = state[0];
  v[y_start] = state[1];
  v[psi_start] = state[2];
  v[v_start] = state[3];
  v[cte_start] = state[4];
  v[epsi_start] = state[5];

}

void MPC::setupVarBounds(Dvector& lower, Dvector& upper) {
  for (size_t i=0; i<delta_start; i++) {
    lower[i] = numeric_limits<double>::min();
    upper[i] = numeric_limits<double>::max();
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
