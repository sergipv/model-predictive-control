#ifndef MPC_H
#define MPC_H

#include <cppad/cppad.hpp>
#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

typedef CPPAD_TESTVECTOR(double) Dvector;

const size_t N = 10;
const double dt = 0.1;
const size_t x_start = 0;
const size_t y_start = N;
const size_t psi_start = y_start + N;
const size_t v_start = psi_start + N;
const size_t cte_start = v_start + N;
const size_t epsi_start = cte_start + N;
const size_t delta_start = epsi_start + N;
const size_t a_start = delta_start + N - 1;

const double ref_v = 90;
const double Lf = 2.67;

// N number of timesteps. 6 elements in state. 2 actuators and N-1 actuations.
const size_t n_vars =  N * 6  + (N-1) * 2;
const size_t n_constraints = N * 6;


class FG_eval {

 public:
  Eigen::VectorXd coeffs;
  typedef CPPAD_TESTVECTOR(CppAD::AD<double>) ADvector;
    
  FG_eval(Eigen::VectorXd coeffs) {
    this->coeffs = coeffs;
  }
    
  void operator()(ADvector& fg, const ADvector& vars) {
		fg[0] = 0;
    updateCost(fg, vars);
		updateConstraints(fg, vars);
	}

 private:
  void updateCost(ADvector& fg, const ADvector& vars);

  void addStateErrorToCost(ADvector& fg, const ADvector& vars);

  void addActuatorUseToCost(ADvector& fg, const ADvector& vars);

  // Cost term added to minimize sudden changes on the actuation between two
  // timestamps. 
  void addSequentialActuatorToCost(ADvector& fg, const ADvector& vars);

  void updateConstraints(ADvector& fg, const ADvector& vars);

  void updatePrediction(ADvector& fg, const ADvector& vars, int t); 
};

class MPC {
 public:
  MPC() {}

  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuations.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);

 private:
  void addState(Dvector& v, const Eigen::VectorXd state);
  std::string getOptions();
  void initializeVars(Dvector& vars, const Eigen::VectorXd state);
  void setupVarBounds(Dvector& lower, Dvector& upper);
  void setupConstraintBounds(Dvector& lower, Dvector& upper, const Eigen::VectorXd state);
};

#endif /* MPC_H */
