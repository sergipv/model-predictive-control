#ifndef MPC_H
#define MPC_H

#include <cppad/cppad.hpp>
#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

typedef CPPAD_TESTVECTOR(double) Dvector;

const size_t N = 20;
const double dt = 0.05;
const size_t x_start = 0;
const size_t y_start = N;
const size_t psi_start = y_start + N;
const size_t v_start = psi_start + N;
const size_t cte_start = v_start + N;
const size_t epsi_start = cte_start + N;
const size_t delta_start = epsi_start + N;
const size_t a_start = delta_start + N - 1;

const double ref_v = 70;
const double Lf = 2.67;

// N number of timesteps. 6 elements in state. 2 actuators and N-1 actuations.
const size_t n_vars =  N * 6  + (N-1) * 2;
const size_t n_constraints = N * 6;

class MPC {
 private:

  class FG_eval {
   private:
    double ref_v;

   public:
    Eigen::VectorXd coeffs;
    typedef CPPAD_TESTVECTOR(CppAD::AD<double>) ADvector;
    
    FG_eval(Eigen::VectorXd coeffs) {
      this->coeffs = coeffs;
    }
    
    void operator()(ADvector& fg, const ADvector& vars) {
      updateCost(fg, vars);
      updateConstraints(fg, vars);
    }

   private:
    void updateCost(ADvector& fg, const ADvector& vars) {
      // The cost is stored in fg[0]
      fg[0] = 0;
      
      addStateErrorToCost(fg, vars);
      addActuatorUseToCost(fg, vars);
      addSequentialActuatorToCost(fg, vars);
    }

    void addStateErrorToCost(ADvector& fg, const ADvector& vars) {
      for (size_t t=0; t<N; t++) {
        fg[0] += CppAD::pow(vars[cte_start+t], 2);
        fg[0] += CppAD::pow(vars[epsi_start+t], 2);
        fg[0] += CppAD::pow(vars[v_start+t] - ref_v, 2);
      }
    }

    void addActuatorUseToCost(ADvector& fg, const ADvector& vars) {
      for (size_t t=0; t<N-1; t++) {
        fg[0] += CppAD::pow(vars[delta_start+t], 2);
        fg[0] += CppAD::pow(vars[a_start+t], 2);
      }
    }

    // Cost term added to minimize sudden changes on the actuation between two
    // timestamps. 
    void addSequentialActuatorToCost(ADvector& fg, const ADvector& vars) {
      for (size_t t=0; t<N-1; t++) {
        fg[0] += CppAD::pow(vars[delta_start+t+1] - vars[delta_start+t], 2);
        fg[0] += CppAD::pow(vars[a_start+t+1] - vars[a_start+t], 2);
      }
    }

    void updateConstraints(ADvector& fg, const ADvector& vars) {
      fg[1 + x_start] = vars[x_start];
      fg[1 + y_start] = vars[y_start];
      fg[1 + psi_start] = vars[psi_start];
      fg[1 + v_start] = vars[v_start];
      fg[1 + cte_start] = vars[cte_start];
      fg[1 + epsi_start] = vars[epsi_start];

      for (size_t t=1; t<N; t++) {
        updatePrediction(fg, vars, t);
      }
    }

    void updatePrediction(ADvector& fg, const ADvector& vars, int t) {
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

      CppAD::AD<double> f0 = coeffs[0] + coeffs[1] * x0;
      CppAD::AD<double> psides0 = CppAD::atan(coeffs[1]);

			fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[1 + psi_start + t] = psi1 - (psi0 + v0 * delta0 / Lf * dt);
      fg[1 + v_start + t] = v1 - (v0 + a0 * dt);
      fg[1 + cte_start + t] =
          cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
      fg[1 + epsi_start + t] =
          epsi1 - ((psi0 - psides0) + v0 * delta0 / Lf * dt);
    }
  };

 public:
  MPC() {}

  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuations.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);

 private:
  void addState(Dvector& v, const Eigen::VectorXd state);
  void initializeVars(Dvector& vars, const Eigen::VectorXd state);
  void setupVarBounds(Dvector& lower, Dvector& upper);
  void setupConstraintBounds(Dvector& lower, Dvector& upper, const Eigen::VectorXd state);
};

struct State {
  double x;
  double y;
  double psi;
  double v;
  double cte;
  double epsi;
};


void initialState(State& s, const Eigen::VectorXd v);
void updateState(State& s, double dt);

#endif /* MPC_H */
