#include <iostream>
#include "tools.h"

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
	VectorXd rmse(4);
	rmse << 0, 0, 0, 0;

	// Check for data parity
	if (estimations.size() != ground_truth.size() || estimations.size() == 0) {
		std::cout << "Invalid data" << std::endl;
		return rmse;
	}
	
	// Calculate sum of squared residuals
	for (unsigned int i = 0, i < estimations.size(), ++i) {
		VectorXd residual = estimations[i] - ground_truth[i];

		residual = residual.array() * residual.array();
		rmse += residual;
	}

	// Calculate mean
	rmse = rmse / estimations.size();

	// Calculate square root
	rmse = rmse.array().sqrt();

	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
	// State parameters
	float px = x_state[0];
	float py = x_state[1];
	float vx = x_state[2];
	float vy = x_state[3];

	// Temporary variables for calculations
	float px2 = pow(px, 2);
	float py2 = pow(py, 2);
	float px2_plus_py2 = px2 + py2;
	float px2py2_sqrt = sqrt(px2_plus_py2);
	float px2py2_1_5 = pow(px2_plus_py2, 1.5);

	// Create Jacobian matrix
	MatrixXd Hj = MatrixXd::Zero(3, 4);

	// Check for division by zero
	if (fabs(px2_plus_py2) < 0.00001) {
		std::cout << "CalculateJacobian() Error: Division by zero" << std::endl;
		return Hj;
	}

	Hj << px / px2py2_sqrt, 0, 0,
		-py / px2_plus_py2, px / px2_plus_py2, 0, 0,
		py*(vx*py - vy*px) / px2py2_1_5, px*(vy*px - vx*py) / px2py2_1_5, px / px2py2_sqrt, py / px2py2_sqrt;

	return Hj;
}
