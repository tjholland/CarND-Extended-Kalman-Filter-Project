#include "kalman_filter.h"

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
	MatrixXd F_t = F_.transpose();

	x_ = F_ * x_;

	P_ = F_ * P_ * F_t + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
	VectorXd z_p = H_ * x_;
	VectorXd y = z - z_p;

	MatrixXd H_t = H_.transpose();
	MatrixXd PH_t = P_ * H_t;
	MatrixXd S = H_ * PH_t + R_;
	MatrixXd S_i = S.inverse();
	MatrixXd K = PH_t * S_i;

	x_ = x_ + K * y;
	long xn = x_.size();
	MatrixXd I = MatrixXd::Identity(xn, xn);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
	float rho_p = sqrt(pow(x_[0], 2) + pow(x_[1], 2));
	float phi_p = 0.0;
	if (fabs(x_[0]) > 0.001) {
		phi_p = atan2(x_[1], x_[0]);
	}

	float rhodot_p = 0.0;
	if (fabs(rho_p) > 0.001) {
		rhodot_p = (x_[0] * x_[2] + x_[1] * x_[3]) / rho_p;
	}

	VectorXd z_p(3);
	z_p << rho_p, phi_p, rhodot_p;

	VectorXd y = z - z_p;

	MatrixXd H_t = H_.transpose();
	MatrixXd PH_t = P_ * H_t;
	MatrixXd S = H_ * PH_t + R_;
	MatrixXd S_i = S.inverse();
	MatrixXd K = PH_t * S_i;

	x_ = x_ + (K * y);
	long xn = x_.size();
	MatrixXd I = MatrixXd::Identity(xn, xn);
	P_ = (I - K * H_) * P_;
}
