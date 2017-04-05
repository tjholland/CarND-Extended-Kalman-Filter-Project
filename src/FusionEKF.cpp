#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  /**
  TODO:
    * Finish initializing the FusionEKF.
  */

  // measurement covariance
  R_laser_ << .0225, 0,
			  0, .0225;

  R_radar_ << .09, 0, 0,
			  0, .0009, 0,
			  0, 0, .09;

  // measurement matrix
  H_laser_ << 1, 0, 0, 0,
			  0, 1, 0, 0;

  Hj_ << 0, 0, 0, 0,
		 0, 0, 0, 0,
		 0, 0, 0, 0;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
	
    // lidar variables
	float px, py, vx, vy;
	// radar variables
	float rho, phi, rhodot;

	if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

	// initialize state covariance matrix P
	ekf_.P_ = MatrixXd(4, 4);
	ekf_.P_ << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1000, 0,
		0, 0, 0, 1000;

	// initialize transition matrix F_
	ekf_.F_ = MatrixXd(4, 4);

	// initialize Q 
	ekf_.Q_ = MatrixXd(4, 4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
		// calculate x_ from raw measurements
		rho = measurement_pack.raw_measurements_(0);
		phi = measurement_pack.raw_measurements_(1);
		rhodot = measurement_pack.raw_measurements_(2);

		px = rho * cos(phi);
		py = rho * sin(phi);
		vx = rhodot * cos(phi);
		vy = rhodot * sin(phi);

		ekf_.x_ << px, py, vx, vy;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
		// grab x_ from raw measurements
		ekf_.x_ << measurement_pack.raw_measurements_(0), measurement_pack.raw_measurements_(1), 0, 0;
    }

	// timestamp
	previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
   */

	// temporary variables
	float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
	float dt2, dt3, dt4;
	dt2 = pow(dt, 2);
	dt3 = pow(dt, 3);
	dt4 = pow(dt, 4);

	previous_timestamp_ = measurement_pack.timestamp_;

	// update F_
	ekf_.F_ << 1, 0, dt, 0,
		0, 1, 0, dt,
		0, 0, 1, 0,
		0, 0, 0, 1;

	// update Q_
	float noise_ax = 9;
	float noise_ay = 9;
	
	ekf_.Q_ << dt4 / 4 * noise_ax, 0, dt3 / 2 * noise_ax, 0,
		0, dt4 / 4 * noise_ay, 0, dt3 / 2 * noise_ay,
		dt3 / 2 * noise_ax, 0, dt2 * noise_ax, 0,
		0, dt3 / 2 * noise_ay, 0, dt2 * noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
	  rho = measurement_pack.raw_measurements_(0);
	  phi = measurement_pack.raw_measurements_(1);
	  rhodot = measurement_pack.raw_measurements_(2);

	  VectorXd z_radar(3);
	  z_radar << rho, phi, rhodot;

	  // get Jacobian
	  Tools tools;
	  ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
	  ekf_.R_ = R_radar_;

	  // update
	  ekf_.UpdateEKF(z_radar);
  } else {
    // Laser updates
	  VectorXd z_laser(2);
	  z_laser << measurement_pack.raw_measurements_(0), measurement_pack.raw_measurements_(1);

	  ekf_.H_ = H_laser_;
	  ekf_.R_ = R_laser_;

	  // update
	  ekf_.Update(z_laser);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
