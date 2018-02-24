#include <iostream>
#include "ukf.h"

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  //check if UKF is successfully initialised
  is_initialized_ = false;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.6; //9 or 0.9 or 0.09 or 3 or 0.2

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.2; //3 or 0.3 or 0.03 or 0.5 or 0.2

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.1; //0.07

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.1; //0.07

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  // timestamp of last measurement for delta_t
  previous_timestamp_ = 0;

  //set state dimension
  n_x_ = 5;

  //set augmented dimension
  n_aug_ = 7;

  //define spreading parameter
  lambda_ = 3 - n_aug_;

  //set measurement dimension, radar can measure p_x and p_y
  n_y_ = 2;

  //set measurement dimension, radar can measure r, phi, and r_dot
  n_z_ = 3;

  //set example state
  x_ = VectorXd(n_x_);

  //set example covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  //Initial values
  P_ <<   0.0043, -0.0013, 0.0030, -0.0022, -0.0020,
          -0.0013, 0.0077, 0.0011, 0.0071, 0.0060,
          0.0030, 0.0011, 0.0054, 0.0007, 0.0008,
          -0.0022, 0.0071, 0.0007, 0.0098, 0.0100,
          -0.0020, 0.0060, 0.008, 0.0100, 0.0123;

  //create matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  //set vector for weights_
  weights_ = VectorXd((2*n_aug_+1));

  //initialise weights
  //weights_(0)= lambda_/(lambda_ + n_aug_);
  weights_(0)= (-4.)/(3.);

  for (int i=1; i<(2*n_aug_+1); i++) {
    weights_(i) = 0.5/(n_aug_+lambda_);
  }

  //set process noise covariance matrix
  Q_ = MatrixXd(2, 2);
  Q_ <<   std_a_ * std_a_, 0,
          0, std_yawdd_ * std_yawdd_;

  //add measurement noise covariance matrix
  R_radar_ = MatrixXd(n_z_,n_z_);

  R_radar_ << std_radr_ * std_radr_, 0, 0,
          0, std_radphi_ * std_radphi_, 0,
          0, 0,std_radrd_ * std_radrd_;

  R_laser_ = MatrixXd(2, 2);

  R_laser_ << std_laspx_ * std_laspx_, 0,
          0, std_laspy_ * std_laspy_;

  //create sigma point matrix
  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //set augmented state
  x_aug_ = VectorXd(n_aug_);

  //set augmented covariance matrix
  P_aug_ = MatrixXd(n_aug_, n_aug_);

  //calculate square root of P_aug
  A_ = P_aug_.llt().matrixL();

  //mean predicted measurement
  z_pred_ = VectorXd(n_z_);

  //measurement covariance matrix S
  S_z_ = MatrixXd(n_z_, n_z_);

  //create matrix for sigma points in measurement space
  Zsig_ = MatrixXd(n_z_, 2 * n_aug_ + 1);

  //mean predicted measurement
  y_pred_ = VectorXd(n_y_);

  //measurement covariance matrix S
  S_y_ = MatrixXd(n_y_, n_y_);

  //create matrix for sigma points in measurement space
  Ysig_ = MatrixXd(n_y_, 2 * n_aug_ + 1);

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  /*****************************************************************************
 *  Initialization
 ****************************************************************************/
  if (!is_initialized_) {

    //set timestamp
    previous_timestamp_ = meas_package.timestamp_;

    /**
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * convert radar from polar to cartesian coordinates.
    */

    // first measurements
    std::cout << "UKF Initialisation: " << std::endl;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      std::cout << "Unscented Kalman Filter Initialization with RADAR" << std::endl;

      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */

      double rho_ = meas_package.raw_measurements_[0];
      double phi_ = meas_package.raw_measurements_[1];
      double rho_dot_ = meas_package.raw_measurements_[2];

      /**
       * Initialize state.
       */

      //set the state with the initial location and zero velocity
      x_ << rho_ * cos(phi_), rho_ * sin(phi_), 0, 0, 0;

    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      std::cout << "Unscented Kalman Filter Initialization with LIDAR" << std::endl;

      /**
      Initialize state.
      */

      //set the state with the initial location and zero velocity directly
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }

    // check if initialisation was successful
    if ((x_[0] != 0.) and (x_[1] != 0.)) {
      is_initialized_ = true;
    }

    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  //get time passed since last measurement
  dt_ = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;

  //https://discussions.udacity.com/t/numerical-instability-of-the-implementation/230449/2
  //do prediction step
  //Prediction(dt_);

  while (dt_ > 0.2)
  {
      std::cout << "Large dt_ value " << dt_ << std::endl;
      const double dt_step = 0.1;
      Prediction(dt_step);
      dt_ -= dt_step;
  }

  Prediction(dt_);

  //set previous timestamp
  previous_timestamp_ = meas_package.timestamp_;


  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {

    /**
     * Radar update
     * measurement update
     */

    UpdateRadar(meas_package);

  } else {
    /**
     * Laser update
     * measurement update
     */

    UpdateLidar(meas_package);

  }

}
/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  /*****************************************************************************
  * Prediction
  * Generate sigma points X a,k|k => Augmented State vector(px, py, v, yaw, yaw_dot, nu_a, nu_yaw_dot_dot) -> 7x15 dimension
  * Predict sigma points X k+1|k => Predicted state Vector (px,py, v, yaw, yaw_dot) -> 5x15 dimension
  * Predict mean and covariance
  ****************************************************************************/

  //create sigma points
  GenerateAugmentedSigmaPoints(x_, P_, Xsig_aug_);

  // predict sigma points
  SigmaPointPrediction(Xsig_aug_, Xsig_pred_, delta_t);

  // Predict new mean and covariance
  PredictMeanAndCovariance(Xsig_pred_, x_, P_);
}

/**
 * Updates the state and the state covariance matrix using A_laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

/**
  Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  Calculate the lidar NIS.
  */


  /*****************************************************************************
 *  Update
 ****************************************************************************/

  /**
     * Use Lidar to perform the update step.
     * Update the state and covariance matrices.
   */

  if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    PredictUpdateMeasurement(meas_package, y_pred_, S_y_, Xsig_pred_, Ysig_, x_, P_);
  }

  // print the output
  std::cout << "x_ = " << x_ << std::endl;
  std::cout << "P_ = " << P_ << std::endl;

  /**
   * expected result:
   * x_ = 5
   * P_ = 5x5
   */
}

/**
 * Updates the state and the state covariance matrix using A_radar measurement.
 * @param {MeasurementPackage} meas_package
 */

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**

  Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.
  Calculate the radar NIS.

  */

  /*****************************************************************************
 *  Update
 ****************************************************************************/

  /**
     * Use the radar to perform the update step.
     * Update the state and covariance matrices.
   */

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {

    PredictUpdateMeasurement(meas_package, z_pred_, S_z_, Xsig_pred_, Zsig_, x_, P_);
  }

  // print the output
  std::cout << "x_ = " << x_ << std::endl;
  std::cout << "P_ = " << P_ << std::endl;

  /**
   * expected result:
   * x_ = 5
   * P_ = 5x5
   */

}

void UKF::GenerateAugmentedSigmaPoints(const VectorXd &x, const MatrixXd &P, MatrixXd &Xsig_aug) {

/*******************************************************************************
 * adapted Udacity Sample Code
 ******************************************************************************/

  // create augmented state
  x_aug_ << x, 0, 0;

  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(n_x_, n_x_) = P;
  P_aug_.bottomRightCorner(Q_.rows(), Q_.cols()) = Q_;

  // calculate square root of P
  A_ = P_aug_.llt().matrixL();

  //create sigma points
  Xsig_aug.col(0) = x_aug_;

  for (int i = 0; i < n_aug_; i++)
  {
    Xsig_aug.col(i+1)     = x_aug_ + sqrt(lambda_+n_aug_) * A_.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_+n_aug_) * A_.col(i);
  }

  //print result
  std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;

  /**
   * expected result:
   *  Xsig_aug = 7 x 15
   */

}

void UKF::SigmaPointPrediction(const MatrixXd &Xsig_aug, MatrixXd &Xsig_pred, double delta_t) {

  /*******************************************************************************
 * adapted Udacity sample code
 ******************************************************************************/

  //predict sigma points
  for (int i = 0; i < (2*n_aug_+1); i++) {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
      //std::cout << "no division by zero" << std::endl;
      px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
      py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
      //std::cout << "division by zero" << std::endl;
      px_p = p_x + v*delta_t*cos(yaw);
      py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred(0,i) = px_p;
    Xsig_pred(1,i) = py_p;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;
  }

  std::cout << "Xsig_pred_ = " << std::endl << Xsig_pred_ << std::endl;

  /**
   * expected result:
   * Xsig_pred = 5 x 15
   */
}

void UKF::PredictMeanAndCovariance(const MatrixXd &Xsig_pred, VectorXd &x, MatrixXd &P) {

  /*******************************************************************************
 * adapted Udacity sample code
 ******************************************************************************/

  //predicted state mean
  x.fill(0.0);
  for (int i = 0; i < (2 * n_aug_ + 1); i++) {  //iterate over sigma points
    x = x + weights_(i) * Xsig_pred.col(i);
  }

  //angle normalization
  x(3) = x(3) - ceil((x(3) - M_PI) / (2. * M_PI)) * 2. * M_PI;
  x(4) = x(4) - ceil((x(4) - M_PI) / (2. * M_PI)) * 2. * M_PI;

  //predicted state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < (2 * n_aug_ + 1); i++) {  //iterate over sigma points

    // state difference
    //VectorXd x_diff = Xsig_pred.col(i) - x;

    //https://discussions.udacity.com/t/numerical-instability-of-the-implementation/230449/14
    VectorXd x_diff = Xsig_pred.col(i) - Xsig_pred.col(0);

    //angle normalization
    x_diff(3) = x_diff(3) - ceil((x_diff(3) - M_PI) / (2. * M_PI)) * 2. * M_PI;
    x_diff(4) = x_diff(4) - ceil((x_diff(4) - M_PI) / (2. * M_PI)) * 2. * M_PI;

    P = P + weights_(i) * x_diff * x_diff.transpose() ;
  }

  //print result
  std::cout << "Predicted state" << std::endl;
  std::cout << x_ << std::endl;
  std::cout << "Predicted covariance matrix" << std::endl;
  std::cout << P_ << std::endl;

  /**
   * expected result:
   * x_ = 5
   * P = 5x5
   */

}


void UKF::PredictUpdateMeasurement(const MeasurementPackage meas_package, VectorXd &z_pred, MatrixXd &S, const MatrixXd &Xsig_pred, MatrixXd &Zsig, VectorXd &x, MatrixXd &P) {


/*******************************************************************************
 * adapted Udacity sample code
 * Put everything in one function for Laser AND Radar to avoid code redundancy
 ******************************************************************************/

  // Laser measurement model

  if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    //transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points

      // extract values for better readibility
      double p_x = Xsig_pred(0, i);
      double p_y = Xsig_pred(1, i);

      // measurement model
      Zsig(0, i) = p_x;
      Zsig(1, i) = p_y;
    }
  }

  // Radar measurement model

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    //transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

      // extract values for better readibility
      double p_x = Xsig_pred(0, i);
      double p_y = Xsig_pred(1, i);

      double rho = sqrt(p_x * p_x + p_y * p_y);

      // measurement model
      Zsig(0, i) = rho;       //rho

      if (fabs(rho) > 0.001) {
        double v = Xsig_pred(2, i);
        double yaw = Xsig_pred(3, i);

        double v1 = cos(yaw) * v;
        double v2 = sin(yaw) * v;

        double phi = atan2(p_y, p_x);
        double rho_dot = (p_x * v1 + p_y * v2) / std::max(eps_, rho);

        Zsig(1, i) = phi;       //phi
        Zsig(2, i) = rho_dot;   //r_dot
      } else {
        Zsig(1, i) = 0.0;       //phi
        Zsig(2, i) = 0.0;       //r_dot
      }
    }
  }

  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points
    //residual
    VectorXd y_diff = Zsig.col(i) - z_pred;

    //angle normalization
    y_diff(1) = y_diff(1) - ceil((y_diff(1) - M_PI) / (2. * M_PI)) * 2. * M_PI;

    S = S + weights_(i) * y_diff * y_diff.transpose();
  }

  if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    S = S + R_laser_;
  }

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    S = S + R_radar_;
  }

  VectorXd z = meas_package.raw_measurements_;
  int n_z = z_pred.rows();
  int n_x = x.rows();

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    z_diff(1) = z_diff(1) - ceil((z_diff(1) - M_PI) / (2. * M_PI)) * 2. * M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x;

    //angle normalization
    x_diff(3) = x_diff(3) - ceil((x_diff(3) - M_PI) / (2. * M_PI)) * 2. * M_PI;
    x_diff(1) = x_diff(1) - ceil((x_diff(1) - M_PI) / (2. * M_PI)) * 2. * M_PI;

    // calculate cross correlation Tc_
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  z_diff(1) = z_diff(1) - ceil((z_diff(1) - M_PI) / (2. * M_PI)) * 2. * M_PI;

  //update state mean and covariance matrix
  x = x + K * z_diff;

  //angle normalisation
  x(3) = x(3) - ceil((x(3) - M_PI) / (2. * M_PI)) * 2. * M_PI;
  x(4) = x(4) - ceil((x(3) - M_PI) / (2. * M_PI)) * 2. * M_PI;


  // calculate P
  P = P - K * S * K.transpose();

  if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    //Calculate NIS
    //residual
    VectorXd z_diff = z - z_pred;
    NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
    std::cout << "NIS_laser_ = " << NIS_laser_ << std::endl;
  }

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    //Calculate NIS
    //residual
    VectorXd z_diff = z - z_pred;
    NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
    std::cout << "NIS_radar_ = " << NIS_radar_ << std::endl;
  }

  //print result
  std::cout << "Updated state x: " << std::endl << x_ << std::endl;
  std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;

  /**
   * expected result P:
   * x_ = 5
   * P_ = 5x5
   */


}