#ifndef UKF_H
#define UKF_H
#include "Eigen/Dense"
#include "measurement_package.h"
#include "ground_truth_package.h"
#include <vector>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* State dimension
  int n_x_;

  ///* set measurement dimension, lidar can measure p_x and p_y
  int n_y_;

  ///* Augmented state dimension
  int n_aug_;

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z_;

  ///* Sigma point spreading parameter
  int lambda_;

  ///*previous timestamp
  long previous_timestamp_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* augmented state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate nu_a nu_yaw_dd]
  VectorXd x_aug_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* augmented state covariance matrix
  MatrixXd P_aug_;

  ///* process noise covariance matrix
  MatrixXd Q_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* calculate square root of P
  MatrixXd A_;

  ///* create example sigma point matrix
  MatrixXd Xsig_aug_;

  ///* mean predicted measurement
  VectorXd y_pred_;

  ///* create matrix for sigma points in measurement space
  MatrixXd Ysig_;

  ///* mean predicted measurement
  VectorXd z_pred_;

  ///* create matrix for sigma points in measurement space
  MatrixXd Zsig_;

  ///* measurement covariance matrix S
  MatrixXd S_z_;

  ///* measurement covariance matrix S
  MatrixXd S_y_;

  ///* add measurement noise covariance matrices
  MatrixXd R_laser_;
  MatrixXd R_radar_;

  ///* delta t between measurements
  float dt_;

  ///* minimum value
  const double eps_ = 0.00001;

  ///* the current NIS for radar
  double NIS_radar_;

  ///* the current NIS for laser
  double NIS_laser_;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   * @param gt_package The ground truth of the state x at measurement time
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   */
  void UpdateRadar(MeasurementPackage meas_package);

  /**
  * Generate augemented sigma points
  */
  void GenerateAugmentedSigmaPoints(const VectorXd &x, const MatrixXd &P, MatrixXd &Xsig_aug);

  /**
  * Predicts sigma points in new step
  */
  void SigmaPointPrediction(const MatrixXd &Xsig_aug, MatrixXd &Xsig_pred, double delta_t);

  /**
  * Predicts new mean and covariance for distribution of predicted sigma points
  */
  void PredictMeanAndCovariance(const MatrixXd &Xsig_pred, VectorXd &x, MatrixXd &P);

  /**
  * Applies update step using either a radar or laser measurement
  */
  void PredictUpdateMeasurement(const MeasurementPackage meas_package, VectorXd &z_pred, MatrixXd &S, const MatrixXd &Xsig_pred, MatrixXd &Zsig, VectorXd &x, MatrixXd &P);
};

#endif /* UKF_H */
