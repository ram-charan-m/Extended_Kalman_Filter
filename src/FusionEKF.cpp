// This class inititalizes the kalman filter when a first measurement is made.
// It calls the predict and update methods of kalman_filter class instance ekf_.
#define _USE_MATH_DEFINES
#include <math.h>
#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;
using std::pow;

// CONSTRUCTOR
FusionEKF::FusionEKF() {
  is_initialized_ = false;
  previous_timestamp_ = 0;

  // initializing EKF matrices

  // Update Matrices
  // Laser - measures in cartesian and only px and py
  MatrixXd H_laser_(2, 4); // Measurement transformation matrix of Lidar
  MatrixXd R_laser_(2, 2); // Measurement Covariance matrix of Lidar
  
  // Radar - measures in radial and range rho, bearing phi, radial velocity rho dot  
  MatrixXd R_radar_(3, 3); // Measurement covariance matrix of Radar  
  
  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;
   
  H_laser_ << 1, 0, 0, 0, // Setting Laser Measurement transformation matrix
              0, 1, 0, 0;    
}

// DESTRUCTOR
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  // Process measurement function uses KalmanFilter class instance ekf_
  // Initialization 
  if (!is_initialized_) { // For first measurement this holds true and the following loop is executed

    //Only initializing the state ekf_.x_ with the first measurement without prediction and update
    //Creating the covariance matrix.
    
    // First measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;
    
    ekf_.P_ << 1, 0, 0, 0, // setting covariance matrix
               0, 1, 0, 0,
               0, 0, 1000, 0,
               0, 0, 0, 1000;    

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Cnverting radar from polar to cartesian coordinates 
      // Normalizing phi
      cout<< " First Measurement - RADAR" << endl;
      float rho = measurement_pack.raw_measurements_(0);
      float phi = measurement_pack.raw_measurements_(1);
      float v_rad = measurement_pack.raw_measurements_(2);
      while (phi<-M_PI){
        phi = phi + 2*M_PI;
      }
      while (phi>M_PI){
        phi = phi - 2*M_PI;
      }
      float tanphi_term = pow((1+pow((tan(phi)),2)),0.5);
      float py = rho/tanphi_term;
      float px = py*tan(phi);
      float vx = v_rad*cos(phi);
      float vy = v_rad*sin(phi);
      ekf_.x_ << px, py, vx, vy;
      
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      cout << " First Measurement - LIDAR " << endl;
      ekf_.x_ << measurement_pack.raw_measurements_(0), measurement_pack.raw_measurements_(1), 0, 0;
    }
    
    previous_timestamp_ = measurement_pack.timestamp_; // Updating previous timestamp
    is_initialized_ = true; // Updating flag
    cout<< " First Measurement Complete." << endl;
    return;
  }

  //Prediction
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_; // Updating previous_timestamp
  
  ekf_.F_ << 1,0,dt,0, // setting state transition matrix
             0,1,0,dt,
             0,0,1,0,
             0,0,0,1;
  
  // Setting Process Covariance Matrix
  float noise_ax = 9;
  float noise_ay = 9;
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;
  ekf_.Q_ << 0.25*(dt_4)*noise_ax, 0, 0.5*dt_3*noise_ax, 0,
             0, 0.25*(dt_4)*noise_ay, 0, 0.5*dt_3*noise_ay,
             0.5*dt_3*noise_ax, 0, dt_2*noise_ax, 0,
             0, 0.5*dt_3*noise_ay, 0, dt_2*noise_ay;
    
  ekf_.Predict();

  // Measurement Update
  // 1. Using the sensor type to perform the update step. 2.Updating the state and covariance matrices.  

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates    
    ekf_.R_ = R_radar_;
    ekf_.H_ = tools_.CalculateJacobian(ekf_.x_);  
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);  
  } 
  else {
    // TODO: Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "New Measurement" << endl;
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
