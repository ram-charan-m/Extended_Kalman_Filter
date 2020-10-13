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

  // initializing matrices that vary with sensor type (Jacobian is directly assigned to ekf_.H_ below)
  H_laser_ = MatrixXd(2, 4); // Measurement transformation matrix of Lidar
  R_laser_ = MatrixXd(2, 2); // Measurement Covariance matrix of Lidar 
  R_radar_ = MatrixXd(3, 3); // Measurement covariance matrix of Radar  
  R_laser_ << 0.0225, 0, //measurement covariance matrix - laser
              0, 0.0225;
  R_radar_ << 0.09, 0, 0, //measurement covariance matrix - radar
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
    
    ekf_.P_ = MatrixXd(4,4);
    ekf_.P_ << 1, 0, 0, 0, // setting covariance matrix
               0, 1, 0, 0,
               0, 0, 1000, 0,
               0, 0, 0, 1000;    

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Cnverting radar from polar to cartesian coordinates 
      // Normalizing phi
      cout<< " First Measurement Sensor: RADAR" << endl;
      float rho = measurement_pack.raw_measurements_(0);
      float phi = measurement_pack.raw_measurements_(1);
//       cout<<"Phi = "<<phi<<endl;
      float v_rad = measurement_pack.raw_measurements_(2);
//       while (phi<-M_PI){
//         phi = phi + 2*M_PI;
//       }
//       while (phi>M_PI){
//         phi = phi - 2*M_PI;
//       }
//       cout<<"Normalized phi = "<<phi<<endl;
      
      float px = cos(phi)*rho;
      float py = sin(phi)*rho;
      float vx = v_rad*cos(phi);
      float vy = v_rad*sin(phi);
      if (px < 0.0001){ // Redundancy to avoid NaNs
        px = 0.0001;
      }
      if (py < 0.0001){
        py = 0.0001;
      }
      ekf_.x_ << px, py, vx, vy;
      
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      cout << " First Measurement Sensor: LIDAR " << endl;
      ekf_.x_ << measurement_pack.raw_measurements_(0), measurement_pack.raw_measurements_(1), 0, 0;
    }
    
    previous_timestamp_ = measurement_pack.timestamp_; // Updating previous timestamp
    is_initialized_ = true; // Updating flag
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
    cout<< " First Measurement Complete." << endl;
    return;
  }

  //Prediction
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_; // Updating previous_timestamp
  
  ekf_.F_ = MatrixXd(4,4);
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
  ekf_.Q_ = MatrixXd(4,4);
  ekf_.Q_ << 0.25*(dt_4)*noise_ax, 0, 0.5*dt_3*noise_ax, 0,
             0, 0.25*(dt_4)*noise_ay, 0, 0.5*dt_3*noise_ay,
             0.5*dt_3*noise_ax, 0, dt_2*noise_ax, 0,
             0, 0.5*dt_3*noise_ay, 0, dt_2*noise_ay;
    
  ekf_.Predict();
  cout << "New Measurement" << endl;

  // Measurement Update
  // 1. Using the sensor type to perform the update step. 2.Updating the state and covariance matrices.  

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    //Radar updates    
    cout<< "Sensor Type: RADAR" << endl;
    ekf_.R_ = R_radar_;
    ekf_.H_ = tools_.CalculateJacobian(ekf_.x_);  
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);  
  } 
  else {
    //Laser updates
    cout<< "Sensor Type: LIDAR" << endl;
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // Print the output  
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
