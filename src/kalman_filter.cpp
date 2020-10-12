#include "kalman_filter.h"
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;
using std::pow;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

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
  // predict the state
  x_ = F_*x_;
  P_ = F_*P_*F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  // LIDAR
  // Updating the state by using Kalman Filter equations
  VectorXd y = z - H_*x_;
  MatrixXd S = H_*P_*H_.transpose() + R_;
  MatrixXd K = P_*H_.transpose()*S.inverse();
  x_ = x_ + K*y;
  MatrixXd I = MatrixXd::Identity(4, 4);
  P_ = (I - K*H_)*P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  // RADAR
  // Updating the state by using Extended Kalman Filter equations

  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);
  float phii = py/px;
  VectorXd Hf(3);

  if ((px && py)!=0){
    float H_radar0 = pow((pow(px,2) + pow(py,2)),0.5);
    while (phii < -M_PI){
      phii = phii + 2*M_PI;
    }
    while (phii > M_PI) {
      phii = phii - 2*M_PI;
    }

    float H_radar1 = atan(phii);
    float H_radar2 = (px*vx + py*vy)/H_radar0;
    Hf << H_radar0, H_radar1, H_radar2;
  }
  else{
    cout<<"ERROR in H-function calculation. NaNs detected!"<<endl;
  }    
    
  VectorXd y = z - Hf;
  MatrixXd S = H_*P_*H_.transpose() + R_;
  MatrixXd K = P_*H_.transpose()*S.inverse();
  x_ = x_ + K*y;
  MatrixXd I = MatrixXd::Identity(4, 4);
  P_ = (I - K*H_)*(P_);
}
