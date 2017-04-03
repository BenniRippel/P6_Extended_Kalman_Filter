#include "kalman_filter.h"
#include "tools.h"
#include "FusionEKF.h"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

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

void KalmanFilter::Predict(float delta_T) {
  /**
  TODO:
    * predict the state
  */
    cout<<"R: "<<R_<<endl;

    F_ << 1, 0, delta_T, 0,
        0, 1, 0, delta_T,
        0, 0, 1, 0,
        0, 0, 0, 1;

  int noise_ax = 9;
  int noise_ay =9;

  float dt_p4 = pow(delta_T, 4);
  float dt_p3 = pow(delta_T, 3);
  float dt_p2 = pow(delta_T, 2);
  Q_ << (noise_ax*dt_p4/4), 0, (noise_ax*dt_p3/2), 0,
          0, (noise_ay*dt_p4/4), 0, (noise_ay*dt_p3/2),
          (noise_ax*dt_p3/2), 0, (noise_ax*dt_p2), 0,
          0, (noise_ay*dt_p3/2), 0, (noise_ay*dt_p2);

  x_=F_*x_;
  P_=F_*P_*F_.transpose()+Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */

  VectorXd z_pred;    //h(x) - vorherige pos in pol koord (3 eintraege)
  float c1 = z[0]*z[0]+z[1]*z[1];
  float c2 = pow(c1, 0.5);
  z_pred<<c2, atan2(z[1],z[0]), (z[0]*z[2]+z[1]*z[3])/c2;

  /////////////////////////////
  Tools tools;
  MatrixXd Hj = tools.CalculateJacobian(x_);

  VectorXd y = z - z_pred;
  MatrixXd Ht = Hj.transpose();
  MatrixXd S = Hj * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;


}
