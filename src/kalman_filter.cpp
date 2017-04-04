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
    // predict the state

    // update F_ and Q_
    F_(0, 2) = delta_T;
    F_(1, 3) = delta_T;

    double dt_p4 = pow(delta_T, 4);
    double dt_p3 = pow(delta_T, 3);
    double dt_p2 = pow(delta_T, 2);
    double nax_p2 = pow(noise_ax, 1);
    double nay_p2 = pow(noise_ay, 1);
    Q_ << (nax_p2*dt_p4/4), 0, (nax_p2*dt_p3/2), 0,
            0, (nay_p2*dt_p4/4), 0, (nay_p2*dt_p3/2),
            (nax_p2*dt_p3/2), 0, (nax_p2*dt_p2), 0,
            0, (nay_p2*dt_p3/2), 0, (nay_p2*dt_p2);

    x_=F_*x_;
    MatrixXd Ft = F_.transpose();
    P_=F_*P_*Ft+Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
    // update the state by using Kalman Filter equation
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
  // update the state by using Extended Kalman Filter equations
    double c1 = x_[0]*x_[0]+x_[1]*x_[1];
    double c2 = pow(c1, 0.5);

    Eigen::Vector3d z_pred(c2, atan2(x_[1],x_[0]), (x_[0]*x_[2]+x_[1]*x_[3])/c2);    //h(x) - vorherige pos in pol koord (3 eintraege)

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
