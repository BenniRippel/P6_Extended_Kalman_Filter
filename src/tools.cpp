#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using namespace std;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  Calculation of the RMSE.
  */
    VectorXd rmse(4);
    rmse << 0,0,0,0;
    // check the validity of the inputs

    if ((estimations.size()==0) || (estimations.size()!=ground_truth.size())){
        cout<<"Error in CalculateRMSE(). Size of inputs are suspicious"<<endl;
    }


    //accumulate squared residuals
    for(int i=0; i < estimations.size(); ++i){
        // ... your code here
        VectorXd res = estimations[i] - ground_truth[i];
        res = res.array()*res.array();
        rmse += res;
    }

    //calculate the mean
    rmse = rmse.array() * 1.0/estimations.size();
    //calculate the squared root
    rmse=rmse.array().sqrt();
    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
    Calculation of Jacobian Matrix
  */

    MatrixXd Hj(3,4);
    //recover state parameters
    double px = x_state(0);
    double py = x_state(1);
    double vx = x_state(2);
    double vy = x_state(3);
    // init denominators
    double den = px*px+py*py;
    double den_sqrt = pow(den, 1.0/2.0);
    double den_3_2 = pow(den, 3.0/2.0);

    //check division by zero
    if (fabs(den)< 0.0000001){
        cout<<"Division by 0 when calculating the Jacobian Matrix. Returning all zeros..."<<endl;
        return Hj;
    }

    //compute the Jacobian matrix

    Hj(0,0) = px / den_sqrt;
    Hj(0,1) = py / den_sqrt;
    Hj(0,2) = 0;
    Hj(0,3) = 0;

    Hj(1,0) = -py / den;
    Hj(1,1) = px/den;
    Hj(1,2) = 0;
    Hj(1,3) = 0;

    Hj(2,0) = py*(vx*py-vy*px)/den_3_2;
    Hj(2,1) = px*(vy*px-vx*py)/den_3_2;
    Hj(2,2) = px / den_sqrt;
    Hj(2,3) = py / den_sqrt;

    return Hj;
}
