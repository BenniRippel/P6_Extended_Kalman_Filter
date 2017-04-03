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

    // TODO: YOUR CODE HERE

    // check the validity of the following inputs:
    //  * the estimation vector size should not be zero
    //  * the estimation vector size should equal ground truth vector size
    // ... your code here
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
    // ... your code here
    rmse = rmse.array() * 1.0/estimations.size();
    //calculate the squared root
    rmse=rmse.array().sqrt();
    // ... your code here

    //return the result
    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
    Calculation of Jacobian Matrix
  */

    MatrixXd Hj(3,4);
    //recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);
    // init denominators
    float den = px*px+py*py;
    float den_sqrt = pow(den, 1.0/2.0);
    float den_3_2 = pow(den, 3.0/2.0);

    //check division by zero
    if (fabs(den)< 0.0000001){
        cout<<"Division by 0 when calculating the Jacobian Matrix. Returning all zeros..."<<endl;
        return Hj;
    }

    //compute the Jacobian matrix

    Hj(0,0) = px / den_sqrt;
    Hj(0,1) = py / den_sqrt;
    Hj(1,0) = -py / den;
    Hj(1,1) = px/den;
    Hj(2,0) = py*(vx*py-vy*px)/den_3_2;
    Hj(2,1) = px*(vy*px-vx*py)/den_3_2;
    Hj(2,2) = px / den_sqrt;
    Hj(2,3) = py / den_sqrt;

    return Hj;
}
