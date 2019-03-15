#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
    VectorXd rmse(4);
   rmse << 0, 0, 0, 0;
  
  	if (estimations.size() == 0 || estimations.size() != ground_truth.size()  ) 
	{
		cout << "something wrong with the vector size" << endl;
		return rmse;
	}

	//accumulate squared residuals
	for (unsigned int i = 0; i < estimations.size(); ++i) {

		VectorXd residual = estimations[i] - ground_truth[i];

		//coefficient-wise multiplication
		residual = residual.array()*residual.array();
		rmse += residual;
	}

	rmse = rmse / estimations.size();
	rmse = rmse.array().sqrt();
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3, 4);
  //recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  //pre-compute a set of terms to avoid repeated calculation

  float rho  = sqrt(pow(px, 2) + pow(py, 2));
  float rho2 = pow( rho, 2 );
  float rho3 = pow( rho, 3 );

  //check division by zero
  if (fabs(rho) < 0.00001){
      std::cout << "CalculateJacobian - Error - Division by Zero" << std::endl;
      return Hj;
  }

  //compute the Jacobian matrix
  Hj << px / rho, py / rho, 0, 0,
  -py / rho2, px / rho2, 0, 0,
  py*(vx*py - vy*px) / rho3, px*(vy*px - vx*py) / rho3, px / rho, py / rho;

  return Hj;
}
