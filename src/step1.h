/*
 Jiaxing Lin
 Aug 25 2016
 */


/*
 function to compute the feasible seach direction for the 
 general gradient projection method.
 */


VectorXd feaDirec(MatrixXd &hessian, VectorXd &gradient, MatrixXd &X)
{
  // compute the inversion of hessian matrix, which seem unevidable
  double delta = 0.05;
  MatrixXd eye = MatrixXd::Identity(hessian.rows(),hessian.cols());
  MatrixXd W = -hessian + delta*eye;
  MatrixXd invHessian = W.inverse();


  VectorXd fdirection = (eye - invHessian*X.transpose()*(X*invHessian*X.transpose()).inverse()*X)*invHessian*gradient;
  
  // set small value to zero.
  for(int i = 0; i < fdirection.size(); i++)
    if(abs(fdirection(i)) < 1e-7/5)
      fdirection(i) = 0;

  return fdirection;
}










