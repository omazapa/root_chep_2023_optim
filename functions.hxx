
double RosenBrock(const double *xx )
{
  const Double_t x = xx[0];
  const Double_t y = xx[1];
  const Double_t tmp1 = y-x*x;
  const Double_t tmp2 = 1-x;
  return 100*tmp1*tmp1+tmp2*tmp2;
}

double RosenBrockGrad(const double *x, unsigned int ipar)
{
   if (ipar == 0)
      return -2 * (1 - x[0]) + 200 * (x[1] - x[0] * x[0]) * (-2 * x[0]);
   else
      return 200 * (x[1] - x[0] * x[0]);
}

bool RosenBrockHessian(const std::vector<double> &xx, double *hess)
{
   const double x = xx[0];
   const double y = xx[1];
   
   hess[0] = 1200*x*x - 400*y + 2;
   hess[1] = -400*x;
   hess[2] = -400*x;
   hess[3] = 200;
   
   return true;
}
