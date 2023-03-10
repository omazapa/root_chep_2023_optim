#include "Math/ScipyMinimizer.h"
#include "Math/MultiNumGradFunction.h"
#include "Math/Functor.h"
#include <string>
#include "Math/MinimizerOptions.h"
#include "TStopwatch.h"
#include "functions.hxx"


// methods that requires hessian to work "dogleg", "trust-ncg","trust-exact","trust-krylov"
using namespace std;
int scipy()
{ 
   
   std::string methods[]={"Nelder-Mead","L-BFGS-B","Powell","CG","BFGS","TNC","COBYLA","SLSQP","trust-constr","Newton-CG", "dogleg", "trust-ncg","trust-exact","trust-krylov"};
   //Wood4GradientFunction wgf;
   TStopwatch t;
   for(const std::string &text : methods)
   {
   ROOT::Math::Experimental::ScipyMinimizer minimizer(text.c_str());
   minimizer.SetMaxFunctionCalls(1000000);
   minimizer.SetMaxIterations(100000);
   minimizer.SetTolerance(1e-3);
   minimizer.SetExtraOption("gtol",1e-3);
   ROOT::Math::GradFunctor f(&RosenBrock,&RosenBrockGrad,2); 
   double step[2] = {0.01,0.01};
   double variable[2] = { -1.0,1.0};
 
   minimizer.SetFunction(f);
   minimizer.SetHessianFunction(RosenBrockHessian);
   //minimizer.SetFunction(rgf);
   //minimizer.SetFunction(wgf);
 
   // Set the free variables to be minimized!
   minimizer.SetVariable(0,"x",variable[0], step[0]);
   minimizer.SetVariable(1,"y",variable[1], step[1]);
   t.Reset();
   t.Start();
   minimizer.Minimize(); 
   t.Stop();
   const double *xs = minimizer.X();
   cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): " 
        << RosenBrock(xs) << endl;
   cout << "Cpu Time (sec) = " << t.CpuTime() <<endl<< "Real Time (sec) = " << t.RealTime() << endl;
   cout << endl << "===============" << endl;
   }
   return 0;
}

int main()
{
  return scipy();
}