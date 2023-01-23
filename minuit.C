#include "Minuit2/Minuit2Minimizer.h"
#include "Math/MultiNumGradFunction.h"
#include "Math/Functor.h"
#include <string>
#include "Math/MinimizerOptions.h"
#include "TStopwatch.h"
#include "functions.hxx"

using namespace std;
int minuit()
{
  map<string, vector<string>> methods_map;

  methods_map["Minuit2"] = {"Migrad", "Simplex", "Combined", "Scan", "Fumili2"};
  methods_map["GSLMultiMin"] = {"ConjugateFR", "ConjugatePR", "BFGS", "BFGS2", "SteepestDescent"};
//methods_map["Fumili"] = {""}; // this is for likehood fitting
  methods_map["GSLSimAn"] = {""};
  methods_map["Genetic"] = {""};

  TStopwatch t;
  for (auto const &[minName, algoNames] : methods_map)
  {
    cout << minName << endl; // string (key)

    for (auto const &algoName : algoNames)
    {
      ROOT::Math::Minimizer *minimizer =
          ROOT::Math::Factory::CreateMinimizer(minName, algoName);
      if (!minimizer)
      {
        std::cerr << "Error: cannot create minimizer \"" << minName
                  << "\". Maybe the required library was not built?" << std::endl;
        return 1;
      }
      minimizer->SetMaxFunctionCalls(1000000);
      minimizer->SetMaxIterations(100000);
      minimizer->SetTolerance(1e-3);
      ROOT::Math::GradFunctor f(&RosenBrock, &RosenBrockGrad, 2);
      double step[2] = {0.01, 0.01};
      double variable[2] = {-1.0, 1.0};

      minimizer->SetFunction(f);
      minimizer->SetHessianFunction(RosenBrockHessian); // RMinimizer doesn't have support for heesian.

      // Set the free variables to be minimized!
      minimizer->SetVariable(0, "x", variable[0], step[0]);
      minimizer->SetVariable(1, "y", variable[1], step[1]);
      t.Reset();
      t.Start();
      minimizer->Minimize();
      t.Stop();
      const double *xs = minimizer->X();
      cout << minName << "   " << algoName << endl;
      cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): "
           << RosenBrock(xs) << endl;
      cout << "Cpu Time (sec) = " << t.CpuTime() << endl
           << "Real Time (sec) = " << t.RealTime() << endl;
      cout << endl
           << "===============" << endl;
    }
  }
  return 0;
}

int main()
{
  return minuit();
}