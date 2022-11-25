#include <iostream>
#include <cmath>
double FuncToMin(double x); //function prototypes
int PerformGoldenSearch(double(*FuncToMin)(double x),double a, double b, double
    tol, double& minimiser, int& numberIts);

int main(int argc, char* argv[])
{
    double a=0.0,b=1.0,tol=1e-10; //declare as given in question
    int IterationNumber,minFound;
    double xmin;
    minFound=PerformGoldenSearch(FuncToMin,a,b,tol,xmin,IterationNumber);
    if (minFound==0) //minimum found by function
    {
    std::cout << "Minimum at x=" << xmin << " with f(x)=" << FuncToMin(xmin) <<
    " after " << IterationNumber << " iterations." << std::endl;
    }
    else
    {
       std::cout << "Warning: Max. number of iterations reached, best guess is"
       "x=" << xmin << " with f(x)=" << FuncToMin(xmin) << std::endl;
    }
    return 0;
}

double FuncToMin(double x)
{
    double fx; //declare variable to be outputted
    fx=x*(x-1.0)*exp(x); //use function given in question
    return fx;
}

int PerformGoldenSearch(double(*FuncToMin)(double x),double a, double b, double
    tol, double& minimiser, int& numberIts)
{
  const double phi=(1.0+sqrt(5.0))/2.0; // this is golden ratio
  double x1=a,x2=b,x3,x4;
  numberIts=0;
  for (int i=1;i<=100;i++) // this goes up to max. number of iterations
  {
  if (fabs(x2-x1)<tol)
  {
      minimiser=(x1+x2)/2.0;
      break; //leave for loop
  }
  else
  {
     x3=(x2+phi*x1)/(1.0+phi);
     x4=(x1+phi*x2)/(1.0+phi);

     if (FuncToMin(x3)>std::max(FuncToMin(x1),FuncToMin(x2)) && FuncToMin(x4)>
         std::max(FuncToMin(x1),FuncToMin(x2)))
     {
         if (FuncToMin(x1)<FuncToMin(x2))
         {
             minimiser=x1; //following algorithm, minimiser on the boundary
             break; //exit the for loop
         }
         else
         {
             minimiser=x2;
             break;
         }
     }
     else
     {
         if (FuncToMin(x3)<=FuncToMin(x4))
         {
             x2=x4; //following algorithm, update values if condition holds
         }
         else
         {
             x1=x3;
         }
     }
     numberIts++; // add 1 to the number of iterations
  }
  }
  if (numberIts==100) //maximum number of iterations reached
  {
      return 1; //minimum not found
  }
  else
  {
      return 0; //minimum found
  }
}
