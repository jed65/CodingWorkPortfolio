#include <iostream>
#include <cmath>
#include <cassert>
double FuncToMin(double x); //function prototype

int main(int argc, char* argv[])
{
 double a,b;
 std::cout << "Enter endpoints of interval [a,b]:" << std::endl;
 std::cin >> a >> b; // read in the inputs
 int n;
 std::cout << "Enter positive integer number of points n:" << std::endl;
 std::cin >> n; //read in number of points
 assert(n>0); //check that positive number of points inputted, error if not
 assert(b>a); //check that the endpoint b is greater than a, error if not
 double CurrentMin, fx, xi,xmin;
 for (int i=0;i<n;i++)
 {
     xi=a+(double)i*((b-a)/((double)n-1.0)); //Calculate discretisation point
     fx=FuncToMin(xi); //Compute f(xi)
     //CHECK
     //std::cout << xi << " " << fx << std::endl;
     if (i==0)
     {
         CurrentMin=fx; //Initialise the current minimum
         xmin=xi;
     }
     else
     {
         if (CurrentMin-fx>=0.0)
         {
             CurrentMin=fx; //if current minimum greater/equal than fx, update
             xmin=xi; //also update the x value
         } //if its not bigger, it'll stay on the same value
     }
 }
 std::cout << "The minimum value is " << CurrentMin << " at x=" << xmin <<
 std::endl;
 return 0;
}

double FuncToMin(double x)
{
    double fx; //declare variable to be outputted
    fx=x*(x-1.0)*exp(x); //use function given in question
    return fx;
}
