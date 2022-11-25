#include <iostream>
#include <cmath>
double FuncToMin(double x); //function prototype

int main(int argc, char* argv[])
{
 double x,fx;
 std::cout << "Enter a double x:" << std::endl; //ask user for input
 std::cin >> x; // read input in
 fx=FuncToMin(x); //call function to compute f(x)
 std::cout << "Value of f(x) at x=" << x << " is " << fx << std::endl;
 return 0;
}

double FuncToMin(double x)
{
    double fx; //declare variable to be outputted
    fx=x*(x-1.0)*exp(x); //use function given in question
    return fx;
}
