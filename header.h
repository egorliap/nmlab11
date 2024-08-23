double f(double x);
double chebyshev_method(double(*f)(double), double a1, double b1, int n, double** tab);
double find_approx_integral(double(*f)(double), double a, double b, int n, int power, int* counter, double* length, int richardson);
double f_modified(double x);