#define _USE_MATH_DEFINES
#include <math.h>
int sgn(double val) {
	return (0 < val) - (val < 0);
}
double f(double x) {
	return (pow(x, 5) - 2.9 * pow(x, 3) + 6.5 * pow(x, 2) - 7 * x)/cos(2*x);
}
double f_modified(double x) {
	return sgn(pow(x, 5) - 2.9 * pow(x, 3) + 6.5 * pow(x, 2) - 7 * x  + M_PI/2) / cos(2 * x);
}
