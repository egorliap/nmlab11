#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "header.h"

int main() {
	int N = 0;

	double a = 0, b = M_PI/5;
	FILE* file;
	int counter = 0;
	char path[100];
	double S, length = 0;
	int n = 3;

		for (int power = 0;power <= 10;power++) {
			for (int richardson = 0;richardson <= 1;richardson++) {
				counter = 0;
				S = find_approx_integral(f, a, b, n, power, &counter, &length,richardson);
				sprintf(path, "res/%d_%d.txt", richardson, power);
				file = fopen(path, "w+");
				fprintf(file, "%.15lf %d %.16lf", S, counter, length);
				fclose(file);

			}
		}
	

	
	return 0;
}