#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "header.h"



double polyval(double* coeffs, double x, int n) {
	double ans = 0;
	for (int i = 0;i < n + 1;i++) {
		ans += coeffs[i] * pow(x, n - i);
	}
	return ans;
}
double* polyder(double* coeffs, int n) {
	double* der = malloc(sizeof(double) * (n));
	for (int i = 0;i < n ;i++) {
		der[i] = coeffs[i] * (n - i);
	}
	return der;
}
double newton_solve(double* coeffs, int n, double a, double b) {
	double x_next=b,x_prev = b;
	double* der = polyder(coeffs, n);
	do{	
		x_prev = x_next;
		x_next = x_prev - polyval(coeffs, x_prev, n) / polyval(der, x_prev, n - 1);
	} while (fabs(x_next - x_prev) >= 10e-10);
	free(der);
	return x_next;
}

double* solve_poly(double* coeffs, int n,double a,double b) {
	double* xx = calloc(n,sizeof(double));
	int c = 0;
	double x1 = a, x2 = a;
	while (x1 < b) {
		while (x2 < x1) {
			if (polyval(coeffs, x1, n) * polyval(coeffs, x2, n) < 0) {
				xx[c] = newton_solve(coeffs, n, x2, x1);
				c++;
				x1 = x2;
				x2 += (b - a) / 10000;
				break;
			}
			x2 += (b - a) / 10000;
		}
		x1 += (b - a) / 10000;
	}
	return xx;
}
double** table_xx() {
	double** ret = malloc(sizeof(double*) * 4);
	
	for (int n = 2;n <= 5;n++) {
		double* coeffs = malloc(sizeof(double) * (n + 1));
		double* sk = malloc(sizeof(double) * (n + 1));
		if (coeffs == NULL || sk == NULL) {
			exit(16);
		}
		coeffs[0] = 1;
		
		double A = 2.0 / n;
		for (int i = 1;i < (n + 1);i++) {
			sk[i] = (1.0 / (i + 1) - pow(-1, i + 1) / (i + 1)) / A;
			coeffs[i] = sk[i];
			for (int j = 1;j < i;j++) {
				coeffs[i] += coeffs[j] * sk[i - j];
			}
			coeffs[i] *= -1.0 / i;
		}

		ret[n-2] = solve_poly(coeffs, n, -1, 1);

		free(sk);
		free(coeffs);
	}
	return ret;
}
double chebyshev_method(double(*f)(double), double a1, double b1, int n, double** xx_tab) {
	double s = 0, a, b, ans;

	a = a1;
	b = b1;
	double* xx = xx_tab[n - 2];
	double A = 2.0 / n;

	ans = 0;
	for (int k = 0;k < n;k++) {
		ans += f((a + b) / 2 + ((b - a) / 2) * xx[k]);
	}
	s = ((b - a) / 2) * ans * A;



	return s;
}


double adaptive_div(double(*f)(double),double** xx_tab,double a,double b,int n,int power, double S_n,double S_n_plus_1,int *N,double* len) {

	if (fabs(S_n_plus_1 - S_n) < pow(10,-power)/(*N)) {
		return S_n;
	}
	else {
		double S_n_plus_1_r = chebyshev_method(f, a, (a + b) / 2, n+1, xx_tab);
		double S_n_r = chebyshev_method(f, a, (a + b) / 2, n, xx_tab);
		double S_n_plus_1_l = chebyshev_method(f, (a + b) / 2, b, n+1, xx_tab);
		double S_n_l = chebyshev_method(f, (a + b) / 2, b,n, xx_tab);
		*N += 1;
		
		if (*len > (b - a) / (*N)) {
			*len = (b-a)/(*N);
		}

		return (adaptive_div(f, xx_tab, a, (a + b) / 2,  n, power, S_n_l, S_n_plus_1_l,N,len) + adaptive_div(f, xx_tab, (a + b) / 2, b,  n, power, S_n_r, S_n_plus_1_r,N,len));
	}
	
}

double adaptive_div_Richardson(double** xx_tab, double a, double b, double(*f)(double), int n, int power, double S_prev, double S_new, int* N, double* len) {

	if (fabs(S_new - S_prev)/(pow(2,n)-1) < pow(10, -power) / ( (*N))) {
		return S_new;
	}
	else {
		double S_newl = chebyshev_method(f, a, (a + b) / 2, n + 1, xx_tab);
		double S_prevl = chebyshev_method(f, a, (a + b) / 2, n, xx_tab);
		double S_newr = chebyshev_method(f, (a + b) / 2, b, n + 1, xx_tab);
		double S_prevr = chebyshev_method(f, (a + b) / 2, b, n, xx_tab);
		*N += 1;

		if (*len > (b - a) / (*N)) {
			*len = (b - a) / (*N);
		}



		return (adaptive_div_Richardson(xx_tab, a, (a + b) / 2, f, n, power, S_prevl, S_newl, N, len) + adaptive_div_Richardson(xx_tab, (a + b) / 2, b, f, n, power, S_prevr, S_newr, N, len));
	}

}

double find_approx_integral(double(*f)(double), double a, double b, int n, int power, int* counter,double* length,int richardson) {
	int N = 1;
	double** xx_table = table_xx(a, b);
	
	double S2;
	double len = 10000;
	if (richardson) {
		S2 = adaptive_div_Richardson(xx_table, a, b, f, n, power, 1000000, 100, &N, &len);
	}
	else {
		S2 = adaptive_div(f,xx_table, a, b,  n, power, 1000000, 100, &N, &len);
	}
	printf("N = %ld, n = %d, S = %.15lf, eps = %e, len=%.16lf\n", N,n, S2, pow(10, -power),len);
	*counter = N;
	*length = len;
	return S2;
}


