#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <iomanip>
using namespace std;
void GAUS(double m[2][2], double b[2], double x[2]);
double sqo(double x[7], double y[7], double a[2], int N, int m);
void main()
{
	double t[] = { 60, 70, 80, 90, 100, 110, 120 };
	double h1[] = { 0.0148, 0.0124, 0.0102, 0.0085, 0.0071, 0.0059, 0.0051 };
	double h[7];
	int N = 7, m = 1;

	for (int i = 0; i < N; i++)
		h[i] = log(h1[i]);

	double Po[3] = { 0 };
	Po[0] = N;

	for (int k = 1; k <= 2 * m; k++)
		for (int i = 0; i < N; i++)
			Po[k] += pow(t[i], k);

	double S[2][2] = { 0 };

	for (int i = 0; i < m + 1; i++)
		for (int j = 0; j < m + 1; j++)
			S[i][j] = Po[i + j];

	double pr[3] = { 0 };

	for (int l = 0; l < m + 1; l++)
		for (int i = 0; i < N; i++)
			pr[l] += h[i] * pow(t[i], l);

	double koef[2];

	GAUS(S, pr, koef);

	double a = pow(M_E, koef[0]);
	cout << " a = " << a << endl;
	double b = -koef[1];
	cout << " b = " << b << endl;

	cout << "sigma = " << sqo(t, h, koef, N, m) << endl;

}

void GAUS(double m1[2][2], double b1[2], double x1[2])
{
	int p[2] = { 0,1 };

	for (int i = 0; i < 2; i++)
	{
		int m = p[i];
		for (int j = i + 1; j < 2; j++)
			if (m1[p[m]][i] < m1[p[j]][i])
				m = j;
		if (m != p[i])
		{
			p[i] = m;
			p[m] = i;
		}

		b1[p[i]] /= m1[p[i]][i];
		for (int j = 1; j >= i; j--)
			m1[p[i]][j] /= m1[p[i]][i];

		for (int j = i + 1; j < 2; j++)
		{
			b1[p[j]] -= b1[p[i]] * m1[p[j]][i];
			for (int k = 1; k >= i; k--)
				m1[p[j]][k] -= m1[p[i]][k] * m1[p[j]][i];
		}
	}

	for (int i = 1; i >= 0; i--)
	{
		x1[i] = b1[p[i]];
		for (int j = 1; j >= i + 1; j--)
			x1[i] -= m1[p[i]][j] * x1[j];
	}
}

double sqo(double x[7], double y[7], double a[2], int N, int m)
{
	double S = 0, s;
	for (int i = 0; i < N; i++)
	{
		s = y[i];
		for (int j = 0; j < m + 1; j++)
			s -= a[j] * pow(x[i], j);
		S += s * s;
	}
	S /= N - m - 1;
	return sqrt(S);
}