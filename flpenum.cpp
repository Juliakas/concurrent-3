#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <mpi.h>
#include <unistd.h>

using namespace std;

int numDP = 1000; // Vietoviu skaicius (demand points, max 10000)
int numPF = 10;	  // Esanciu objektu skaicius (preexisting facilities)
int numCL = 25;	  // Kandidatu naujiems objektams skaicius (candidate locations)
int numX = 5;	  // Nauju objektu skaicius

double **demandPoints; // Geografiniai duomenys

//=============================================================================

double getTime();
void loadDemandPoints();
double HaversineDistance(double *a, double *b);
double evaluateSolution(int *X);
int increaseX(int *X, int index, int maxindex);
int kthCombination(int *choices, int k, int r, int n, int *result);
int factorial(int n, int limit);
int getnCr(int n, int r);

//=============================================================================

int main(int argc, char **argv)
{
	double ts = getTime(); // Algoritmo vykdymo pradzios laikas

	loadDemandPoints(); // Nuskaitomi duomenys

	// Sudarom pradini sprendini: [0, 1, 2, 3, ...]
	int *X = new int[numX];
	int *bestX = new int[numX];
	int *choices = new int[numCL];
	int rank;
	int size;
	if (argc > 1)
	{
		numX = stoi(argv[1]);
	}
	for (int i = 0; i < numX; i++)
	{
		X[i] = i;
		bestX[i] = i;
	}
	for (int i = 0; i < numCL; i++)
	{
		choices[i] = i;
	}

	// n = numCL, r = numX
	int nCr = getnCr(numCL, numX);
	double tp = getTime();
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	int current = rank * (nCr / size);
	int end; // Non inclusive
	if (rank == size - 1)
	{
		end = nCr - 1;
	}
	else
	{
		end = (rank + 1) * (nCr / size);
	}
	kthCombination(choices, current++, numX, numCL, X);
	int tag = 0;
	double u = evaluateSolution(X);
	double bestU = u;
	int r;

	//----- Pagrindinis ciklas ------------------------------------------------
	while (current++ < end)
	{
		increaseX(X, numX - 1, numCL);
		u = evaluateSolution(X);

		if (u > bestU)
		{
			bestU = u;
			for (int i = 0; i < numX; i++)
				bestX[i] = X[i];
		}
	}

	double UBuf[1];
	int *XBuf = new int[numX];
	if (rank == 0)
	{
		for (int i = 1; i < size; i++)
		{
			MPI_Status status;
			MPI_Recv(UBuf, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			MPI_Recv(XBuf, numX, MPI_INT, status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if (UBuf[0] > bestU)
			{
				bestU = UBuf[0];
				for (int i = 0; i < numX; i++)
				{
					X[i] = XBuf[i];
				}
			}
		}
	}
	else
	{
		for (int i = 0; i < numX; i++)
		{
			XBuf[i] = X[i];
		}
		MPI_Send(UBuf, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
		MPI_Send(XBuf, numX, MPI_INT, 0, rank, MPI_COMM_WORLD);
	}

	delete XBuf;

	//----- Rezultatu spausdinimas --------------------------------------------
	if (rank == 0)
	{
		double tf = getTime(); // Skaiciavimu pabaigos laikas

		cout << "Geriausias sprendinys: ";
		for (int i = 0; i < numX; i++)
			cout << bestX[i] << " ";
		cout << "(" << bestU << ")" << endl;
		cout << "Skaiciavimo trukme: " << tf - ts << endl;
		cout << "Nuoseklioji dalis: " << tp - ts << endl;
		cout << "Lygiagrecioji dalis: " << tf - tp << endl;
	}
	MPI_Finalize();
}

//=============================================================================

void loadDemandPoints()
{
	FILE *f;
	f = fopen("demandPoints.dat", "r");
	demandPoints = new double *[numDP];
	for (int i = 0; i < numDP; i++)
	{
		demandPoints[i] = new double[3];
		fscanf(f, "%lf%lf%lf", &demandPoints[i][0], &demandPoints[i][1], &demandPoints[i][2]);
	}
	fclose(f);
}

//=============================================================================

double HaversineDistance(double *a, double *b)
{
	double dlon = fabs(a[0] - b[0]);
	double dlat = fabs(a[1] - b[1]);
	double aa = pow((sin((double)dlon / (double)2 * 0.01745)), 2) + cos(a[0] * 0.01745) * cos(b[0] * 0.01745) * pow((sin((double)dlat / (double)2 * 0.01745)), 2);
	double c = 2 * atan2(sqrt(aa), sqrt(1 - aa));
	double d = 6371 * c;
	return d;
}

//=============================================================================

double getTime()
{
	struct timeval laikas;
	gettimeofday(&laikas, NULL);
	double rez = (double)laikas.tv_sec + (double)laikas.tv_usec / 1000000;
	return rez;
}

//=============================================================================

double evaluateSolution(int *X)
{
	double U = 0;
	int bestPF;
	int bestX;
	double d;
	for (int i = 0; i < numDP; i++)
	{
		bestPF = 1e5;
		for (int j = 0; j < numPF; j++)
		{
			d = HaversineDistance(demandPoints[i], demandPoints[j]);
			if (d < bestPF)
				bestPF = d;
		}
		bestX = 1e5;
		for (int j = 0; j < numX; j++)
		{
			d = HaversineDistance(demandPoints[i], demandPoints[X[j]]);
			if (d < bestX)
				bestX = d;
		}
		if (bestX < bestPF)
			U += demandPoints[i][2];
		else if (bestX == bestPF)
			U += 0.3 * demandPoints[i][2];
	}
	return U;
}

//=============================================================================

int increaseX(int *X, int index, int maxindex)
{
	if (X[index] + 1 < maxindex - (numX - index - 1))
	{
		X[index]++;
	}
	else
	{
		if ((index == 0) && (X[index] + 1 == maxindex - (numX - index - 1)))
		{
			return 0;
		}
		else
		{
			if (increaseX(X, index - 1, maxindex))
				X[index] = X[index - 1] + 1;
			else
				return 0;
		}
	}
	return 1;
}

int kthCombination(int *choices, int k, int r, int n, int *result)
{
	if (r == 0)
	{
		return 0;
	}
	else if (n == r)
	{
		for (int i = 0; i < n; i++)
		{
			result[i] = choices[i];
		}
		return n;
	}
	else
	{
		int nCr = getnCr(n - 1, r - 1);
		if (k < nCr)
		{
			int *tmpResult = new int[r - 1];
			kthCombination(choices + 1, k, r - 1, n - 1, tmpResult);
			result[0] = choices[0];
			for (int i = 1; i < r; i++)
			{
				result[i] = tmpResult[i - 1];
			}
			delete tmpResult;
			return r;
		}
		else
		{
			kthCombination(choices + 1, k - nCr, r, n - 1, result);
			return r;
		}
	}
}

int factorial(int n, int limit)
{
	if (limit == 0 || n == 0)
		return 1;
	return n * factorial(n - 1, limit - 1);
}

int getnCr(int n, int r)
{
	int limit = min(r, n - r);
	return factorial(n, limit) / factorial(limit, limit);
}