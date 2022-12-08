#include "experiment.h"
#include "sbml/SBMLTypes.h"
#include <stdio.h>
#define rep(i, n) for (int i = 0; i < n; i++)
#define max(p,q)((p)>(q)?(p):(q))
const int INF = 1001001001;
const int N = 1000, M = 1000000;
const char *IDs1 = "s1", *IDs2 = "s2";
double dk = 0.01, dt = 0.01, k_max = 1.01, k = 0;
double S1_SIM[M], S2_SIM[M], time[M];

void init(double initS1, double initS2) {
	S1_SIM[0] = initS1, S2_SIM[0] = initS2, time[0] = 0;
	return;
}

double f(double k, double x) { return k*x; }

void RungeKutta4D(double k) {
  rep(i,N) {
    time[i+1] = time[i]+dt;
    double d1 = dt * f(k,S1_SIM[i]);
    double d2 = dt * f(k, S1_SIM[i] - d1 / 2.0);
    double d3 = dt * f(k, S1_SIM[i] - d2 / 2.0);
		double d4 = dt * f(k, S1_SIM[i] - d3);
		double dif = (d1 + 2.0 * d2 + 2.0 * d3 + d4) / 6.0;
    S1_SIM[i+1] = S1_SIM[i]-dif;
    S2_SIM[i+1] = S2_SIM[i]+dif;
  }
}

void DumpSimulate(FILE *fp) { rep(i,N+1) { fprintf(fp,"%.10lf,%.10lf\n", (double)i*dt,S1_SIM[i]); } }

void DumpReal(FILE *fp) { rep(i,N+1) { fprintf(fp,"%.10lf,%.10lf\n", (double)i*dt, S1[i]); } }

double likelihood_simple(double initS1, double initS2, double k) {
	init(initS1, initS2);
  RungeKutta4D(k);
  double likelihood = 0;
  //10/0.5 = 20
  //0.5秒後の実データとシミュレートの実データで尤度を計算
  rep(i,21) likelihood -= pow((S1[i]-S1_SIM[i*50]),2);
  return likelihood;
}

int main() {
	SBMLDocument_t *d = readSBML("simple.xml");
	Model_t *m = SBMLDocument_getModel(d);
	// ListOf_t *ListOfSpecies = Model_getListOfSpecies(m);
    FILE *fp1 = fopen("res.csv", "w");
    FILE *fp2 = fopen("bestK_res.csv", "w");
    FILE *fp3 = fopen("real.csv", "w");
	Species_t *s1 = Model_getSpeciesById(m, IDs1);
	Species_t *s2 = Model_getSpeciesById(m, IDs2);

	double initS1 = Species_getInitialAmount(s1);
	double initS2 = Species_getInitialAmount(s2);
    double bestK = -1;
    double tmp = -INF;
	while (k <= k_max) {
    double likelihood = likelihood_simple(initS1, initS2, k);
    fprintf(fp1, "%.10lf,%.10lf\n", k, likelihood);
    if(max(tmp,likelihood) == likelihood) { bestK = k, tmp = likelihood; }
		k += dk;
	}
    printf("The best k is %.10lf\n", bestK);
    init(initS1, initS2);
    //variable K' Simulate & DumpRes
    RungeKutta4D(bestK);
    DumpSimulate(fp2);
    DumpReal(fp3);
	SBMLDocument_free(d);
	return 0;
}
