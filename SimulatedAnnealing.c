#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include <stdbool.h>
#include"sbml/SBMLTypes.h"
#include"experiment.h"
#define rep(i, n) for (int i = 0; i < n; i++)
#define max(p,q)((p)>(q)?(p):(q))
double unif_param(void);
double unif_update(void);
#define r unif_param()
#define pp unif_update()
const int INF = 1001001001;
const double DINF = 1001001001;
const int N = 1000, M = 1000000;
const char *IDs1 = "s1", *IDs2 = "s2";
double dk = 0.01, dt = 0.01, k_max = 1.00, k = 0;
double S1_SIM[M], S2_SIM[M], time[M];

void init(double initS1, double initS2) {
	S1_SIM[0] = initS1, S2_SIM[0] = initS2, time[0] = 0;
	return;
}

double f(double k, double x) { return k*x; }

void RungeKutta4D(double k) {
    rep(i,N) {
        time[i+1] = time[i]+dt;
        double d1 = dt * f(k, S1_SIM[i]);
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

int main(void){
    srand(100); /* 乱数のシード設定 */
    double k = 0.8;
    double T = 1000;
    double alpha = 0.01;
    int Nt = 10;
    double rt = 0.8;
    double epsilon = pow(10,-8);
    double difference; /* 更新前後の尤度の差を格納する */
    /* simple.xmlを読み込み最急勾配法に必要な情報を取得する */
    SBMLDocument_t *d = readSBML("simple.xml");
    Model_t *m = SBMLDocument_getModel(d);
    FILE *fp = fopen("res.csv", "w");
	Species_t *s1 = Model_getSpeciesById(m, IDs1);
	Species_t *s2 = Model_getSpeciesById(m, IDs2);
	double initS1 = Species_getInitialAmount(s1);
	double initS2 = Species_getInitialAmount(s2);

    double l_ans = -DINF, k_ans = 0, max_l = 0, min_l = 0, k_lmax = 0;

    while(1) {
        int cnt = Nt;
        bool initialized = false;
        while(1) {
            double kp = k+alpha*r;
            double l = likelihood_simple(initS1, initS2, k);
            double lp = likelihood_simple(initS1, initS2, kp);
      
            if(!initialized) {
                max_l = -DINF, min_l = DINF;
                initialized = true;
            }
            double p = exp((lp-l)/T);
            if(p >= pp) {
                k = kp;
                cnt--;        
                if(lp > max_l) { max_l = lp, k_lmax = kp; }
                if(lp < min_l) { min_l = lp; }
                fprintf(fp,"%.10lf,%.10lf\n",kp,lp);
            } 
            if(cnt == 0) break;
        }
        if(l_ans < max_l) { l_ans = max_l, k_ans = k_lmax; }
        difference = fabs(max_l - min_l);
        if(difference <= epsilon) break;
        T *= rt;
    }
    printf("The best K is %.10lf, linklihood = %.10lf\n", k_ans, l_ans);
    fprintf(fp,"The best K is %.10lf, linklihood = %.10lf\n", k_ans, l_ans);
    SBMLDocument_free(d);
    return 0;
}
