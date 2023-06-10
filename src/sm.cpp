/**
 * @file sm.cpp
 * @brief SM phase diagram
 * @author So Chigusa (Lawrence Berkeley National Laboratory)
 * @date Created on: 2023/06/07
 */

#include "elvas.h"
#include "qedqcd.h"
#include "sm.h"
#include <cfloat>
#include <fstream>
#include <iostream>
#include <sstream>

// physical parameters 2023.6
// #define ALPHAS 0.1179
// #define MTPOLE 172.5
// #define MW 80.377
// #define MH 125.25
#define ALPHAS 0.1181
#define MTPOLE 173.1
#define MW 80.377
#define MH 125.09

using namespace std;

int main(int argc, char **argv)
{
  // SM parameters @ Mt
  StandardModel sm;
  QedQcd qq;
  qq.setAlphas(ALPHAS);
  qq.setPoleMt(MTPOLE);
  qq.setMW(MW);
  qq.toMt();
  // sm.flagThreeLoop(false);
  sm.matchQedQcd(qq, MW, MTPOLE, MH, ALPHAS);

  // RGE flow settings
  int npts = 1000;
  vector<double> vec_mu(npts + 1), vec_lambda(npts + 1), vec_yt(npts + 1), vec_gW(npts + 1), vec_gZ(npts + 1);
  double muI = sm.displayMu();
  double muF = 1.e40;
  double logmuI = log(muI);
  double logmuF = log(muF);
  double dlogmu = (logmuF - logmuI) / npts;

  // take RGE flow data
  int nI = -1, nF = -1;
  double mu;
  for (int i = 0; i <= npts; ++i)
  {
    mu = exp(logmuI + i * dlogmu);
    sm.runto(mu);

    double g1 = sm.displayGaugeElement(1);
    double g2 = sm.displayGaugeElement(2);
    vec_mu[i] = mu;
    vec_lambda[i] = sm.displayLambda();
    vec_yt[i] = sm.displayYukawaElement(yukawaID::YU, 3, 3);
    vec_gW[i] = 0.5 * g2;
    vec_gZ[i] = 0.5 * sqrt(0.6 * g1 * g1 + g2 * g2);
    if (nI == -1 && sm.displayLambda() < 0)
      nI = i;
    if (nI >= 0 && nF == -1 && sm.displayLambda() > 0)
      nF = i;
  }

  // output RGE flow
  ofstream ofs("output/RGEFlow.dat");
  for (int i = 0; i <= npts; ++i)
  {
    ofs << vec_mu[i] << "\t" << vec_lambda[i] << endl;
  }

  // decay rate calculation
  int n;
  vector<pair<double, double>> lndgam(nF - nI);
  for (int i = 0; i < nF - nI; ++i)
  {
    n = i + nI;
    double lambdaAbs = fabs(vec_lambda[n]);
    double B0 = Elvas::instantonB(lambdaAbs);
    double mlnAh = Elvas::higgsQC(lambdaAbs, 0.);
    double mlnAt = Elvas::fermionQC(vec_yt[n], lambdaAbs, 0.);
    double mlnAgW = Elvas::gaugeQC(vec_gW[n] * vec_gW[n], lambdaAbs, 0.);
    double mlnAgZ = Elvas::gaugeQC(vec_gZ[n] * vec_gZ[n], lambdaAbs, 0.);
    double B1 = mlnAh + 3 * mlnAt - log(2 * M_PI * M_PI) + 2 * mlnAgW + mlnAgZ;

    double lndgamdRinv;
    double lnRinv = log(vec_mu[n]); // take mu = 1/R
    bool pert1 = (fabs(B1 / B0) < 0.8);
    bool perth = (fabs(mlnAh / B0) < 0.8);
    bool pertt = (fabs(3 * mlnAt / B0) < 0.8);
    bool pertgW = (fabs(2 * mlnAgW / B0) < 0.8);
    bool pertgZ = (fabs(mlnAgZ / B0) < 0.8);
    if (pert1 && perth && pertt && pertgW && pertgZ)
      lndgamdRinv = -B0 - B1 + 4 * log(vec_mu[n]);
    else
      lndgamdRinv = -DBL_MAX;
    lndgam[i] = pair<double, double>(lnRinv, lndgamdRinv);
  }
  cout << vec_mu[nI] << "\t" << vec_mu[nF] << "\t" << endl;
  double lngamma = Elvas::getLnGamma(lndgam, log(8.5e10), log(3e28));
  cout << lngamma / log(10) + (log10(1.83) + 164) << endl;

  // output differential rate
  ofstream ofs2("output/differential_rate.dat");
  for (int i = 0; i < nF - nI; ++i)
  {
    n = i + nI;
    ofs2 << exp(lndgam[i].first) << "\t" << lndgam[i].second << endl;
  }

  return 0;
}
