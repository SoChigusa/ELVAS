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
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

// physical parameters as of 2023.6
#define DATE_LABEL "202306"
#define ALPHAS 0.1179
#define MTPOLE 172.5
#define MW 80.377
#define MH 125.25

// physical parameters used for [1803.03902]
// just for check
// #define ALPHAS 0.1181
// #define MTPOLE 173.1
// #define MW 80.377
// #define MH 125.09

using namespace std;

double calcLog10gamma(double arg_alphas, double arg_mtpole, double arg_mw, double arg_mh, const string &arg_ofprefix = "")
{
  // SM parameters @ Mt
  // use three-loop RGE by default
  StandardModel sm;
  QedQcd qq;
  qq.setAlphas(arg_alphas);
  qq.setPoleMt(arg_mtpole);
  qq.setMW(arg_mw);
  qq.toMt();
  sm.matchQedQcd(qq, arg_mw, arg_mtpole, arg_mh, arg_alphas);

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

  // absolute stability
  if (nI == -1)
    return -DBL_MAX;

  // decay rate calculation
  int n, mI = -1, mF = -1;
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
    {
      lndgamdRinv = -B0 - B1 + 4 * log(vec_mu[n]);
      if (mI == -1)
        mI = i;
    }
    else
    {
      lndgamdRinv = -DBL_MAX;
      if (mI != -1 && mF == -1)
        mF = i - 1;
    }
    lndgam[i] = pair<double, double>(lnRinv, lndgamdRinv);
  }

  // absolute stability
  // (near criticality such that no integration range obtained)
  if (mI == -1)
    return -DBL_MAX;

  double lngamma = Elvas::getLnGamma(lndgam, log(vec_mu[mI + nI]), log(vec_mu[mF + nI]));

  // save RGE data
  if (arg_ofprefix != "")
  {
    // RGE flow
    ofstream ofs("output/" + arg_ofprefix + "_RGEFlow.dat");
    ofs << scientific << setprecision(6);
    for (int i = 0; i <= npts; ++i)
    {
      ofs << vec_mu[i] << "\t" << vec_lambda[i] << endl;
    }

    // output differential rate
    ofstream ofs2("output/" + arg_ofprefix + "_differential_rate.dat");
    ofs2 << scientific << setprecision(6);
    for (int i = 0; i < nF - nI; ++i)
    {
      n = i + nI;
      ofs2 << exp(lndgam[i].first) << "\t" << lndgam[i].second << endl;
    }
  }

  // log10(gamma / Gyr / Gpc^3)
  double log10GeV4inGyrGpc3 = log10(1.83) + 164;
  return lngamma / log(10) + log10GeV4inGyrGpc3;
}

int main(int argc, char **argv)
{
  // current central value
  calcLog10gamma(ALPHAS, MTPOLE, MW, MH, DATE_LABEL);

  // ----- contour plot -----
  // grid setting
  int nPts = 40;
  double mt_min = 170.;
  double mt_max = 180.;
  double mh_min = 120.;
  double mh_max = 130.;
  double dmt = (mt_max - mt_min) / nPts;
  double dmh = (mh_max - mh_min) / nPts;

  // run and save
  ofstream ofs("output/" + string(DATE_LABEL) + ".dat");
  ofs << scientific << setprecision(6);
  double mt, mh, log10gamma;
  for (int i = 0; i <= nPts; ++i)
  {
    mt = mt_min + dmt * i;
    cout << "start running mt = " << mt << " GeV..." << endl;
    for (int j = 0; j <= nPts; ++j)
    {
      mh = mh_min + dmh * j;
      log10gamma = calcLog10gamma(ALPHAS, mt, MW, mh);
      ofs << mt << "\t" << mh << "\t" << log10gamma << endl;
    }
  }

  return 0;
}
