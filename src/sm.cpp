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
// #define DATE_LABEL "202306"
// #define ALPHAS 0.1179
// #define ALPHAS_ERR 0.0009
// #define MTPOLE 172.69 // MC mass
// #define MT_ERR 0.3
// #define MW 80.377
// #define MH 125.25
// #define MH_ERR 0.17

// physical parameters used for [1803.03902]
#define DATE_LABEL "201803"
#define ALPHAS 0.1181
#define ALPHAS_ERR 0.0011
#define MTPOLE 173.1 // MC mass
#define MT_ERR 0.6
#define MW 80.379
#define MH 125.09
#define MH_ERR 0.24

// contour plot settings
#define MT_MIN 170.
#define MT_MAX 180.
#define MH_MIN 120.
#define MH_MAX 130.

using namespace std;

double fakeRateAbsoluteStability(double mh, double mt)
{
  // rate varies from norm -- 3*norm
  // double norm = -1e4;
  // return norm * (1 + (mh - MH_MIN) / (MH_MAX - MH_MIN) + (MT_MAX - mt) / (MT_MAX - MT_MIN));
  return -DBL_MAX;
}

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
    return fakeRateAbsoluteStability(arg_mh, arg_mtpole);

  // decay rate calculation
  int n;
  muI = muF = 0.;
  vector<pair<double, double>> lndgam;
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
    double phiC = vec_mu[n] * sqrt(8 / lambdaAbs);
    double cutoff_phiC = 2.4e18;
    double thresh = 0.8;
    bool pert1 = (fabs(B1 / B0) < thresh);
    bool perth = (fabs(mlnAh / B0) < thresh);
    bool pertt = (fabs(3 * mlnAt / B0) < thresh);
    bool pertgW = (fabs(2 * mlnAgW / B0) < thresh);
    bool pertgZ = (fabs(mlnAgZ / B0) < thresh);
    if (pert1 && perth && pertt && pertgW && pertgZ)
    {
      // cutoff the integration at \bar{phi}_C = Mpl
      if (phiC < cutoff_phiC)
      {
        lndgamdRinv = -B0 - B1 + 4 * log(vec_mu[n]);
        lndgam.emplace_back(pair<double, double>(lnRinv, lndgamdRinv));

        if (muI == 0.)
          muI = vec_mu[n];
      }
      else if (muF == 0.)
        muF = vec_mu[n - 1];
    }
  }

  // absolute stability
  // (near criticality such that no integration range obtained)
  if (lndgam.size() < 3)
    return fakeRateAbsoluteStability(arg_mh, arg_mtpole);

  // for debug of integration
  // cout << nI << "\t" << nF << "\t" << muI << "\t" << muF << endl;
  // cout << exp(lndgam[0].first) << "\t" << lndgam[0].second << endl;
  // cout << exp(lndgam[-1].first) << "\t" << lndgam[-1].second << endl;

  // integration
  // +- 3 to mI/mF to ensure the smooth interpolation within the range
  double lngamma = Elvas::getLnGamma(lndgam, log(muI), log(muF));

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
    for (auto itr = lndgam.begin(); itr != lndgam.end(); ++itr)
    {
      ofs2 << exp(itr->first) << "\t" << itr->second << endl;
    }
  }

  // log10(gamma / Gyr / Gpc^3)
  double log10GeV4inGyrGpc3 = log10(1.83) + 164;
  return lngamma / log(10) + log10GeV4inGyrGpc3;
}

int main(int argc, char **argv)
{
  // current value and error bars
  double center = calcLog10gamma(ALPHAS, MTPOLE, MW, MH, DATE_LABEL);
  double mh_plus = calcLog10gamma(ALPHAS, MTPOLE, MW, MH + MH_ERR, DATE_LABEL);
  double mh_minus = calcLog10gamma(ALPHAS, MTPOLE, MW, MH - MH_ERR, DATE_LABEL);
  double mt_plus = calcLog10gamma(ALPHAS, MTPOLE + MT_ERR, MW, MH, DATE_LABEL);
  double mt_minus = calcLog10gamma(ALPHAS, MTPOLE - MT_ERR, MW, MH, DATE_LABEL);
  double alphas_plus = calcLog10gamma(ALPHAS + ALPHAS_ERR, MTPOLE, MW, MH, DATE_LABEL);
  double alphas_minus = calcLog10gamma(ALPHAS - ALPHAS_ERR, MTPOLE, MW, MH, DATE_LABEL);
  cout << "log_10 gamma = " << int(center)
       << " +" << int(mh_minus - center) << "-" << int(center - mh_plus)
       << " +" << int(mt_plus - center) << "-" << int(center - mt_minus)
       << " +" << int(alphas_minus - center) << "-" << int(center - alphas_plus) << endl;

  // for point specific debug
  // cout << calcLog10gamma(ALPHAS, 173.50, MW, 129.75, DATE_LABEL) << endl;

  // ----- contour plot -----
  // grid setting
  int nPts = 300;
  double dmt = (MT_MAX - MT_MIN) / nPts;
  double dmh = (MH_MAX - MH_MIN) / nPts;

  // run and save
  ofstream ofs("output/" + string(DATE_LABEL) + ".dat");
  ofs << scientific << setprecision(6);
  double mt, mh;
  vector<double> log10gamma(3);
  for (int i = 0; i <= nPts; ++i)
  {
    mt = MT_MIN + dmt * i;
    cout << "start running mt = " << mt << " GeV..." << endl;
    for (int j = 0; j <= nPts; ++j)
    {
      mh = MH_MIN + dmh * j;
      log10gamma[0] = calcLog10gamma(ALPHAS, mt, MW, mh);
      log10gamma[1] = calcLog10gamma(ALPHAS + ALPHAS_ERR, mt, MW, mh);
      log10gamma[2] = calcLog10gamma(ALPHAS - ALPHAS_ERR, mt, MW, mh);
      ofs << mt << "\t" << mh << "\t" << log10gamma[0] << "\t" << log10gamma[1] << "\t" << log10gamma[2] << endl;
    }
  }

  return 0;
}
