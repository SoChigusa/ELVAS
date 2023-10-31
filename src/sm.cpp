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
#define DATE_LABEL "2306"
#define ALPHAS 0.1179
#define ALPHAS_ERR 0.0009
#define MTPOLE 172.69 // MC mass
#define MT_ERR 0.3
#define MW 80.377
#define MH 125.25
#define MH_ERR 0.17
int nPts = 100;

// // physical parameters used for [1803.03902]
// #define DATE_LABEL "1803"
// #define ALPHAS 0.1181
// #define ALPHAS_ERR 0.0011
// #define MTPOLE 173.1 // MC mass
// #define MT_ERR 0.6
// #define MW 80.379
// #define MH 125.09
// #define MH_ERR 0.24
// int nPts = 50;

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

double calcLog10gamma(double arg_alphas, double arg_mtpole, double arg_mw, double arg_mh,
                      ofstream &arg_ofs, bool arg_save = true)
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
  int npts = 189; // just a convention
  vector<double> vec_mu(npts + 1), vec_lambda(npts + 1), vec_yt(npts + 1), vec_yb(npts + 1), vec_g2(npts + 1), vec_g1(npts + 1);
  double muI = 240.;          // just a convention
  double muF = 1.5142976e+40; // just a convention
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

    vec_mu[i] = mu;
    vec_lambda[i] = sm.displayLambda();
    vec_yt[i] = sm.displayYukawaElement(yukawaID::YU, 3, 3);
    vec_yb[i] = sm.displayYukawaElement(yukawaID::YD, 3, 3);
    vec_g2[i] = sm.displayGaugeElement(2);
    vec_g1[i] = sm.displayGaugeElement(1);
    if (nI == -1 && sm.displayLambda() < 0)
      nI = i;
    if (nI >= 0 && nF == -1 && sm.displayLambda() > 0)
      nF = i;
  }

  // save RGE data
  if (arg_save)
  {
    arg_ofs << scientific << setprecision(7);

    // ELVAS input header
    arg_ofs << "[DATASET] (" << arg_mh << " " << arg_mtpole << ")" << endl;

    // RGE flow
    for (int i = 0; i <= npts; ++i)
    {
      // {Q, g2, g1, yt, yb, lambda}
      arg_ofs << vec_mu[i] << "\t" << vec_g2[i] << "\t" << vec_g1[i] << "\t" << vec_yt[i] << "\t" << vec_yb[i] << "\t" << vec_lambda[i] << endl;
    }

    // ELVAS input footer
    arg_ofs << endl;

    // output differential rate
    // ofstream ofs2("output/" + arg_ofprefix + "_differential_rate.dat");
    // ofs2 << scientific << setprecision(6);
    // for (auto itr = lndgam.begin(); itr != lndgam.end(); ++itr)
    // {
    //   ofs2 << exp(itr->first) << "\t" << itr->second << endl;
    // }
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
    double mlnAgW = Elvas::gaugeQC(0.25 * vec_g2[n] * vec_g2[n], lambdaAbs, 0.);
    double mlnAgZ = Elvas::gaugeQC(0.25 * (0.6 * vec_g1[i] * vec_g1[i] + vec_g2[i] * vec_g2[i]), lambdaAbs, 0.);
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

  // log10(gamma / Gyr / Gpc^3)
  double log10GeV4inGyrGpc3 = log10(1.83) + 164;
  return lngamma / log(10) + log10GeV4inGyrGpc3;
}

double calcLog10gamma(double arg_alphas, double arg_mtpole, double arg_mw, double arg_mh)
{
  ofstream ofs_tmp("tmp");
  return calcLog10gamma(arg_alphas, arg_mtpole, arg_mw, arg_mh, ofs_tmp, false);
}

int main(int argc, char **argv)
{
  // current value and error bars
  ofstream ofs_center_ELVAS("output/" + string(DATE_LABEL) + "_center_ELVAS_input.dat");
  double center = calcLog10gamma(ALPHAS, MTPOLE, MW, MH, ofs_center_ELVAS);
  double mh_plus = calcLog10gamma(ALPHAS, MTPOLE, MW, MH + MH_ERR, ofs_center_ELVAS);
  double mh_minus = calcLog10gamma(ALPHAS, MTPOLE, MW, MH - MH_ERR, ofs_center_ELVAS);
  double mt_plus = calcLog10gamma(ALPHAS, MTPOLE + MT_ERR, MW, MH, ofs_center_ELVAS);
  double mt_minus = calcLog10gamma(ALPHAS, MTPOLE - MT_ERR, MW, MH, ofs_center_ELVAS);
  double alphas_plus = calcLog10gamma(ALPHAS + ALPHAS_ERR, MTPOLE, MW, MH, ofs_center_ELVAS);
  double alphas_minus = calcLog10gamma(ALPHAS - ALPHAS_ERR, MTPOLE, MW, MH, ofs_center_ELVAS);
  cout << "log_10 gamma = " << int(center)
       << " +" << int(mh_minus - center) << "-" << int(center - mh_plus)
       << " +" << int(mt_plus - center) << "-" << int(center - mt_minus)
       << " +" << int(alphas_minus - center) << "-" << int(center - alphas_plus) << endl;

  // for point specific debug
  // cout << calcLog10gamma(ALPHAS, 173.50, MW, 129.75) << endl;

  // ----- contour plot -----
  // grid setting
  double dmt = (MT_MAX - MT_MIN) / nPts;
  double dmh = (MH_MAX - MH_MIN) / nPts;

  // run and save
  ofstream ofs("output/" + string(DATE_LABEL) + ".dat");
  ofstream ofs_ELVAS("output/" + string(DATE_LABEL) + "_ELVAS_input.dat");
  ofstream ofs_ELVAS_alphap1s("output/" + string(DATE_LABEL) + "_alphap1s_ELVAS_input.dat");
  ofstream ofs_ELVAS_alpham1s("output/" + string(DATE_LABEL) + "_alpham1s_ELVAS_input.dat");
  ofs << scientific << setprecision(7);
  double mt, mh;
  vector<double> log10gamma(3);
  for (int i = 0; i <= nPts; ++i)
  {
    mt = MT_MIN + dmt * i;
    cout << "start running mt = " << mt << " GeV..." << endl;
    for (int j = 0; j <= nPts; ++j)
    {
      mh = MH_MIN + dmh * j;

      // calculate decay rate
      log10gamma[0] = calcLog10gamma(ALPHAS, mt, MW, mh, ofs_ELVAS);
      log10gamma[1] = calcLog10gamma(ALPHAS + ALPHAS_ERR, mt, MW, mh, ofs_ELVAS_alphap1s);
      log10gamma[2] = calcLog10gamma(ALPHAS - ALPHAS_ERR, mt, MW, mh, ofs_ELVAS_alpham1s);
      ofs << mt << "\t" << mh << "\t" << log10gamma[0] << "\t" << log10gamma[1] << "\t" << log10gamma[2] << endl;
    }
  }

  return 0;
}
