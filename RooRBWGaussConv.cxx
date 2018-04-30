//
// Compile and load the class using:
// root> .L RooRBWGaussConv.cxx+
//

#include "RooFit.h"

#include "Riostream.h"
#include <math.h>

#include "RooMath.h"
#include "RooRealConstant.h"
#include "RooRealVar.h"
#include "TMath.h"
#include "TRegexp.h"

#include "RooRBWGaussConv.h"

ClassImp(RooRBWGaussConv)

    double RooRBWGaussConv::D0MASS = 1.8645;
double RooRBWGaussConv::MPI        = 0.13957;

RooRBWGaussConv::RooRBWGaussConv(const char *name,
                                 const char *title,
                                 RooAbsReal &_x,
                                 RooAbsReal &_bwmean,
                                 RooAbsReal &_bwwidth,
                                 RooAbsReal &_gmean,
                                 RooAbsReal &_gsigma,
                                 MassType _massType)
    : RooAbsPdf(name, title)
    , x("x", "Dependent", this, _x)
    , bw_mean("bw_mean", "BW Mean", this, _bwmean)
    , bw_width("bw_width", "BW Width", this, _bwwidth)
    , g_mean("g_mean", "Gaussian mean", this, _gmean)
    , g_sigma("g_sigma", "Gaussian width", this, _gsigma)
    , massType(_massType) {}

RooRBWGaussConv::RooRBWGaussConv(const RooRBWGaussConv &other, const char *name)
    : RooAbsPdf(other, name)
    , x("x", this, other.x)
    , bw_mean("bw_mean", this, other.bw_mean)
    , bw_width("bw_width", this, other.bw_width)
    , g_mean("g_mean", this, other.g_mean)
    , g_sigma("g_sigma", this, other.g_sigma)
    , massType(other.massType) {}

double RooRBWGaussConv::breitwigner(double val) const {
    double newGam = bw_mean * bw_width;
    double newM   = val * val - bw_mean * bw_mean;
    return (newGam) / (newM * newM + newGam * newGam);
}

double RooRBWGaussConv::gaussian(double val) const {
    Double_t arg = val - g_mean;
    return exp(-0.5 * arg * arg / (g_sigma * g_sigma));
}
// Blatt-Weisskopf form-factor
// Treat as square well potential where r is effective radius in which teh D-pi scatter,
// p is momentum of either daugther in Dpi center of mass and l=J is ang momentum of Dpi cm.
//

double RooRBWGaussConv::evalbwff(double pval) const {
    static double r = 1.6; // for spin 1 particle r=1.6 (GeV/c)^-1 (~3fm)
    return 1.0 / sqrt(1.0 + r * r * pval * pval);
    // return 1.0;
}

double RooRBWGaussConv::evalMass(double val) const { return D0MASS + val; }
double RooRBWGaussConv::evalphsp(double val) const {
    // process width Gamma_(D*+ -> D0 pi+) = phspfactor * Gamma_0
    //    ==> we assume Gamma_total = Gamma_(D*+->D0pi+)
    double realM0;
    double realM;

    if(massType == DeltaMass) {
        realM0 = evalMass(bw_mean);
        realM  = evalMass(val);
    } else {
        realM0 = bw_mean;
        realM  = val;
    }

    // make sure denominator in phspfactor not zero (evalbwff never 0)
    if(evalP(realM0) == 0)
        return 0.0;

    double phspfactor = pow(evalP(realM) / evalP(realM0), 3) * pow(evalbwff(evalP(realM)) / evalbwff(evalP(realM0)), 2);
    double phspM      = realM0 * realM0 - realM * realM; // m_0^2 - m^2
    double phspGam    = bw_width * phspfactor;

    //  phspfactor * (m0 Gam_0)^2 / ((m0^2-m^2)^2+(m0 Gam_Total)^2)
    return (phspfactor * realM0 * realM0 * bw_width * bw_width) / (phspM * phspM + realM0 * phspGam * realM0 * phspGam);
}
Double_t RooRBWGaussConv::evalP(double mass) const {
    // masses are passed into here **not** delta m
    double lambda = (mass * mass - MPI * MPI - D0MASS * D0MASS) * (mass * mass - MPI * MPI - D0MASS * D0MASS)
                    - 4 * MPI * MPI * D0MASS * D0MASS;
    if(lambda > 0) {
        return sqrt(lambda) / (2 * mass);
    } else {
        return 0;
    }
}
// To evaluate numeric method
Double_t RooRBWGaussConv::evaluate() const {
    double interval = 10 * std::max(bw_width, g_sigma);
    double step     = interval * 0.001;
    double ret      = 0;
    for(double val = x - interval; val < x + interval; val += step) {
        // std::cout<<"at val ("<<val<<") : ";
        // std::cout<<"G: "<<gaussian(val-x)<<"  |  ";
        // std::cout<<"BW: "<<evalphsp(val)<<std::endl;
        ret += evalphsp(val) * gaussian(val - x);
    }
    return ret;
}
