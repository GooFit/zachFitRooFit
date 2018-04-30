
#ifndef ROORBWGAUSSCONV
#define ROORBWGAUSSCONV

#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooRealProxy.h"

// class RooRealVar;
// class RooAbsReal;

class RooRBWGaussConv : public RooAbsPdf {
  public:
    enum MassType { DeltaMass, Mass };
    RooRBWGaussConv(){};
    RooRBWGaussConv(const char *name,
                    const char *title,
                    RooAbsReal &_x,
                    RooAbsReal &_bwmean,
                    RooAbsReal &_bwwidth,
                    RooAbsReal &_gmean,
                    RooAbsReal &_gsigma,
                    MassType _massType = DeltaMass);

    RooRBWGaussConv(const RooRBWGaussConv &other, const char *name = 0);
    virtual TObject *clone(const char *newname) const { return new RooRBWGaussConv(*this, newname); }
    virtual ~RooRBWGaussConv() {}

    static double D0MASS;
    static double MPI;

  protected:
    RooRealProxy x;
    RooRealProxy bw_mean;
    RooRealProxy bw_width;
    RooRealProxy g_mean;
    RooRealProxy g_sigma;
    MassType massType;

    double evalMass(double val) const;
    double breitwigner(double val) const;
    double gaussian(double val) const;
    double evalbwff(double pval) const;
    double evalphsp(double val) const;

    Double_t evaluate() const;
    Double_t evalP(double mass) const;

    ClassDef(RooRBWGaussConv, 1); // Argus background shape generic PDF (for lower and upper thresholds)
};

#endif
