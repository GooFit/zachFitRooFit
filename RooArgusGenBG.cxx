/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: RooDKDalitz
 *    File: $Id: RooArgusGenBG.cc,v 1.2 2007/11/08 14:46:13 martinee Exp $
 * Authors:
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *   FMV, Fernando Martinez-Vidal, IFIC-Valencia, martinef@slac.stanford.edu *
 * History:
 *     - Apr 20: creation
 *     - Nov 11, 2007, adapted for D*->D0 bkg (lower threshold)
 *****************************************************************************/

//
// Compile and load the class using:
// root> .L RooArgusGenBG.cxx+
//

#include "RooFit.h"

#include "Riostream.h"
#include <math.h>

#include "RooMath.h"
#include "RooRealConstant.h"
#include "RooRealVar.h"
#include "TMath.h"
#include "TRegexp.h"

#include "RooArgusGenBG.h"

ClassImp(RooArgusGenBG)

    RooArgusGenBG::RooArgusGenBG(const char *name,
                                 const char *title,
                                 RooAbsReal &_m,
                                 RooAbsReal &_m0,
                                 RooAbsReal &_c,
                                 ThresholdType _thresholdType)
    : RooAbsPdf(name, title)
    , m("m", "Mass", this, _m)
    , m0("m0", "Resonance mass", this, _m0)
    , c("c", "Slope parameter", this, _c)
    , c2("c2", "Slope parameter2", this, (RooRealVar &)RooRealConstant::value(0.))
    , p("p", "Power", this, (RooRealVar &)RooRealConstant::value(0.5))
    , thresholdType(_thresholdType) {}

RooArgusGenBG::RooArgusGenBG(const char *name,
                             const char *title,
                             RooAbsReal &_m,
                             RooAbsReal &_m0,
                             RooAbsReal &_c,
                             RooAbsReal &_p,
                             ThresholdType _thresholdType)
    : RooAbsPdf(name, title)
    , m("m", "Mass", this, _m)
    , m0("m0", "Resonance mass", this, _m0)
    , c("c", "Slope parameter", this, _c)
    , c2("c2", "Slope parameter2", this, (RooRealVar &)RooRealConstant::value(0.))
    , p("p", "Power", this, _p)
    , thresholdType(_thresholdType) {}
RooArgusGenBG::RooArgusGenBG(const char *name,
                             const char *title,
                             RooAbsReal &_m,
                             RooAbsReal &_m0,
                             RooAbsReal &_c,
                             RooAbsReal &_c2,
                             RooAbsReal &_p,
                             ThresholdType _thresholdType)
    : RooAbsPdf(name, title)
    , m("m", "Mass", this, _m)
    , m0("m0", "Resonance mass", this, _m0)
    , c("c", "Slope parameter", this, _c)
    , c2("c2", "Slope parameter2", this, _c2)
    , p("p", "Power", this, _p)
    , thresholdType(_thresholdType) {}

RooArgusGenBG::RooArgusGenBG(const RooArgusGenBG &other, const char *name)
    : RooAbsPdf(other, name)
    , m("m", this, other.m)
    , m0("m0", this, other.m0)
    , c("c", this, other.c)
    , c2("c2", this, other.c2)
    , p("p", this, other.p)
    , thresholdType(other.thresholdType) {}
Double_t RooArgusGenBG::evaluate() const {
    if(thresholdType != BrianLowerThreshold) {
        Double_t t = m / m0;
        if((t >= 1 && thresholdType == UpperThreshold) || (t <= 1 && thresholdType == LowerThreshold))
            return 0;
        Double_t u = thresholdType == UpperThreshold ? 1 - t * t : t * t - 1;
        // cout << "c = " << c << " result = " << m*TMath::Power(u,p)*exp(c*u) << endl ;
        assert(u >= 0);
        return m * TMath::Power(u, p) * exp(c * u);
    } else {
        // u^a e^{-b*x+c*x^2)

        Double_t t = m / m0;
        if(t <= 1)
            return 0;
        Double_t u = t * t - 1;
        assert(u >= 0);
        return m * TMath::Power(u, p) * exp(c * u + c2 * u * u);
    }
}
