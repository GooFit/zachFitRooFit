/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: RooDKDalitz
 *    File: $Id: RooArgusGenBG.rdl,v 1.2 2007/11/08 14:46:17 martinee Exp $
 *  See .cc file
 *****************************************************************************/

#ifndef ROOARGUSGENBG
#define ROOARGUSGENBG

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"

//class RooRealVar;
//class RooAbsReal;

class RooArgusGenBG : public RooAbsPdf {
public:

  enum ThresholdType { UpperThreshold, LowerThreshold, BrianLowerThreshold };

  RooArgusGenBG() {};
  RooArgusGenBG(const char *name, const char *title,
		RooAbsReal& _m, RooAbsReal& _m0, RooAbsReal& _c,
		ThresholdType _thresholdType=UpperThreshold);

  RooArgusGenBG(const char *name, const char *title,
		RooAbsReal& _m, RooAbsReal& _m0, RooAbsReal& _c,
		RooAbsReal& _p, ThresholdType _thresholdType=UpperThreshold);

  RooArgusGenBG(const char *name, const char *title,
		RooAbsReal& _m, RooAbsReal& _m0, RooAbsReal& _c,
		RooAbsReal& _c2, RooAbsReal& _p, ThresholdType _thresholdType=BrianLowerThreshold);

  RooArgusGenBG(const RooArgusGenBG& other,const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooArgusGenBG(*this,newname); }
  virtual ~RooArgusGenBG() = default;

  //Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  //Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

protected:
  RooRealProxy m ;
  RooRealProxy m0 ;
  RooRealProxy c ;
  RooRealProxy c2 ;
  RooRealProxy p ;
  ThresholdType thresholdType;

  Double_t evaluate() const ;
  
  // Argus background shape generic PDF (for lower and upper thresholds)
  ClassDef(RooArgusGenBG,1);
};

#endif
