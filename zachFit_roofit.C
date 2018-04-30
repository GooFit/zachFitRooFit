//////////////////////////////////////////////////////////////////////////
// ROOT version:
// root>.L RooArgusGenBG.cxx+         // Compile and load created class
// root>.L RooRBWGaussConv.cxx+       // Compile and load created class
// root>.x zachFit_roofit.C           // run compiled code
// Or compile with cmake.
/////////////////////////////////////////////////////////////////////////

#ifndef __CLING__
#include "RooGlobalFunc.h"
#endif
#include "RooAddPdf.h"
#include "RooClassFactory.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooMinuit.h"
#include "RooNLLVar.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TStyle.h"

#include "RooArgusGenBG.h"
#include "RooRBWGaussConv.h"

#include <CLI11.hpp>
#include <rang.hpp>
#include <thread>

using namespace RooFit;

// 0 for Kpi, 1 for K3pi
void zachFit_roofit_compute(int fitType = 0, int ncpu = 0) {
    TString MC_TXTNAME, RD_TXTNAME, plot_name;

    if(ncpu == 0)
        ncpu = std::thread::hardware_concurrency();
    if(ncpu == 0)
        ncpu = 1;

    std::cout << rang::fg::magenta << "Number of threads: " << ncpu << rang::style::reset << std::endl;

    if(fitType == 0) {
        MC_TXTNAME = "DstarWidthAnalysis_Data/DstarWidth_D0ToKpi_deltaM_MC.txt";
        RD_TXTNAME = "DstarWidthAnalysis_Data/DstarWidth_D0ToKpi_deltaM_Data.txt";
        plot_name  = "D0ToKpi";
    } else if(fitType == 1) {
        MC_TXTNAME = "DstarWidthAnalysis_Data/DstarWidth_D0ToK3pi_deltaM_MC.txt";
        RD_TXTNAME = "DstarWidthAnalysis_Data/DstarWidth_D0ToK3pi_deltaM_Data.txt";
        plot_name  = "D0ToK3pi";
    }

    RooRealVar *dM = new RooRealVar("dm", "", 0.1395, 0.1665);
    dM->setBins(540); // 50keV bins

    RooDataSet *MC_filedata = RooDataSet::read(MC_TXTNAME, *dM);
    RooAbsData *MC_data     = MC_filedata->reduce(RooFit::CutRange("fitrange"));
    RooDataHist MC_datahist("MC_datahist", "", RooArgSet(*dM), *MC_data);

    RooRealVar MC_mean1("MC_mean1", "", 0.14542, 0.143, 0.148);
    RooRealVar MC_mean2("MC_mean2", "", 0.145415, 0.145, 0.1465);
    RooRealVar MC_mean3("MC_mean3", "", 0.14542, 0.144, 0.147);

    RooRealVar MC_sigma1("MC_sigma1", "", 0.0001, 0.00005, 0.002);
    RooRealVar MC_sigma2("MC_sigma2", "", 0.0005, 0.00001, 0.005);
    RooRealVar MC_sigma3("MC_sigma3", "", 0.0002, 0.000005, 0.001);

    RooGaussian MC_gauss1("MC_gauss1", "", *dM, MC_mean1, MC_sigma1);
    RooGaussian MC_gauss2("MC_gauss2", "", *dM, MC_mean2, MC_sigma2);
    RooGaussian MC_gauss3("MC_gauss3", "", *dM, MC_mean3, MC_sigma3);

    RooRealVar pionmass("pionmass", "", 0.13957);
    RooRealVar MC_aslope("MC_aslope", "", -20.0, -100.0, 10.0);
    RooRealVar MC_apower("MC_apower", "", 1.8, 0.1, 10.0); // 2.5

    RooArgusGenBG MC_argus("MC_argus", "", *dM, pionmass, MC_aslope, MC_apower, RooArgusGenBG::LowerThreshold);

    RooRealVar MC_gfrac1("MC_gfrac1", "", 0.75, 0.0, 0.9);
    RooRealVar MC_gfrac2("MC_gfrac2", "", 0.02, 0.0, 0.3);
    RooRealVar MC_afrac("MC_afrac", "", 0.006, 0.0, 0.1);

    RooAddPdf MC_signal("MC_signal", "", RooArgList(MC_gauss1, MC_gauss2, MC_gauss3), RooArgList(MC_gfrac1, MC_gfrac2));
    RooAddPdf MC_total("MC_total", "", RooArgList(MC_argus, MC_signal), RooArgList(MC_afrac));

    RooNLLVar *MC_nll    = new RooNLLVar("MC_nll", "", MC_total, MC_datahist);
    RooMinuit *MC_minuit = new RooMinuit(*MC_nll);

    MC_minuit->setStrategy(2);
    MC_minuit->hesse();
    MC_minuit->migrad();
    MC_minuit->hesse();

    // fit to the data
    MC_mean1.setConstant(true);
    MC_mean2.setConstant(true);
    MC_mean3.setConstant(true);
    MC_sigma1.setConstant(true);
    MC_sigma2.setConstant(true);
    MC_sigma3.setConstant(true);
    MC_gfrac1.setConstant(true);
    MC_gfrac2.setConstant(true);
    MC_afrac.setConstant(true);
    MC_aslope.setConstant(true);
    MC_apower.setConstant(true);

    RooDataSet *rdfiledata = RooDataSet::read(RD_TXTNAME, *dM);
    RooAbsData *rddata     = rdfiledata->reduce(RooFit::CutRange("fitrange"));
    RooDataHist rddhist("rddhist", "", RooArgSet(*dM), *rddata);

    RooRealVar delta("RD_delta", "", 0.0000001, -0.00005, 0.0001);
    RooRealVar epsilon("RD_epsilon", "", 0.05, 0.0, 0.4);
    RooRealVar width_bw("RD_width_bw", "", 0.000085, 0.000060, 0.0001);
    RooRealVar dummyzero("RD_dummyzero", "", 0);
    RooFormulaVar rdmean1("RD_mean1", "", "@0+@1", RooArgList(MC_mean1, delta));
    RooFormulaVar rdmean2("RD_mean2", "", "@0+@1", RooArgList(MC_mean2, delta));
    RooFormulaVar rdmean3("RD_mean3", "", "@0+@1", RooArgList(MC_mean3, delta));

    // sigma^{rd}_i = sigma^{mc}_i * (1+epsilon) + deltaSigma_i = mc parameter *(1 + scale factor) + extra shift used in
    // chi^2
    RooFormulaVar rdsigma1("RD_sigma1", "", "(1+@0)*@1", RooArgList(epsilon, MC_sigma1));
    RooFormulaVar rdsigma2("RD_sigma2", "", "(1+@0)*@1", RooArgList(epsilon, MC_sigma2));
    RooFormulaVar rdsigma3("RD_sigma3", "", "(1+@0)*@1", RooArgList(epsilon, MC_sigma3));

    RooRBWGaussConv convo_1("RD_convo_1", "", *dM, rdmean1, width_bw, dummyzero, rdsigma1);
    RooRBWGaussConv convo_2("RD_convo_2", "", *dM, rdmean2, width_bw, dummyzero, rdsigma2);
    RooRBWGaussConv convo_3("RD_convo_3", "", *dM, rdmean3, width_bw, dummyzero, rdsigma3);

    RooRealVar slope("RD_slope", "", -3.0, -10.0, 10.0);
    RooArgusGenBG mcargus("mcargus", "", *dM, pionmass, MC_aslope, MC_apower, RooArgusGenBG::LowerThreshold);

    RooRealVar rpower("rpower", "", 0.5);
    RooArgusGenBG rdargus("rdargus", "", *dM, pionmass, slope, rpower, RooArgusGenBG::LowerThreshold);
    RooRealVar bkg_frac("RD_bkg_frac", "", 0.04, 0.0, 0.4);
    RooRealVar nsig("RD_nsig", "", 138500, 100000, 0.99 * rddhist.sumEntries());
    RooRealVar nbkg("RD_nbkg", "", 3000, 1000, 50000);
    RooRealVar nmcarg("RD_nmcarg", "", (MC_afrac.getVal() * nsig.getVal()));
    RooRealVar nmcg1("RD_nmcg1", "", (MC_gfrac1.getVal() * nsig.getVal()));
    RooRealVar nmcg2("RD_nmcg2", "", (MC_gfrac2.getVal() * nsig.getVal()));
    RooRealVar nmcg3(
        "RD_nmcg3", "", (1.0 - MC_afrac.getVal() - MC_gfrac1.getVal() - MC_gfrac2.getVal()) * nsig.getVal());

    RooAddPdf signal(
        "signal", "", RooArgList(mcargus, convo_1, convo_2, convo_3), RooArgList(nmcarg, nmcg1, nmcg2, nmcg3));
    RooAddPdf RD_total("RD_total", "", RooArgList(rdargus, signal), RooArgList(nbkg, nsig));
    RooNLLVar *rdnll    = new RooNLLVar("rdnll", "", RD_total, rddhist, RooFit::NumCPU(ncpu), RooFit::Extended());
    RooMinuit *rdMinuit = new RooMinuit(*rdnll);
    TStopwatch _timer2;
    std::cout << rang::fg::magenta << "Starting fit timer" << rang::fg::green << std::endl;
    _timer2.Start();

    rdMinuit->setStrategy(1);
    rdMinuit->hesse();
    rdMinuit->migrad();
    rdMinuit->hesse();

    std::cout << rang::fg::magenta << "Time at the end of fit = " << _timer2.RealTime() << " (real) "
              << _timer2.CpuTime() << " (cpu) seconds" << rang::fg::magenta << std::endl;

    TStyle *_gStyle = new TStyle();

    _gStyle->SetCanvasBorderMode(0);
    _gStyle->SetCanvasColor(10);
    _gStyle->SetFrameFillColor(10);
    _gStyle->SetFrameBorderMode(0);
    _gStyle->SetPadColor(0);
    _gStyle->SetStatColor(0);
    _gStyle->SetFillColor(0);
    _gStyle->SetFuncWidth(1);
    _gStyle->SetLineWidth(1);
    _gStyle->SetLineColor(1);
    _gStyle->SetPalette(1, 0);
    _gStyle->SetPadRightMargin(0.15);
    _gStyle->SetOptStat();

    TCanvas *canvas = new TCanvas("c1");

    RooPlot *mframe = dM->frame();
    rddhist.plotOn(mframe);
    RD_total.plotOn(mframe);
    mframe->Draw();
    mframe->GetXaxis()->SetTitle("#Delta m [GeV]");
    mframe->GetXaxis()->SetLabelSize(0.03);
    mframe->GetYaxis()->SetLabelSize(0.03);

    canvas->SetLogy(true);
    mframe->SetMinimum(0.1);
    mframe->GetYaxis()->SetTitle("Events / 50 keV");
    mframe->GetXaxis()->SetTitle("#Delta m [GeV]");
    mframe->GetXaxis()->SetLabelSize(0.03);
    mframe->GetYaxis()->SetLabelSize(0.03);
    canvas->SaveAs("zachsFit_" + plot_name + "_roofit_semilog.pdf");
}

/// Called from ROOT
void zachFit_roofit() {
    gROOT->ProcessLineSync(".L RooArgusGenBG.cxx+");
    gROOT->ProcessLineSync(".L RooRBWGaussConv.cxx+");
    zachFit_roofit_compute();
}

#ifndef __CLING__
#include <csignal>

void signal_handler(int s) {
    std::cout << std::endl << rang::style::reset << rang::fg::red << rang::style::bold;
    std::cout << "zachFit: Control-C detected, exiting..." << rang::style::reset << std::endl;
    std::exit(1); // will call the correct exit func, no unwinding of the stack though
}

/// Called from command line
int main(int argc, char **argv) {
    // Nice exit
    std::atexit([]() { std::cout << rang::style::reset; });

    // Nice Control-C
    struct sigaction sigIntHandler;
    sigIntHandler.sa_handler = signal_handler;
    sigemptyset(&sigIntHandler.sa_mask);
    sigIntHandler.sa_flags = 0;
    sigaction(SIGINT, &sigIntHandler, nullptr);

    CLI::App app{"zachFit"};

    int fit_type = 0, ncpu = 0;
    app.add_option("-f,--fit-type,fit-type", fit_type, "0 for Kpi, 1 for K3pi", true);
    app.add_option("-n,--numcpu,ncpu", ncpu, "Number of CPUs to use", true);

    try {
        app.parse(argc, argv);
    } catch(const CLI::ParseError &e) {
        std::cout << (e.get_exit_code() == 0 ? rang::fg::blue : rang::fg::red);
        return app.exit(e);
    }

    zachFit_roofit_compute(fit_type, ncpu);
}
#endif
