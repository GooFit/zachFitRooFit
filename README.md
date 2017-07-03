# RooFit version of ZachFit

There is some extremely basic visualization code to look at the final fit. I used a simple integer number "fitType" to switch between Kpi and K3pi.

First we need to compile and link the classes we wrote ourselves and then can run the macro:
```
root>.L RooArgusGenBG.cxx+
root>.L RooRBWGaussConv.cxx+
root>.x zachFit_roofit.C
```

The "NumCPU" option in the real data fit can be modified to suit your individual system.
