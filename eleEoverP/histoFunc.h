#ifndef histoFunc_h
#define histoFunc_h

#include "TH1.h"


// copied from https://raw.githubusercontent.com/lbrianza/ECALELF/newMaster/EOverPCalibration/CommonTools/histoFunc.h
// this class provides a histogram template to fit E/P distribution
// see also https://github.com/lbrianza/ECALELF/blob/newMaster/ZFitter/bin/ZFitter.cpp#L1309-L2700


class histoFunc
{
 public:
  
  
  //! ctor
  histoFunc(TH1F* histo)
  {
    histo_p = histo;
  };
  
  
  //! dtor
  ~histoFunc()
  {};
  
  //norm histo
  double GetIntegral(){
    double nn = histo_p -> Integral();
    return(nn);
  }

  //! operator()
  double operator()(double* x, double* par)
  {
    double xx = par[1] * (x[0] - par[2]);
    
    double xMin = histo_p -> GetBinCenter(1);
    double xMax = histo_p -> GetBinCenter(histo_p -> GetNbinsX());
    
    
    
    if( (xx < xMin) || (xx >= xMax) )
      return 1.e-10;
    
    else
    {  
      int bin = histo_p -> FindBin(xx);
      int bin1 = 0;
      int bin2 = 0;
      
      if(xx >= histo_p -> GetBinCenter(bin))
      {
        bin1 = bin;
        bin2 = bin+1;
      }
      
      else
      {
        bin1 = bin-1;
        bin2 = bin;
      }
      
      
      double x1 = histo_p -> GetBinCenter(bin1);
      double y1 = histo_p -> GetBinContent(bin1);
      
      double x2 = histo_p -> GetBinCenter(bin2);
      double y2 = histo_p -> GetBinContent(bin2);
      
      double m = 1. * (y2 - y1) / (x2 - x1);
      
      
      
      if( (y1 + m * (xx - x1)) < 1.e-10)
        return 1.e-10;
      
      
      return par[0] * par[1] * (y1 + m * (xx - x1));
    }
    
    return 1.e-10;  
  }
  
  
  
 private:
  
  TH1F* histo_p;
};

#endif
