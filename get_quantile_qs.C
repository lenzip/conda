#include "TRandom3.h"
#include "TH1F.h"
#include "TMath.h"


TRandom3 rand3;

double get_quantile_qs(double integral, double b, double s, bool fc=false){
  TH1F h_qs("h_qs", "h_qs", 10000, 0., 10.);
  //cout << "b " << b << " s " << s << endl;
  for (unsigned int i = 0; i < 10000; ++i){
    //we throw poissonianly around s+b
    int N = rand3.Poisson(s+b);
    double qs=-2*log(TMath::Poisson(N, s+b)/TMath::Poisson(N,N));
    if (fc){
      qs=-2*log(TMath::Poisson(N, s+b)/TMath::Poisson(N, max(N-b, 0.)+b));
      //cout << qs << endl; 
    }
    h_qs.Fill(qs);
  }
  for (unsigned int i = 1; i <= 10000; ++i){
    h_qs.GetXaxis()->SetRangeUser(0, h_qs.GetXaxis()->GetBinCenter(i));
    double this_integral = float(h_qs.Integral())/h_qs.GetEntries();
    //cout << i << " " << h_qs.Integral() << endl;
    //if (abs(this_integral-integral)/integral<0.05){
    if (this_integral>integral){
      //cout << "s: " << s << " " << i << " quantile: " << h_qs.GetXaxis()->GetBinCenter(i) << endl;
      return h_qs.GetXaxis()->GetBinCenter(i);
    }
  }
  return 0;
}
