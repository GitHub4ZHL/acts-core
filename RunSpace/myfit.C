#include "TF1.h"
#include "TH1F.h"

Double_t myfunction(Double_t *x, Double_t *p)
{
   Double_t xx =x[0];
   //Double_t f = p[0]*std::exp(-0.5*std::pow((xx-p[1])/p[2],2)) + p[3] + p[4]*xx + p[5]*std::pow(xx,2);
   Double_t f = p[0]*std::exp(-0.5*std::pow((xx-p[1])/p[2],2)) + std::exp(-0.5*std::pow((xx-p[3])/p[4],2));
   return f;
}



void myfit()
{
   auto h1 = new TH1F("h1","test",100,-5,5);
   for(int i=0; i<1000; ++i){
      h1->Fill(gRandom->Gaus(0, 1));
   } 


   TF1 *f1 = new TF1("myfunc",myfunction,-5,5,6);
   f1->SetParameter(1,0);
   f1->SetParameter(2,1);
 
   //f1->SetParameters(800,1);
   h1->Fit("myfunc","QS0");
   h1->Draw();
   f1->Draw("same");
}
