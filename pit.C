#include <iostream>
#include <math.h>
#include <TF1.h>

double psi(double *x, double *par) {
  return pow((2/TMath::Pi()), 0.25)*pow(par[0], -0.5)*exp(-pow(x[0], 2)/pow(par[0], 2));
}

// double psi(double *x, double *par) {
//   return exp(-fabs(x[0])/par[0])/(x[0]/par[0] + par[0]/x[0]);
// }

double scalar_expr(double *x, double *par){
  double A = par[0];
  double U;
  if (fabs(x[0]) < 10){
        U = -0.5;
  }
  else { 
      U = 0;
  }
  TF1 Psi("psi", psi, -10.*A, 10.*A, 1);
  Psi.SetParameter(0, A);
  double result = Psi.Eval(x[0])*(-3.8*Psi.Derivative2(x[0]) + U*Psi.Eval(x[0]));
  return result;
}


double scal_prod(double *x, double *par){
  TF1 scal_exp("sc", scalar_expr, -10.*x[0], 10.*x[0], 1);
  scal_exp.SetParameter(0, x[0]);
  return scal_exp.Integral(-10.*x[0], 10.*x[0]);
}


void pit(double min, double max){  
  //ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("GAUSS");
  gStyle->SetLabelSize(0.05,"xy");
  gStyle->SetTitleSize(0.05,"xy");
  TCanvas *s = new TCanvas("s", "Canvas", 900, 600);
  s->Divide(2, 1);
  TF1 *scalar_pr = new TF1("sp", scal_prod, min, max, 0);
  TF1 *wavefunc = new TF1("wave", psi, -max, max, 1);
  double a = scalar_pr->GetMinimumX(min, max);
  std::cout << "min: " << a << std::endl;
  wavefunc->SetParameter(0, a);

  s->cd(1);

  scalar_pr->SetTitle("<#psi|H|#psi> as a function of a");
  scalar_pr->GetXaxis()->SetTitle("a");
  scalar_pr->GetYaxis()->SetTitle("<#psi|H|#psi>");
  scalar_pr->Draw();

  s->cd(2);

  wavefunc->SetTitle("#psi(x)");
  wavefunc->GetXaxis()->SetTitle("x");
  wavefunc->GetYaxis()->SetTitle("#psi(x)");
  wavefunc->Draw();
}