#include "TF1.h"

// Here is a set of corrections for b jets for the 14 TeV delphes samples. The functions are for
// different ranges of pseudorapidity and parametrized as a function of pT. They have been derived 
// and tested for configuration 4. 


// define the TF1 instances that will be used to find the factor based on the eta value of the bJet

TF1 *f0 = new TF1("f0", "((x<=70)*((0.194143+(0.0137527*x))+(-8.1324e-05*(x*x))))+((x>70)*((0.524402+(0.00391276*x))+(-7.79334e-06*(x*x))))", 20, 200);

TF1 *f1 = new TF1("f1", "((x<=70)*((0.18875+(0.0137748*x))+(-7.97683e-05*(x*x))))+((x>70)*((0.512526+(0.00419686*x))+(-1.05184e-05*(x*x))))", 20, 200);

TF1 *f2 = new TF1("f2", "((x<=70)*((1.06413+(0.0171371*x))+(-0.000174651*(x*x))))+((x>70)*((1.85248+(-0.00881212*x))+(3.30379e-05*(x*x))))", 20, 200);  

TF1 *f3 = new TF1("f3", "((x<=70)*((1.6544+(-0.0360747*x))+(0.000421802*(x*x))))+((x>70)*((7.27617+(-0.0747203*x))+(-1.96959e-05*(x*x))))", 20, 200);

TF1 *f4 = new TF1("f4", "((x<=70)*((0.21795+(0.0131077*x))+(-8.50514e-05*(x*x))))+((x>70)*((0.411553+(0.00570935*x))+(-1.57514e-05*(x*x))))", 20, 200);

TF1 *f5 = new TF1("f5", "((x<=90)*((0.197611+(0.0136277*x))+(-8.01258e-05*(x*x))))+((x>90)*((0.419862+(0.00502803*x))+(-1.18909e-05*(x*x))))", 20, 200);

TF1 *f6 = new TF1("f6", "((x<=80)*((0.296067+(0.0107349*x))+(-5.23645e-05*(x*x))))+((x>80)*((0.496073+(0.00495706*x))+(-1.49215e-05*(x*x))))", 20, 170);

TF1 *f7 = new TF1("f7", "((x<=70)*((0.404015+(0.0107033*x))+(-5.17765e-05*(x*x))))+((x>70)*((0.823724+(0.000758248*x))+(2.37591e-07*(x*x))))", 20, 200);

TF1 *f8 = new TF1("f8", "((x<=70)*((0.362559+(0.0216821*x))+(-0.000176645*(x*x))))+((x>70)*((1.03738+(-0.000227365*x))+(1.29881e-06*(x*x))))", 20, 170);

TF1 *f9 = new TF1("f9", "((x<=70)*((1.06128+(0.00719554*x))+(-6.51501e-05*(x*x))))+((x>70)*((1.54663+(-0.00601546*x))+(2.11055e-05*(x*x))))", 20, 200);

TF1 *f10 = new TF1("f10", "((x<=90)*((1.40908+(0.00121756*x))+(-2.20873e-05*(x*x))))+((x>90)*((2.23614+(-0.0131383*x))+(4.14771e-05*(x*x))))", 20, 200);

TF1 *f11 = new TF1("f11", "((x<=60)*((1.93773+(0.000415024*x))+(-0.000138073*(x*x))))+((x>60)*((2.03251+(-0.00906801*x))+(2.45387e-05*(x*x))))", 20, 200);

TF1 *f12 = new TF1("f12", "1", 20, 200);

// returns the appropriate TF1 from above
TF1 * f(Double_t eta){

    if (abs(eta) <=  0.3) return f0;

    if (0.3  <= abs(eta) <=  0.6) return f1;
    if (3.0  <= abs(eta) <=  4.0) return f2;
    //if 4.0  <= abs(eta) return f3
    if (0.6  <= abs(eta) <=  0.9) return f4;
    if (0.9  <= abs(eta) <=  1.2) return f5;
    if (1.2  <= abs(eta) <=  1.5) return f6;
    if (1.5  <= abs(eta) <=  1.8) return f7; 
    if (1.8  <= abs(eta) <=  2.1) return f8;
    if (2.1  <= abs(eta) <=  2.4) return f9;
    if (2.4  <= abs(eta) <=  2.7) return f10;
    if (2.7  <= abs(eta) <=  3.0) return f11;
    else return f12;

}

// retrns a factor such that correctedBJetPT = factor(pt,eta) * pt
double factor(Double_t pt, Double_t eta){

    TF1 func = (*f(eta));
    //if (not func) return 1.0;
    
    Double_t xMax = func.GetXmax();

    Double_t xMin = func.GetXmin();
    Double_t x = pt;

    if (xMax < x) x = xMax;

    if (x < xMin) x = xMin;

    Double_t v = func.Eval(x);

    if (v != 0) return 1.0/v;
    else if (v == 0) return 0;

}
