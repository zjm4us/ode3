///
/// @file vterm_energy.cpp
/// @brief Projectile motion without air resistance, RK4 solver with energy check
/// @author Malinda
/// @date 2025
///

#include "RKn.hpp"
#include "TROOT.h"
#include "TApplication.h"
#include "TFile.h"
#include "TCanvas.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <getopt.h>
#include <vector>
#include <cmath>

using namespace std;

struct Params {
    double g;
    double m;
    double air_k;
};

// --- ODEs: no air resistance ---
double f_ri(double x, const vector<double> &y, void *params=0) { return y[1]; }
double f_vi(double x, const vector<double> &y, void *params=0) { (void)x; return 0; }
double f_rj(double x, const vector<double> &y, void *params=0) { return y[3]; }
double f_vj(double x, const vector<double> &y, void *params=0) { 
    Params *p = (Params*)params;
    return -p->g; 
}

// stopping condition
double f_stop(double x, const vector<double> &y, void *params=0) {
    (void)x;
    return (y[2] < 0) ? 1 : 0;
}

int main(int argc, char** argv) {
    // --- defaults ---
    double v0 = 100;
    double theta = 45;
    Params pars;
    pars.g = 9.81;
    pars.m = 1.0;
    pars.air_k = 0.0;
    int nsteps = 200;
    const char* outfile = "vterm_energy.root";

    // --- parse command-line options ---
    int c;
    while ((c = getopt(argc, argv, "v:t:m:s:o:")) != -1) {
        switch(c) {
            case 'v': v0 = atof(optarg); break;
            case 't': theta = atof(optarg); break;
            case 'm': pars.m = atof(optarg); break;
            case 's': nsteps = atoi(optarg); break;  // allow changing step size
            case 'o': outfile = optarg; break;
            case '?':
                fprintf(stderr, "Unknown option `%c'.\n", optopt);
                exit(1);
        }
    }

    void* p_par = (void*)&pars;
    TApplication theApp("App", &argc, argv);

    // --- set up RK4 solver ---
    vector<pfunc_t> v_fun(4);
    v_fun[0] = f_ri;
    v_fun[1] = f_vi;
    v_fun[2] = f_rj;
    v_fun[3] = f_vj;

    vector<double> y(4);
    y[0] = 0; 
    y[1] = v0*cos(theta*M_PI/180.0);
    y[2] = 0;
    y[3] = v0*sin(theta*M_PI/180.0);

    cout << "Initial velocity: " << v0 << " m/s" << endl;
    cout << "Launch angle: " << theta << " deg" << endl;
    cout << "Step size (nsteps): " << nsteps << endl;

    double x = 0;
    double xmax = 20;

    auto tgN = RK4SolveN(v_fun, y, nsteps, x, xmax, p_par, f_stop);

    // --- compute energy ---
    double m = pars.m;
    vector<double> E(tgN[0].GetN());
    double E0 = 0;
    for (int i=0; i<tgN[0].GetN(); i++) {
        double vx = tgN[1].GetY()[i];
        double vy = tgN[3].GetY()[i];
        double y_pos = tgN[2].GetY()[i];
        double Etot = 0.5*m*(vx*vx + vy*vy) + m*pars.g*y_pos;
        if (i==0) E0 = Etot;
        E[i] = Etot - E0;
    }

    // --- canvas for energy ---
    TCanvas *cE = new TCanvas("cE", "Energy vs time", 800, 600);
    cE->SetGrid();
    TGraph *tg_E = new TGraph(tgN[0].GetN());
    for (int i=0; i<tgN[0].GetN(); i++)
        tg_E->SetPoint(i, tgN[0].GetX()[i], E[i]);
    tg_E->SetTitle("Energy deviation vs time; t [s]; E-E0 [J]");
    tg_E->Draw("AL*");
    cE->Update();

    // --- save everything ---
    TFile tf(outfile, "RECREATE");
    for (unsigned i=0; i<v_fun.size(); i++) tgN[i].Write();
    tg_E->Write("EnergyDeviation");
    tf.Close();

    cout << "Max energy deviation: ";
    double max_dev = 0;
    for (auto d : E) if (fabs(d) > max_dev) max_dev = fabs(d);
    cout << max_dev << " J" << endl;

    cout << "Saved ROOT file: " << outfile << endl;
    cout << "Press ^C to exit" << endl;
    theApp.SetIdleTimer(30,".q");
    theApp.Run();
}

