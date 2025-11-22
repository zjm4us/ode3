///
/// @file vterm.cpp
/// @brief Projectile motion without air resistance, RK4 solver
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

// --- ODEs: without air resistance ---
double f_ri(double x, const vector<double> &y, void *params=0) { return y[1]; }
double f_vi(double x, const vector<double> &y, void *params=0) { (void)x; return 0; } // no air
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
    const char* outfile = "vterm.root";

    // --- parse options ---
    int c;
    while ((c = getopt(argc, argv, "v:t:m:k:s:o:")) != -1) {
        switch(c) {
            case 'v': v0 = atof(optarg); break;
            case 't': theta = atof(optarg); break;
            case 'm': pars.m = atof(optarg); break;
            case 'k': pars.air_k = atof(optarg); break;
            case 's': nsteps = atoi(optarg); break;
            case 'o': outfile = optarg; break;
            case '?':
                fprintf(stderr, "Unknown option `%c'.\n", optopt);
                exit(1);
        }
    }

    void* p_par = (void*)&pars;

    // --- initialize ROOT application ---
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

    double x = 0;
    double xmax = 20;

    auto tgN = RK4SolveN(v_fun, y, nsteps, x, xmax, p_par, f_stop);

    // --- optional canvas ---
    TCanvas* c1 = new TCanvas("c1", "Projectile without air", 800, 600);
    tgN[2].Draw("AL*"); // plot y vs t as example
    c1->Update();

    // --- save ROOT file ---
    TFile tf(outfile, "RECREATE");
    for (unsigned i=0; i<v_fun.size(); i++) tgN[i].Write();
    tf.Close();

    cout << "Final velocity = " << sqrt(y[1]*y[1]+y[3]*y[3]) << endl;
    cout << "Saved ROOT file: " << outfile << endl;
    cout << "Press ^C to exit" << endl;

    theApp.SetIdleTimer(30,".q");
    theApp.Run();
}

