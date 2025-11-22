///
/// @file vterm_air.cpp
/// @brief Projectile motion with air resistance; computes terminal velocity
/// @author You
/// @date 2025
///

#include "RKn.hpp"
#include "TROOT.h"
#include "TApplication.h"
#include "TFile.h"
#include "TCanvas.h"
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

// --- parameters ---
struct Params {
    double g;
    double m;
    double air_k;
};

// --- ODEs: projectile with air resistance ---
double f_ri(double x, const vector<double>& y, void* params=0) { return y[1]; }
double f_vi(double x, const vector<double>& y, void* params=0) { 
    Params* p = (Params*)params;
    double v = sqrt(y[1]*y[1] + y[3]*y[3]);
    return -p->air_k * v * y[1] / p->m;
}
double f_rj(double x, const vector<double>& y, void* params=0) { return y[3]; }
double f_vj(double x, const vector<double>& y, void* params=0) { 
    Params* p = (Params*)params;
    double v = sqrt(y[1]*y[1] + y[3]*y[3]);
    return -p->air_k * v * y[3] / p->m - p->g;
}

// stopping condition: projectile hits ground
double f_stop(double x, const vector<double>& y, void* params=0) {
    (void)x;
    return (y[2] < 0) ? 1 : 0;
}

int main(int argc, char** argv) {
    // default parameters
    Params pars{9.81, 10.0, 0.1};
    double v0 = 100.0;    // initial speed
    double theta = 45.0;  // degrees
    int nsteps = 200;

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
    y[1] = v0 * cos(theta*M_PI/180.0);
    y[2] = 0;
    y[3] = v0 * sin(theta*M_PI/180.0);

    cout << "Simulating projectile fall with air resistance...\n";

    double x = 0.0;
    double xmax = 20.0;

    auto tgN = RK4SolveN(v_fun, y, nsteps, x, xmax, p_par, f_stop);

    // --- compute numerical terminal velocity (average last 10% of speed points) ---
    size_t N = tgN[1].GetN();
    vector<double> vx(N), vy(N);
    for(size_t i=0;i<N;i++){
        vx[i] = tgN[1].GetY()[i];
        vy[i] = tgN[3].GetY()[i];
    }
    size_t start = size_t(0.9*N);
    double vt_num = 0;
    for(size_t i=start;i<N;i++){
        vt_num += sqrt(vx[i]*vx[i] + vy[i]*vy[i]);
    }
    vt_num /= (N-start);

    double vt_analytic = sqrt(pars.m * pars.g / pars.air_k);

    cout << "Estimated terminal velocity (numerical) = " << vt_num << " m/s\n";
    cout << "Analytical terminal velocity = " << vt_analytic << " m/s\n";

    // --- save ROOT file ---
    TFile tf("vterm_air.root","RECREATE");
    for(unsigned i=0;i<v_fun.size();i++) tgN[i].Write();
    tf.Close();

    theApp.SetIdleTimer(30,".q");
    theApp.Run();
}

