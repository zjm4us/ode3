///
/// @file vterm_vs_mass.cpp
/// @brief Analytical terminal velocity vs mass
/// @author MALINDA
/// @date 20 Nov 2025
///

#include <TCanvas.h>
#include <TGraph.h>
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

int main() {
    double g = 9.81;       // gravity [m/s^2]
    double air_k = 0.1;    // air drag coefficient

    int npoints = 20;
    vector<double> masses(npoints);
    vector<double> vt_analytic(npoints);

    for(int i=0; i<npoints; i++){
        // log-spaced masses: 0.001 kg -> 10 kg
        masses[i] = 0.001 * pow(10, i * log10(10000.0)/(npoints-1));
        vt_analytic[i] = sqrt(masses[i]*g / air_k);

        cout << "Mass = " << masses[i] << " kg, Analytical vt = " 
             << vt_analytic[i] << " m/s" << endl;
    }

    // Create ROOT graph
    TGraph *tg_a = new TGraph(npoints, &masses[0], &vt_analytic[0]);
    tg_a->SetLineColor(kBlue);
    tg_a->SetLineWidth(2);
    tg_a->SetTitle("Analytical Terminal Velocity vs Mass;Mass [kg];v_t [m/s]");

    TCanvas *c1 = new TCanvas("c1","Analytical Terminal Velocity",800,600);
    c1->SetLogx();
    tg_a->Draw("AL"); // A = axis, L = line
    c1->Update();
    c1->Print("vt_vs_mass.png");

    cout << "Plot saved as vt_vs_mass.png" << endl;
}

