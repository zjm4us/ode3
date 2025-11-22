#!/usr/bin/env python3
"""
Compare energy conservation for different step sizes.
Requires vterm_energy ROOT files saved for multiple step sizes.
"""

import ROOT as r
from math import sqrt

# Files for different step sizes
files = {
    50: "vterm_s50.root",
    200: "vterm_s200.root",
    1000: "vterm_s1000.root"
}

colors = {50: r.kRed, 200: r.kBlue, 1000: r.kGreen+2}

# Create canvas
c = r.TCanvas("c", "Energy Deviation vs Time", 800, 600)
c.SetGrid()
mg = r.TMultiGraph()
mg.SetTitle("Energy deviation vs time; t [s]; E - E0 [J]")

# Loop over files and create graphs
for nsteps, fname in files.items():
    tf = r.TFile(fname)
    if tf.IsZombie():
        print(f"File {fname} not found! Run vterm_energy with -s {nsteps}")
        continue

    tg_E = tf.Get("EnergyDeviation")
    if not tg_E:
        # If your ROOT file uses another name, get TGraph for total energy
        tg_y = tf.Get("xy2")
        tg_vy = tf.Get("xy3")
        tg_vx = tf.Get("xy1")
        tg_x = tf.Get("xy0")
        n = tg_x.GetN()
        tg_E = r.TGraph()
        m = 1.0
        E0 = 0
        for i in range(n):
            t = tg_x.GetX()[i]
            y = tg_y.GetY()[i]
            vx = tg_vx.GetY()[i]
            vy = tg_vy.GetY()[i]
            Etot = 0.5*m*(vx**2 + vy**2) + m*9.81*y
            if i == 0:
                E0 = Etot
            tg_E.SetPoint(i, t, Etot - E0)

    tg_E.SetLineColor(colors[nsteps])
    tg_E.SetMarkerColor(colors[nsteps])
    tg_E.SetMarkerStyle(20)
    tg_E.SetTitle(f"{nsteps} steps")
    mg.Add(tg_E)

# Draw multigraph
mg.Draw("ALP")
c.BuildLegend()
c.Update()
c.Print("energy_vs_steps.png")
print("Energy deviation plot saved as energy_vs_steps.png")

