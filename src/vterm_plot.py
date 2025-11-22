#!/usr/bin/env python3
"""
Plot projectile motion results from vterm.root
- Trajectory: y vs x
- Velocity magnitude vs time
- Energy vs time
"""

import ROOT as r
import sys
from math import sqrt

# -------------------------------
# Open ROOT file
# -------------------------------
tf = r.TFile("vterm.root")
if tf.IsZombie():
    print("vterm.root not found! Run vterm.cpp first.")
    sys.exit(1)

# Load TGraphs (assumes 0:x, 1:vx, 2:y, 3:vy)
tg_x = tf.Get("xy0")
tg_vx = tf.Get("xy1")
tg_y = tf.Get("xy2")
tg_vy = tf.Get("xy3")

n = tg_x.GetN()

# Extract arrays
time = tg_x.GetX()
xval = tg_x.GetY()
vx = tg_vx.GetY()
yval = tg_y.GetY()
vy = tg_vy.GetY()

# -------------------------------
# Prepare graphs
# -------------------------------
tg_traj = r.TGraph(n, xval, yval)
tg_traj.SetTitle("Projectile Trajectory; x [m]; y [m]")

tg_v = r.TGraph()
tg_v.SetTitle("Velocity magnitude vs time; t [s]; v [m/s]")

tg_E = r.TGraph()
tg_E.SetTitle("Energy change vs time; t [s]; E [J/kg]")

# Mass assumed 1 kg for simplicity
m = 1.0
E0 = 0
for i in range(n):
    vmag = sqrt(vx[i]**2 + vy[i]**2)
    Ekin = 0.5 * m * vmag**2
    Epot = m * 9.81 * yval[i]
    Etot = Ekin + Epot
    if i == 0: E0 = Etot
    tg_v.SetPoint(i, time[i], vmag)
    tg_E.SetPoint(i, time[i], Etot - E0)

# -------------------------------
# Draw canvas
# -------------------------------
c = r.TCanvas("c", "Projectile Motion", 1000, 800)
c.Divide(2,2)

c.cd(1)
tg_traj.Draw("AL*")

c.cd(2)
tg_v.Draw("AL*")

c.cd(3)
tg_E.Draw("AL*")

c.cd(4)
r.gPad.SetGrid()
r.gPad.Update()

c.Update()
c.Print("vterm_plots.png")

print("Plots saved as vterm_plots.png")

