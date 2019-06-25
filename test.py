import ROOT
import numpy as np
from matplotlib import pyplot as plt

filename = "/local/scratch1/pipc51/astro/blum/detsim_output_data/user_atmoNC_0.root"

rfile = ROOT.TFile(filename)

rtree_evt = rfile.Get("evt")

rtree_evt.GetEntry(0)

n_photons = int(rtree_evt.GetBranch('nPhotons').GetLeaf('nPhotons').GetValue())

hittime = []

for index in range(n_photons):
    pmt_id = int(rtree_evt.GetBranch('pmtID').GetLeaf('pmtID').GetValue(index))

    if pmt_id == 16085:
        hit_time = float(rtree_evt.GetBranch('hitTime').GetLeaf('hitTime').GetValue(index))
        hittime.append(hit_time)

plt.hist(hittime)
plt.show()

