import ROOT
import numpy as np
from matplotlib import pyplot as plt
import NC_background_functions

filename = "/home/astro/blum/juno/atmoNC/data_Julia/gntp.101.gst.root"

rfile = ROOT.TFile(filename)
rtree = rfile.Get("gst")

num_events = rtree.GetEntries()
# num_events = 100

num = 0
num_i = 0
num_f = 0

for index in range(10000):
    rtree.GetEntry(index)

    nc = int(rtree.GetBranch('nc').GetLeaf('nc').GetValue())

    tgt = int(rtree.GetBranch('tgt').GetLeaf('tgt').GetValue())

    if nc == 1 and tgt == 1000060120:

        # print("------------ event {0:d} ---------".format(index))

        pdgi = []
        pdgf = []
        print_i = False

        num += 1

        ni = int(rtree.GetBranch('ni').GetLeaf('ni').GetValue())

        for index1 in range(ni):

            if int(rtree.GetBranch('pdgi').GetLeaf('pdgi').GetValue(index1)) == 1000060120:
                print_i = True

            pdgi.append(int(rtree.GetBranch('pdgi').GetLeaf('pdgi').GetValue(index1)))

        if print_i:
            print("pdgi = {0}".format(pdgi))

        nf = int(rtree.GetBranch('nf').GetLeaf('nf').GetValue())

        for index2 in range(nf):

            pdgf.append(int(rtree.GetBranch('pdgf').GetLeaf('pdgf').GetValue(index2)))

        if print_i:
            print("pdgf = {0}".format(pdgf))



