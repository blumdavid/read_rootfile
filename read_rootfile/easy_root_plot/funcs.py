import ROOT
import matplotlib.pyplot as plt


class PhotonData:
    def __init__(self):  # this method creates the class object.
        self.index = 0
        self.hit_time = 0
        self.pmt_id = 0
        self.is_cerenkov = 0


class MuonData:
    def __init__(self):  # this method creates the class object.
        self.x_position_init = 0


def read_photon_reco_data(source):
    return_array = []

    rfile = ROOT.TFile(source)
    rtree = rfile.Get("bxtree")

    for event in range(rtree.GetEntries()):
        rtree.GetEntry(event)
        lala = rtree.GetBranch("run").GetLeaf("run").GetValue()
        print(lala)


def read_muon_data(source):
    return_array = []
    rfile = ROOT.TFile(source)
    rtree = rfile.Get("mu")

    for event in rtree:
        for i in range(event.__getattr__("MuMult")):
            muon_data = MuonData()
            muon_data.x_position_init = event.__getattr__("MuInitPosx")[i]

            return_array.append(muon_data)
    return return_array


def hist_from_array(diff_array, set_bins, range_min, range_max, outdir, name, x_label, y_label):

    n, bins, patches = plt.hist(diff_array, bins=set_bins, range=(range_min, range_max))

    '''Nicer font and tex support'''
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.grid(True)

    plt.savefig(outdir + str(name) + ".pdf", bbox_inches='tight')
    plt.savefig(outdir + str(name) + ".png", bbox_inches='tight')

    plt.close()
