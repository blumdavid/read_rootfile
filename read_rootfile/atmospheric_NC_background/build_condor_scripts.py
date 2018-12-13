""" script to build a number of shell scripts, which can be executed by the condor_desc_file on IHEP cluster:

"""

# path, where shell scripts should be saved (string):
path = "/home/astro/blum/juno/atmoNC/data_NC/NC_scripts_cluster/"

# number of shell scripts, that should be created (integer):
num_scripts = 250


# loop over number of shell scripts:
for index in range(num_scripts):

    # file name of the shell scripts (string):
    filename = path + "script_NC_" + str(index) + ".sh"

    outfile = open(filename, 'w')
    outfile.write("# script to run JUNO detector simulation (tut_detsim.py) for atmospheric NC background:\n")
    outfile.write("\n")
    outfile.write("#!/bin/bash\n")
    outfile.write("export JUNO_OFFLINE_OFF=1\n")
    outfile.write("# source latest offline software (global variables like "
                  "$TUTORIALROOT=/afs/ihep.ac.cn/soft/juno/JUNO-ALL-SLC6/Pre-Release/J18v2r1-branch/offline/"
                  "Examples/Tutorial are set):\n")
    outfile.write("source /afs/ihep.ac.cn/soft/juno/JUNO-ALL-SLC6/Pre-Release/J18v2r1-branch/setup.sh\n")
    outfile.write("# go to folder, where output should be saved:\n")
    outfile.write("cd /junofs/users/dblum/work/atmoNC/detsim_output/\n")
    outfile.write("# run tut_detsim.py:\n")
    outfile.write("python $TUTORIALROOT/share/tut_detsim.py --evtmax 1000 --seed 1 "
                  "--output atmoNC_" + str(index) + ".root --user-output user_atmoNC_" + str(index) + ".root "
                  "hepevt "
                  "--file /junofs/users/dblum/work/atmoNC/folder_NC_onlyC12_250000evts_seed1/"
                  "out_gen_NC_onlyC12_1000evts_seed1_" + str(index) + ".txt --volume pTarget --material LS\n")

    outfile.close()
