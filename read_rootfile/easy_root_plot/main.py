import funcs
from os import remove, chdir, path
from glob import glob


input_path = "/home/astro/blum/read_rootfile/"

chdir(input_path)
for file in glob("*.root"):
    read_array = funcs.read_muon_data(file)

print(read_array)


