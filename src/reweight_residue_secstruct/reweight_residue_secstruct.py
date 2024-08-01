import argparse
import os.path
import sys
from os import path
import pandas as pd
import numpy as np
from numba import njit
import multiprocessing as mp
import concurrent.futures
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

'''
Reweighting of GaMD-sampled residue secondary structure propensities.
Secondary structure types are parsed from AmberTools cpptraj secstruct-command output files (secondary structure vs time; DSSP-secondary structure types: ‘ ’, T, G, H, I, S, B, E).
DSSP-secondary structure type definitions:
- loop/irregular ‘ ’ (no letter associated with it; everything that does not match the categories listed below, i.e. none)
- hydrogen-bonded n-turns ‘T’ (3-, 4-, 5-turn)
- n-helices made up of at least two consecutive n-turns (’G’: 3- or 3-10-helix, ’H’: 4- or α-helix, ’I’: 5- or π-helix)
- bend ‘S’ (region defined by geometrical structure, i.e. high curvature with a direction change more than 70˚)
- β-structures ‘B’ and ‘E’:
  *) elementary structures: parallel/anti-parallel β-bridges
  *) compound structures made up of elementary structures: β-ladders → one or more consecutive β-bridges of identical type, β-sheets → one or more β-ladders connected by shared residues
  *) ‘B’: residues that are in isolated/single β-bridge
  *) ‘E’: all remaining residues that participate in β-ladders, i.e. are part of an “extended” β-ladder (continuous stretches of ‘E’ residues are β-strands)
Amber integer naming scheme for DSSP characters: 
‘ ’: 0, T: 6, G: 3, H: 4, I: 5, S: 7, B: 2, E: 1
Associated GaMD boost potentials are parsed from files containing Amber GaMD output (from Amber18), 
which is formatted in "weights.dat" according to input specifications for the corresponding standard GaMD analysis scripts "PyReweighting-*D.py".
Author: Manuel Glaser, manuel.glaser@h-its.org
'''

def parse_command_line():

    parser = argparse.ArgumentParser()

    parser.add_argument("-m", "--mode", type=str, choices=["serial", "parallel"], help="Execution mode (either serial or parallel)", default = "serial")
    
    parser.add_argument("-d", "--dssp", type=str, help="AmberTools cpptraj secstruct-command output file name for secondary structure vs time (DSSP algorithm)", default = "dssp.out")
    
    parser.add_argument("-e", "--energies", type=str, help="File containing GaMD boost potential energies vs time (file format: same as used for the PyReweighting-*D.py scripts)", default = "weights.dat")
   
    parser.add_argument("-T", "--Temperature", type=float, help="Kinetic temperature", default = 300.0)
    
    parser.add_argument("-r", "--residues", type=str, help="List of residues to be reweighted (format: e.g. ASP:1@LYS:2), please make sure you supply existing residues, this is not checked for!")

    parser.add_argument("-o", "--outfiletag", type=str, help="Output file tag", default = "test")
   
    parser.add_argument("-p", "--processes", type=int, help="Number of processes")

    parser.add_argument("-c", "--col_settings", type=str, default = "viridis:0.75,1.0:white,black,black",
                        help="Color settings used for heatmap (':'-separated format, 1st field: heatmap color scheme from https://matplotlib.org/tutorials/colors/colormaps.html," 
                             " 2nd field: threshold values for choosing the text color for heatmap cell annotations (three ranges!), 3nd field: text color)")

    args = parser.parse_args()

    if args.mode == "serial":
        print("Script is executed in serial mode")
    elif args.mode == "parallel":
        print("Script is executed in parallel mode")

    print("Input files/information:")
    if args.dssp:
        print("Dssp output file: %s" % args.dssp)
    if args.energies:
        print("GaMD boost potential energies file: %s" % args.energies)
    if args.col_settings:
        print("Color settings for heatmap: %s" % args.col_settings.split(':'))
    if args.Temperature:
        print("Kinetic temperature: %.2f" % args.Temperature)
    if args.processes:
        try:
            int(args.processes)
        except ValueError:
            print("Input error %s: Please only provide integer-convertable values for number of processes" % i)
            sys.exit(1)
    if args.residues:
        allowed_residues = ["ALA","ASP","CYS","LYS","ARG","GLU","GLN","ASN",
                            "TYR","TRP","GLY","PHE","ILE","SER","THR","VAL",
                            "MET","HID","HIE","HIP","PRO","LEU"]
        for i in args.residues.split('@'):
            residue = i.split(':')[0]
            if not (residue in allowed_residues):
                print("Input error %s: please specify amino acid residues via the three letter code (capital letters)" % residue)
                sys.exit(1)
            try:
                int(i.split(':')[1])
            except ValueError:
                print("Input error %s: Please only provide integer-convertable values as residue numbers" % i)
                sys.exit(1)

    return args

def read_dssp_file(path):

    try:
        dssp_data = pd.read_csv(path,sep='\s+')
        return dssp_data
    except (IOError, TypeError, OSError, ValueError, FileNotFoundError):
        print('Some problems reading the dssp output file %s' % path)
        sys.exit(1)

def read_energies_file(path):

    try:
        energies_data = pd.read_csv(path,sep='\s+',header=None,names=["b_pot_units_kbT", "#Frame", "b_pot_units_kcal/mol"])
        return energies_data
    except (IOError, TypeError, OSError, ValueError, FileNotFoundError):
        print('Some problems reading the boost potential energies file %s, no discrete reweighting is peformed' % path)
        sys.exit(1)

@njit
def compute_b_factor(energy,kin_temp,minus=True):

    beta = 1.0/(0.001987*kin_temp)
    if minus:
        b_factor = np.exp(-beta*energy)
    else:
        b_factor = np.exp(beta*energy)
    return b_factor

def discrete_reweighting(dssp_data,energies_data,kin_temp,res):
  
    print("Performing discrete reweighting of secondary structure propensities for %s" % res) 
 
    rew_data = {"T": {"occs": 0, "av_b_fact": 0.0, "av_ene": 0.0, "std_ene": 0.0, "free_ene": 0.0, "energies": np.array([]), "rew_prop": 0.0, "prop": 0.0, "T":6,"b_factors": np.array([])}, 
                "G": {"occs": 0, "av_b_fact": 0.0, "av_ene": 0.0, "std_ene": 0.0, "free_ene": 0.0, "energies": np.array([]), "rew_prop": 0.0, "prop": 0.0, "G":3,"b_factors": np.array([])},
                "H": {"occs": 0, "av_b_fact": 0.0, "av_ene": 0.0, "std_ene": 0.0, "free_ene": 0.0, "energies": np.array([]), "rew_prop": 0.0, "prop": 0.0, "H":4,"b_factors": np.array([])},
                "I": {"occs": 0, "av_b_fact": 0.0, "av_ene": 0.0, "std_ene": 0.0, "free_ene": 0.0, "energies": np.array([]), "rew_prop": 0.0, "prop": 0.0, "I":5,"b_factors": np.array([])},
                "S": {"occs": 0, "av_b_fact": 0.0, "av_ene": 0.0, "std_ene": 0.0, "free_ene": 0.0, "energies": np.array([]), "rew_prop": 0.0, "prop": 0.0, "S":7,"b_factors": np.array([])},
                "B": {"occs": 0, "av_b_fact": 0.0, "av_ene": 0.0, "std_ene": 0.0, "free_ene": 0.0, "energies": np.array([]), "rew_prop": 0.0, "prop": 0.0, "B":2,"b_factors": np.array([])},
                "E": {"occs": 0, "av_b_fact": 0.0, "av_ene": 0.0, "std_ene": 0.0, "free_ene": 0.0, "energies": np.array([]), "rew_prop": 0.0, "prop": 0.0, "E":1,"b_factors": np.array([])},
                "none": {"occs": 0, "av_b_fact": 0.0, "av_ene": 0.0, "std_ene": 0.0, "free_ene": 0.0, "energies": np.array([]), "rew_prop": 0.0, "prop": 0.0, "none":0,"b_factors": np.array([])}}

    b_factors = np.full(energies_data.shape[0],0.0)
    #print(b_factors.size)

    for i, j in dssp_data.iterrows():
 
        ss_amber_id = dssp_data.iloc[i][res]
        #change this to save compute time, needs only to be done once, not for every single residue again!
        #maybe compute once for a residue and than save to an object and just reload the data for other residues!
        energy = energies_data.iloc[i]["b_pot_units_kcal/mol"]
        b_factor = compute_b_factor(energy,kin_temp,False)
        b_factors[i] = b_factor

        for ss_type in rew_data:

            if ss_amber_id == rew_data[ss_type][ss_type]:
                rew_data[ss_type]["occs"] += 1
                rew_data[ss_type]["b_factors"] = np.append(rew_data[ss_type]["b_factors"],b_factors[i])
                rew_data[ss_type]["energies"] = np.append(rew_data[ss_type]["energies"],energy)
  
    sum_occs = 0      
    for ss_type in rew_data:
    
        if np.sum(rew_data[ss_type]["b_factors"]) != 0.0:
            rew_data[ss_type]["av_b_fact"] = np.mean(rew_data[ss_type]["b_factors"])
            rew_data[ss_type]["av_ene"] = np.mean(rew_data[ss_type]["energies"])
            rew_data[ss_type]["std_ene"] = np.std(rew_data[ss_type]["energies"])

        sum_occs += rew_data[ss_type]["occs"]

    sum_w_props = 0.0
    for ss_type in rew_data:

        sum_w_props += rew_data[ss_type]["av_b_fact"] * (rew_data[ss_type]["occs"] / sum_occs)

    for ss_type in rew_data:

        #maybe just take total sum of av_b_facttors and occs (as sum_av_b_facttors, sum_occs)
        
        rew_data[ss_type]["prop"] = rew_data[ss_type]["occs"] / sum_occs
        w_prop = rew_data[ss_type]["prop"] * rew_data[ss_type]["av_b_fact"]

        rew_data[ss_type]["rew_prop"] = w_prop/sum_w_props

    return rew_data

def sum_up(res_dic,res,rew_key,mode):

    res_sum_prop = 0.0
    res_sum_rew_prop = 0.0
    res_sum_occs = 0
    if mode == "serial":
        for ss_type in res_dic[res][rew_key]:
            res_sum_prop += res_dic[res][rew_key][ss_type]["prop"]
            res_sum_rew_prop += res_dic[res][rew_key][ss_type]["rew_prop"]
            res_sum_occs += res_dic[res][rew_key][ss_type]["occs"]
    else:
        for ss_type in res_dic[res][rew_key].result():
            res_sum_prop += res_dic[res][rew_key].result()[ss_type]["prop"]
            res_sum_rew_prop += res_dic[res][rew_key].result()[ss_type]["rew_prop"]
            res_sum_occs += res_dic[res][rew_key].result()[ss_type]["occs"]

    return (res_sum_prop,res_sum_rew_prop,res_sum_occs)

def write_results(f,res_dic,res,res_sum,prop_key,rew_key,mode,f_tag):

    if mode == "serial":
        f.write(f_tag % (
                res,
                res_dic[res][rew_key]["E"][prop_key],
                res_dic[res][rew_key]["B"][prop_key],
                res_dic[res][rew_key]["G"][prop_key],
                res_dic[res][rew_key]["H"][prop_key],
                res_dic[res][rew_key]["I"][prop_key],
                res_dic[res][rew_key]["T"][prop_key],
                res_dic[res][rew_key]["S"][prop_key],
                res_dic[res][rew_key]["none"][prop_key],
                res_sum))
    else:
        f.write(f_tag % (
                res,
                res_dic[res][rew_key].result()["E"][prop_key],
                res_dic[res][rew_key].result()["B"][prop_key],
                res_dic[res][rew_key].result()["G"][prop_key],
                res_dic[res][rew_key].result()["H"][prop_key],
                res_dic[res][rew_key].result()["I"][prop_key],
                res_dic[res][rew_key].result()["T"][prop_key],
                res_dic[res][rew_key].result()["S"][prop_key],
                res_dic[res][rew_key].result()["none"][prop_key],
                res_sum))

def plot_heatmap(residues,res_dic,rew_key,mode,col_settings):

    #https://matplotlib.org/gallery/images_contours_and_fields/image_annotated_heatmap.html

    ss_types = ["E", "B", "G", "H", "I", "T", "S", "none"]

    props = []
    rew_props = []
    for ss_type in ss_types:
        res_props = []
        res_rew_props = []
        if mode == "serial":
            for res in residues:
                res_props.append(res_dic[res][rew_key][ss_type]["prop"])
                res_rew_props.append(res_dic[res][rew_key][ss_type]["rew_prop"])
        else:
            for res in residues:
                res_props.append(res_dic[res][rew_key].result()[ss_type]["prop"])
                res_rew_props.append(res_dic[res][rew_key].result()[ss_type]["rew_prop"])
        props.append(res_props)
        rew_props.append(res_rew_props)

    props = np.array(props)
    rew_props = np.array(rew_props)

    plt.rc('font', size=7)

    fig, (ax, ax2) = plt.subplots(2,1)

    for axis, data , name in ((ax,props,"GaMD-sampled, biased"),(ax2,rew_props,"GaMD-sampled, unbiased/reweighted")):

        im = axis.imshow(data,cmap=plt.get_cmap(col_settings[0]))

        axis.set_xticks(np.arange(len(residues)))
        axis.set_yticks(np.arange(len(ss_types)))
        axis.set_xticklabels(residues)
        axis.set_yticklabels(ss_types)

        plt.setp(axis.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

        for i in range(len(ss_types)):
            for j in range(len(residues)):

                if 0.0 <= im.norm(data[i, j]) <= float(col_settings[1].split(',')[0]):
                     textcolor = col_settings[2].split(',')[0]
                elif float(col_settings[1].split(',')[0]) < im.norm(data[i, j]) <= float(col_settings[1].split(',')[1]):
                     textcolor = col_settings[2].split(',')[1]
                elif float(col_settings[1].split(',')[1]) < im.norm(data[i, j]) <= 1.0:
                     textcolor = col_settings[2].split(',')[2]

                text = axis.text(j, i, float("{0:.2f}".format(data[i,j])), fontsize=6, ha="center", va="center", color=textcolor)

        axis.set_title(name+" secondary structure propensities")

        divider = make_axes_locatable(axis)
        cax = divider.append_axes('right', size='5%', pad=0.05)

        cbar = axis.figure.colorbar(im, cax=cax, orientation='vertical')
        im.set_clim(0.0, 1.0)

        axis.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
        axis.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
        axis.grid(which="minor", color="w", linestyle='-', linewidth=1)
        axis.tick_params(which="minor", bottom=False, left=False)
        #axis.grid(which="both", color="w", linestyle='-', linewidth=2)

    fig.tight_layout()
    plt.savefig("ss_propensity_heatmap.pdf")

def main():
  
    args = parse_command_line()

    dssp_data = read_dssp_file(args.dssp)
 
    energies_data = read_energies_file(args.energies)

    if args.residues:
        residues = args.residues.split('@')
    else:
        residues = list(dssp_data)
        residues.pop(0)
   
    print("The secondary structure propensity of the following residues is reweighted: ", residues)
    print("Therefore, %i separate reweightings will be performed in %s mode" % (len(residues),args.mode))

    if args.processes:
        processes = args.processes-1
    else:
        processes = mp.cpu_count()-1

    pool = concurrent.futures.ProcessPoolExecutor(max_workers=processes)

    res_dic = {}
    for res in residues:
        
        res_dic.update({res:{"dis_rew": ""}}) 
        if args.mode == "parallel":
            #store "future" object in res_dic dictionary, results can be extracted via future.results() 
            res_dic[res]["dis_rew"] = pool.submit(discrete_reweighting,dssp_data,energies_data,args.Temperature,res)
        elif args.mode == "serial":
            res_dic[res]["dis_rew"] = discrete_reweighting(dssp_data,energies_data,args.Temperature,res)

    #write output
    
    with open(args.outfiletag+"_prop.csv","w+") as f_prop, \
         open(args.outfiletag+"_dis_rew_prop.csv","w+") as f_dis_rew_prop, \
         open(args.outfiletag+"_ene_av.csv","w+") as f_ene_av, \
         open(args.outfiletag+"_ene_std.csv","w+") as f_ene_std, \
         open(args.outfiletag+"_occs.csv","w+") as f_occs:

        f_prop.write("Residue,E,B,G,H,I,T,S,none,sum\n")
        f_dis_rew_prop.write("Residue,E,B,G,H,I,T,S,none,sum\n")
        f_ene_av.write("Residue,E,B,G,H,I,T,S,none,sum\n")
        f_ene_std.write("Residue,E,B,G,H,I,T,S,none,sum\n")
        f_occs.write("Residue,E,B,G,H,I,T,S,none,sum\n")
        for res in residues:
            res_sum_prop, res_sum_dis_rew_prop, res_sum_occs = sum_up(res_dic,res,"dis_rew",args.mode)
            write_results(f_prop,res_dic,res,res_sum_prop,"prop","dis_rew",args.mode,"%s,%f,%f,%f,%f,%f,%f,%f,%f,%f\n")
            write_results(f_dis_rew_prop,res_dic,res,res_sum_dis_rew_prop,"rew_prop","dis_rew",args.mode,"%s,%f,%f,%f,%f,%f,%f,%f,%f,%f\n")
            write_results(f_ene_av,res_dic,res,"","av_ene","dis_rew",args.mode,"%s,%f,%f,%f,%f,%f,%f,%f,%f,%s\n")
            write_results(f_ene_std,res_dic,res,"","std_ene","dis_rew",args.mode,"%s,%f,%f,%f,%f,%f,%f,%f,%f,%s\n")
            write_results(f_occs,res_dic,res,res_sum_occs,"occs","dis_rew",args.mode,"%s,%i,%i,%i,%i,%i,%i,%i,%i,%i\n")
    
    plot_heatmap(residues,res_dic,"dis_rew",args.mode,args.col_settings.split(':'))

if __name__ == "__main__":
    
    main()
