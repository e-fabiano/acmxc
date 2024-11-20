#!/usr/bin/env python3


import acmxclib
import argparse
import tools
from numpy import array as nparray
from numpy import sum as npsum
import acmlib


###################################################################
#   FUNCTIONS
###################################################################

# manage the command line options
def manage_options():
    parser = argparse.ArgumentParser(description='PREREQUISITE: Set up a valid input file to run an ACM calculation with TURBOMOLE or CRYSTAL, then run the script.')
#
    parser.add_argument("-p","--prog","--program",
                        required=True,
                        choices=["turbomole","crystal"],
                        metavar="<string>",
                        help="This is a mandatory option to select the program to be used (turbomole, crystal)")
#
    parser.add_argument("-f","--formula",
                        default="isi",
                        choices=["isi","revisi","spl","lb","spl2","mpacf1","genisi","dpi","hfac24","mp2","hflda"],
                        help="Formula to be used. Options: isi, revisi, genisi, spl, lb, spl2, mpacf1, dpi, hfac24, mp2, hflda",
                        metavar="<string>")
#
    parser.add_argument("-w","-wfunc",
                        default=None,
                        choices=["pc","hpc","mpc","hfpc","lda"],
                        help="W_inf functional. Options: pc, hpc, mpc, hfpc, lda",
                        metavar="<string>")
#
    parser.add_argument("-n","-nthreads",
                        default=1,
                        type=int,
                        help="Number of threads to use (default: 1)",
                        metavar="<int>")
#
    parser.add_argument("-d","-dir",
                        default=None,
                        help="Base path of <program> (default=set from environment variable)",
                        metavar="<string>")
#
    parser.add_argument("-i","--input",
                        default="input",
                        help="root of the crystal file name (not including the .d12 .d3 and .d4 extensions)",
                        metavar="<string>")

#
    parser.add_argument("--metal",
                        action='store_true',
                        default=False,
                        help="Sets E_mp2=-inf. Use for calculations of systems with vanishing gap (default=False)")
#
    parser.add_argument("--int","--interaction",
                        default=None,
                        metavar="<string>",
                        help="Full path of a file with information for interaction eenrgy calculations")
#
    parser.add_argument("--rerun",  
                        action='store_true',             
                        default=False, 
                        help="Rerun calculation (default=False)")
#
    options = vars(parser.parse_args())
#
#
    return options



def manage_iteraction_input(intfile):
    if (intfile == None):
        intlist = ["."]
        intcoeff = [1]
    else:
        with open(intfile,"r") as f:
            intfile_content = f.readlines()
            intlist = [line.split()[1] for line in intfile_content]
            intcoeff = [float(line.split()[0]) for line in intfile_content]
    return intlist, intcoeff
                




####################################################################
#  MAIN CODE
####################################################################

# initial options
options = manage_options()
program = options["prog"]
tdir = options["d"]
prog_input = options["input"]
ncpu = options["n"]
formula = options["formula"]
wfunc = options["w"]
rerun = options["rerun"]
metal_mode = options["metal"]
intfile = options["int"]


# information on systems to be treated
intlist, intcoeff = manage_iteraction_input(intfile)
intcoeff = nparray(intcoeff)


# main loop
acmxc_simul = []
for path in intlist:
    acmxc_simul.append(acmxclib.acmxc(path=path,program=program,tdir=tdir,prog_input=prog_input,ncpu=ncpu,formula=formula,wfunc=wfunc,rerun=rerun,metal_mode=metal_mode,verbose=False))

    
tools.print_header()
print("")

    
i = 1    
for isimul in acmxc_simul:
    if (len(acmxc_simul)>1): print(f"==== System {i} ====")
    isimul.verbose = True
    tools.print_options([isimul.program_name,isimul.tdir,isimul.baseinput,isimul.ncpu,isimul.acm_formula,isimul.wfunc,isimul.rerun,isimul.metal_mode,isimul.w34])
    isimul.run_program()
    isimul.extract_results()
    isimul.compute_acm_xc_energy()
    isimul.print_results()
    i = i + 1
    print("")



# INTERACTION ENERGY ANALYSIS    
if (len(acmxc_simul)>1):
    print("===== INTERACTION ANALYSIS ====")
    print("    (all energies are in Ha)")
    scfene = nparray([isimul.scfene for isimul in acmxc_simul])
    int_scfene = npsum(intcoeff*scfene)
    xene = nparray([isimul.xene for isimul in acmxc_simul])
    int_xene = npsum(intcoeff*xene)
    wene = nparray([isimul.wene for isimul in acmxc_simul])
    int_wene = npsum(intcoeff*wene)
    w1ene = nparray([isimul.w1ene for isimul in acmxc_simul])
    int_w1ene = npsum(intcoeff*w1ene)
    mp2ene = nparray([isimul.mp2ene for isimul in acmxc_simul])
    int_mp2ene = npsum(intcoeff*mp2ene)
    correne = nparray([isimul.correne for isimul in acmxc_simul])
    int_correne = npsum(intcoeff*correne)
    int_xcene = int_xene + int_correne
    print(" %-23s %16.10f" %("Int. SCF energy:",int_scfene))
    print(" %-23s %16.10f" %("Int. exchange energy:",int_xene))
    print(" %-23s %16.10f" %("Int. W_inf energy:",int_wene))
    print(" %-23s %16.10f" %("Int. W1_inf energy:",int_w1ene))
    print(" %-23s %16.10f" %("Int. MP2 energy:",int_mp2ene))
    print(" %-23s %16.10f" %("Int. ACM corr. energy:",int_correne))
    print(" %-23s %16.10f" %("Int. XC energy:",int_xcene))
#   SCC
    inf_scfene = -(int_scfene - intcoeff[0]*scfene[0])
    inf_xene = -(int_xene - intcoeff[0]*xene[0])
    inf_wene = -(int_wene - intcoeff[0]*wene[0])
    inf_w1ene = -(int_w1ene - intcoeff[0]*w1ene[0])
    inf_mp2ene = -(int_mp2ene - intcoeff[0]*mp2ene[0])
    inf_correne = acmlib.compute_acm(formula,inf_xene,inf_wene,inf_w1ene,inf_mp2ene)
    scc_correne = correne[0]*intcoeff[0] - inf_correne 
    scc_xcene = int_xene + scc_correne
    scc = scc_correne - int_correne
#   FINAL OUTPUT
    int_totene = int_scfene + int_correne
    scc_totene = int_scfene + scc_correne
    print("")
    print(" %-28s %16.10f" %("Bare interaction energy:",int_totene))
    print(" %-28s %16.10f" %("Size-consistent correction:",scc))
    print(" %-28s %16.10f" %("SCC interaction energy:",scc_totene))
    print("")
