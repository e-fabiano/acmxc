
import sys
import shutil
import os
from datetime import datetime
import acmlib




# print an error message and exit
def error(stri):
    print()
    print("ERROR!!")
    print(stri)
    print()
    sys.exit(1)


# delete a data group from turbomole control file
def kdg(datagroup):
    stri = "$"+datagroup
    f = open("control","r")
    control_content = f.readlines()
    f.close()
    f = open("control.tmp","w")
    i = 0
    while (i < len(control_content)):
        line = control_content[i]
        if (stri not in line):
            f.write(line)
#            print(line)
            i = i + 1
        else:
            j = 1
            while (i+j < len(control_content)):
                bline = control_content[i+j]
                if ("$" in bline):
                    break
                j = j + 1
            i = i + j
    f.close()
    shutil.move("control.tmp","control")
                

# add a data group to turbomole control file
def adg(datagroup):
    stri = datagroup + "\n"
    kdg("end")
    f = open("control","a")
    f.write(stri)
    f.write("$end\n")
    f.close()


# copy all file (and only files) from one directory to another
def cp_all_files(from_dir,to_dir):
    for ff in os.listdir(from_dir):
        file_full_name = os.path.join(from_dir, ff)
        if os.path.isfile(file_full_name):
            shutil.copy(file_full_name, to_dir)

# remove a non empty directory
def remove_dir(dire):
    for ff in os.listdir(dire):
        file_full_name = os.path.join(dire, ff)
        if os.path.isfile(file_full_name):
            os.remove(file_full_name)
    os.rmdir(dire)


# print header
def print_header():
    data = datetime.now().strftime("%d.%m.%Y - %H:%M:%S")
    print ("")
    print("             +-----------------------------+")
    print("             |            ACMxc            |")
    print("             | Perform HF-ACM calculations |")
    print("             |                             |")
    print("             +-----------------------------+")
    print("")
    print(data)


# print a summary of the options
def print_options(options):
    program = options[0]
    tdir = options[1]
    baseinput = options[2]
    ncpu = f"{options[3]}"
    formula = options[4]
    wfunc = options[5]
    rerun = options[6]
    metal_mode = f"{options[7]}"
    w34 = f"{options[8]}"
    print("")
    print("Program: ",program)
    if (program == "crystal"):
        print ("Input file: ",baseinput)
    print("Program's path: ",tdir)
    print("N. cpu: ",ncpu)
    print("ACM model: ",formula)
    acmlib.print_refs(formula)
    print("W functional: ",wfunc)
    print_w_refs(wfunc)
    if (program == "crystal"):
        print ("Rerun: ",rerun)
    print("W3/4 functional: ",w34)
    print("Metal mode: ",metal_mode)
    print("")


def print_w_refs(wfunc):
    if (wfunc == "pc"):
        refstri = "Phys. Rev. A 62, 012502 (2000)"
    elif (wfunc == "hpc"):
        refstri = "J. Chem. Theory Comput. 18, 5936 (2022)"
    elif (wfunc == "mpc"):
        refstri = "Phys. Rev. B 99, 085117 (2019)"
    elif (wfunc == "hfpc"):
        refstri = "J. Phys. Chem. Lett xx, yyy, (2025)"
    elif (wfunc == "lda"):
        refstri = "Phys. Rev. 140, A1133 (1965)"
    else:
        refstri = ""
    stri = "  Ref: " + refstri
    print(stri)
