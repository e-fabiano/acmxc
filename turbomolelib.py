import os
import tools
import subprocess
import shutil
from numpy import array as nparray
from numpy import sum as npsum
from numpy import power as nppower
from numpy import zeros as npzeros

import sys

class turbomole:

    def __init__(self,tdir,ncpu,wfunc,metal_mode,path):
        self.ncpu = ncpu
        self.wfunc = wfunc
        self.tdir = self.check_turbodir(tdir)
        self.set_turbomole()
        self.check_input(metal_mode,path)
        self.scf_file_name = "scf.out"
        self.w_file_name = "w.out"
        self.mp2_file_name = "mp2.out"


# check the value of TURBODIR and set it if not defined
#   if tdir was not given in command line set it to the default value ($TURBODIR)
    def check_turbodir(self,tdir):
        if (tdir == None):
            envtdir=os.getenv('TURBODIR')
            if (envtdir == None):
                tools.error("The $TURBODIR variable is not set!")
            tdir = envtdir
        return tdir


# set turbomole environment variables and executables paths
    def set_turbomole(self):
        para_arch = "SMP" if (self.ncpu > 1) else ""
        os.putenv("PARA_ARCH",para_arch)
        os.putenv("SMPCPUS",str(self.ncpu))
        sname = subprocess.check_output("sysname", shell=True).decode('utf-8').strip()
        ssname = subprocess.check_output("sysname -s", shell=True).decode().strip()
        if (self.ncpu > 1):
            path = f"{self.tdir}/scripts:{self.tdir}/bin/{sname}:{self.tdir}/bin/{ssname}:/bin:/usr/bin"
        else:
            path = f"{self.tdir}/scripts:{self.tdir}/bin/{sname}:/bin:/usr/bin"
        os.putenv("PATH",path)
        if (self.ncpu > 1):
            self.ridft_exec = f"{self.tdir}/bin/{sname}/ridft_smp"
            self.ricc2_exec = f"{self.tdir}/bin/{sname}/ricc2_smp"
        else:
            self.ridft_exec = f"{self.tdir}/bin/{sname}/ridft"
            self.ricc2_exec = f"{self.tdir}/bin/{sname}/ricc2"


    def check_input(self,metal_mode,path):
        control_file_path = f"{path}/control"        
        if (not os.path.isfile(control_file_path)):
            tools.error("control file not found!")
        with open(control_file_path,"r") as f:
            control_content = f.read()
            if "$rij" not in control_content:
                tools.error("RI option not set in control file")
            if (not metal_mode):
                if "$ricc2" not in control_content:
                    tools.error("$ricc2 datagroup not set in control file")
                if "mp2" not in control_content:
                    tools.error("MP2 option not set in control file")
                

    def check_scf_convergence(self):
        with open(self.scf_file_name,"r") as f:
            scf_content = f.readlines()
        return any("convergence criteria satisfied" in line for line in scf_content)


    def check_mp2_convergence(self):
        with open(self.mp2_file_name,"r") as f:
            mp2_content = f.readlines()
        return any("ricc2 : all done" in line for line in mp2_content)


    def run_ridft(self,f,rerun=False,damp=0.0):
        if (rerun):
            tools.kdg("restart")
            tools.kdg("scfdamp")
            tools.adg(f"$scfdamp   start={damp:.3f}  step=0.050  min=0.100")
            tools.kdg("scfiterlimit")
            tools.adg("$scfiterlimit 200")
        if (self.ncpu > 1):
            subprocess.run([self.ridft_exec,"-nthreads",f"{self.ncpu}"], stdout=f, stderr=subprocess.DEVNULL)
        else:
            subprocess.run([self.ridft_exec], stdout=f, stderr=subprocess.DEVNULL)

            
    def run_scf_and_w(self):      
#       do nothing if the SCF calculation is already there
        if (os.path.isfile(self.scf_file_name) and self.check_scf_convergence()):  
            print("Already done")
        else:
#       SCF part
            self.run_scf()
#       W_inf part
        self.run_winf()

        
    def run_scf(self):
        with open(self.scf_file_name,"w") as f:
            self.run_ridft(f,False)
            if (not self.check_scf_convergence()):
                self.run_ridft(f,True,damp=3)
                if (not self.check_scf_convergence()):
                    self.run_ridft(f,True,damp=16)
                    if (not self.check_scf_convergence()):
                        self.run_ridft(f,True,damp=12)
                        if (not self.check_scf_convergence()):
                            tools.error(f"The SCF calculation did not converge.\n Check the {self.scf_file_name} file")

                            
    def run_winf(self):
#       W calculation is performed in the temporary directory TMP
        if (os.path.isdir("TMP")): shutil.rmtree("TMP") 
        os.mkdir("TMP")
        tools.cp_all_files(".","TMP")
        os.chdir("TMP")
        tools.kdg("scfiterlimit")
        tools.adg("$scfiterlimit 1")
        tools.adg(f"$dft\n   functional {self.wfunc}\n   gridsize 3")
        with open(self.w_file_name,"w") as f:
            if (self.ncpu > 1):
                subprocess.run([self.ridft_exec,"-nthreads",f"{self.ncpu}"], stdout=f, stderr=subprocess.DEVNULL)
            else:
                subprocess.run([self.ridft_exec], stdout=f, stderr=subprocess.DEVNULL)
        os.chdir("..")
        shutil.copy(f"TMP/{self.w_file_name}",".")
        shutil.rmtree("TMP")# tools.remove_dir("TMP")


    def run_w34(self):
#       W34 calculation is performed in the temporary directory TMP
        if (os.path.isdir("TMP")): shutil.rmtree("TMP") 
        os.mkdir("TMP")
        tools.cp_all_files(".","TMP")
        os.chdir("TMP")
#       check if it is an UHF calculation
        uhf = False
        with open("control","r") as f:
            control_content = f.readlines()
            for line in control_content:
                if ("$uhf" in line):
                    uhf = True
                    break        
#       get the number of atoms from file control
        with open("control","r") as f:
            control_content = f.readlines()
            for line in control_content:
                if ("natoms=" in line):
                    natoms = int(line.split("=")[1])
                    break
#       get atomic coordinates and nuclear charges from the scf_file
        atomic_coords = []
        nuclear_charges = []
        with open(self.scf_file_name,"r") as f:
            scf_content = f.readlines()
            for i in range(len(scf_content)):
                line = scf_content[i]
                if ("atomic coordinates            atom    charge  isotop" in line):
                    for j in range(1,natoms+1):
                        linea = scf_content[i+j]
                        xx = float(linea.split()[0])
                        yy = float(linea.split()[1])
                        zz = float(linea.split()[2])
                        ch = float(linea.split()[4])
                        atomic_coords.append((xx,yy,zz))
                        nuclear_charges.append(ch)
                    break
#       compute densities at nuclear positions
        stri = "$pointval dens geo=point\n "
        for i in range(len(atomic_coords)):
            xx = atomic_coords[i][0]
            yy = atomic_coords[i][1]
            zz = atomic_coords[i][2]
            stri = f"{stri} {xx} {yy} {zz}\n "
        stri = stri[:-2]
        tools.adg(stri)
        subprocess.run([self.ridft_exec,"-proper"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        shutil.copy("td.xyz","../densities.dat")
        if (uhf):
            shutil.copy("sd.xyz","../spindensities.dat")
#            
        os.chdir("..")
        shutil.rmtree("TMP")
#       get densities from the file densities.dat
        densities = []
        with open("densities.dat","r") as f:
            densities_content = f.readlines()
            for i in range(len(densities_content)):
                line = densities_content[i]
                if ("cartesian coordinates x,y,z and f(x,y,z)" in line):
                    for j in range(1,natoms+1):
                        linea = densities_content[i+j]
                        densities.append(float(linea.split()[3]))
                    break
#       get spin densities (if uhf)
        if (uhf):
            spindens = []
            with open("spindensities.dat","r") as f:
                spindens_content = f.readlines()
                for i in range(len(spindens_content)):
                    line = spindens_content[i]
                    if ("cartesian coordinates x,y,z and f(x,y,z)" in line):
                        for j in range(1,natoms+1):
                            linea = spindens_content[i+j]
                            spindens.append(float(linea.split()[3]))
#       compute W3/4
        nuclear_charges = nparray(nuclear_charges)
        densities = nparray(densities)
        if (uhf):
            spindens = nparray(spindens)
        else:
            spindens = npzeros(natoms,dtype=float)
        zeta = spindens/densities
        ss = 0.5*(1.+zeta*zeta)
        e34 = -2.002 - 1.588*ss + 0.394*ss*ss
        self.w34ene = 0.470698131888*npsum(e34*nuclear_charges*nppower(densities,(1./4.))) 
        
        
    def run_mp2(self):
#       do nothing if the MP2 calculation is already there
        if (os.path.isfile(self.mp2_file_name)):
           if (self.check_mp2_convergence()):          
               print("Already done")
               return
#       MP2 calculation
        with open(self.mp2_file_name,"w") as f:
            if (self.ncpu > 1):
                subprocess.run([self.ricc2_exec,"-nthreads",f"{self.ncpu}"], stdout=f, stderr=subprocess.DEVNULL)
            else:
                subprocess.run([self.ricc2_exec], stdout=f, stderr=subprocess.DEVNULL)
        if (not self.check_mp2_convergence()):
            tools.error(f"The MP2 calculation did not converge.\n Check the {self.mp2_file_name} file")


# for metallic systems print a simple file with E_mp2=-inf
    def print_metallic_mp2_file(self):
        with open(self.mp2_file_name,"w") as f:
            f.write("MP2 correlation energy (doubles) . . -inf")


    def extract_results(self):
#       scf energy ( = HF energy)
        with open("energy","r") as f:
            raw = f.readlines()
        self.scfene = float(raw[len(raw)-2].split()[1])
#       exchange energy, W_inf and W'_inf        
        with open(self.w_file_name,"r") as f:
            wfile_content = f.readlines()
            found = 0
            for i, line in enumerate(wfile_content):
                if ("exK =" in line):
                    self.xene = float(line.split("=", 1)[1])
                    found += 1
                elif ("RI-K =" in line):
                    self.xene = float(line.split("=", 1)[1])
                    found += 1
                elif ("Strong correlation functionals" in line):
                    self.wene = float(wfile_content[i + 1].split()[2])
                    self.w1ene = float(wfile_content[i + 2].split()[2])
                    found += 1
                if (found == 2):
                    break                
#       MP2 energy
        with open(self.mp2_file_name,"r") as f:
            mp2file_content = f.readlines()
            for line in mp2file_content:
                if ("MP2 correlation energy (doubles)" in line):
                    self.mp2ene = float(line.split()[6])
                    break
