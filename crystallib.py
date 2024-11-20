import os
import tools
import subprocess
import shutil

class crystal:

    allowed_wfunc = ["pc", "hpc", "lda"]

    
    def __init__(self,tdir,ncpu,wfunc,rerun,baseinput,metal_mode,path):
        self.ncpu = ncpu
        self.set_wfunc(wfunc)
        self.tdir = self.check_crystal_dir(tdir)
        self.set_crystal()
        self.baseinput = baseinput
        self.check_input(metal_mode,path)
        self.scf_file_name = "scf.out"
        self.w_file_name = "w.out"
        self.rerun = rerun
        self.guessp = False
        self.mp2_file_name = "mp2.out"


    def set_wfunc(self,wfunc):
        if (wfunc in self.allowed_wfunc):
            self.wfunc = wfunc
        else:
            tools.error(f"The {wfunc} W-functional is not supported in crystal")


# check what is the executable directory for crystal
#   if tdir is not given in command line try to set it from the default value ($CRY23_EXEDIR)
    def check_crystal_dir(self,tdir):
        if (tdir == None):
            envdir = os.getenv('CRY23_UTILS')
            if (envdir == None):
                tools.error("The $CRY23_EXEDIR variable is not set!")
            tdir = envdir
        return tdir


# set crystal environment variables and executables
    def set_crystal(self):
        self.crystal_exec = self.tdir + "/runcry23OMP"
        self.properties_exec = self.tdir + "/runprop23"
        self.cryscor_exec = self.tdir + "/runcryscor"


    def check_input(self,metal_mode,path):
        input_file_path = f"{path}/{self.baseinput}"
        if (not os.path.isfile(f"{input_file_path}.d12")):
            tools.error(f"{input_file_path}.d12 not found!")
        if (not os.path.isfile(f"{input_file_path}.d3")):
            tools.error(f"{input_file_path}.d3 not found!")
        if (not metal_mode):
            if (not os.path.isfile(f"{input_file_path}.d4")):
                tools.error(f"{input_file_path}.d4 not found!")
            with open(self.baseinput +".d12", 'r') as f:
                content = f.read()
                if ("UHF" in content):
                    tools.error("UHF not implemented in cryscor")


    def check_scf_convergence(self):
        with open(self.scf_file_name,"r") as f:
            scf_content = f.readlines()
        return any("SCF ENDED" in line for line in scf_content)


    def check_mp2_convergence(self):
        with open(self.mp2_file_name,"r") as f:
            mp2_content = f.readlines()
        return any("MP2 CONVERGENCE REACHED" in line for line in mp2_content) 

    
    def run_scf_and_w(self):
        if (os.path.isfile(self.scf_file_name) and self.check_scf_convergence() and self.rerun == False):
            print("already done")
        else:
#       SCF part
            self.run_scf()
#       W_inf part
        self.run_winf()


    def run_scf(self):
        with open(self.baseinput +".d12", 'r') as file: 
             content = file.read()     
             if "GUESSP" in content:
                self.guessp = True
        if (self.guessp == False):
            subprocess.run([self.crystal_exec, str(self.ncpu), self.baseinput], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        else:
            subprocess.run([self.crystal_exec, str(self.ncpu), self.baseinput, self.baseinput], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        shutil.move(self.baseinput + ".out", self.scf_file_name)
        if(not self.check_scf_convergence()):
            tools.error(f"The SCF calculation did not converge.\n Check the {self.scf_file_name} file")


    def run_winf(self):             
        subprocess.run([self.properties_exec, self.baseinput, self.baseinput], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        shutil.move(self.baseinput + ".outp", self.w_file_name)
        

    def run_w34(self):
        tools.error("W_3/4 is not yet implemented for Crystal")
        

    def run_mp2(self):
        if(os.path.isfile(self.mp2_file_name) and self.rerun == False):
            print("already done")
            return
        subprocess.run([self.cryscor_exec, self.baseinput], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        shutil.move(self.baseinput + ".outc", self.mp2_file_name)
            

# for metallic systems print a simple file with E_mp2=-inf
    def print_metallic_mp2_file(self):
        with open(self.mp2_file_name,"w") as f:
            f.write("MP2 CORRELATION ENERGY: . -inf")


    def extract_results(self):
#       scf energy (= HF energy) and exchange energy
        with open(self.scf_file_name,"r") as f:
            scffile_content = f.readlines()
            found = 0
            for line in scffile_content:
                if ("::: TOTAL   ENERGY" in line):
                    self.scfene = float(line.split()[3])
                    found = found + 1
                if ("FOCK EXCHANGE ENERGY" in line):
                    self.xene = float(line.split()[3])
                    found = found + 1
                if (found == 2):
                    break
#       W_inf and W'_inf
        with open(self.w_file_name,"r") as f:
            wfile_content = f.readlines()
            if (self.wfunc == "pc"):
                stri1 = "AC-Winf "
                stri2 = "AC-W1inf "
            elif (self.wfunc == "hpc"):
                stri1 = "AC-Winf-HPC"
                stri2 = "AC-W1inf-HPC"
            elif  (self.wfunc == "lda"):          
                stri1 = "PERDEW-WANG LSD"         
                stri2 = "PERDEW-WANG GGA"
            found = 0
            for line in wfile_content:
                if (stri1 in line):
                    nx = 1
                    if (self.wfunc == "lda" ):    
                        nx=2
                    self.wene = float(line.split()[nx])
                    found = found + 1
                if (stri2 in line):
                    nx = 1
                    if (self.wfunc == "lda" ):    
                        nx=2
                    self.w1ene = float(line.split()[nx])
                    found = found + 1
                if (found == 2):
                    break
#       MP2 energy
        with open(self.mp2_file_name,"r") as f:
            mp2file_content = f.readlines()
            for line in mp2file_content:
                if ("MP2 CORRELATION ENERGY:" in line):
                    self.mp2ene = float(line.split()[4])
                    break
