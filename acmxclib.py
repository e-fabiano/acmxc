import tools
import os
import turbomolelib
import crystallib
import acmlib

import sys

class acmxc:

    allowed_programs_list = ["turbomole","crystal"]
    allowed_acm_formulas = ["isi","revisi","spl","lb","spl2","mpacf1","genisi","dpi","hfac24","mp2","hflda"]
    allowed_wfunc = ["pc","hpc","mpc","hfpc","lda"]

    

    def __init__(self, program="turbomole", tdir=None, path=".", prog_input="input", ncpu=1, formula="isi", wfunc=None, rerun=False, metal_mode=False, verbose=True):
        self.verbose = verbose
        if (self.verbose): tools.print_header()
        self.set_tdir(tdir)
        self.set_path(path)
        self.set_prog_input(prog_input)
        self.set_ncpu(ncpu)
        self.set_formula(formula)
        self.set_wfunc(wfunc,formula)
        self.rerun = rerun
        self.set_w34(formula)
        self.set_metal_mode(metal_mode)
        self.set_program(program)
        if (self.verbose): tools.print_options([self.program_name,self.tdir,self.baseinput,self.ncpu,self.acm_formula,self.wfunc,self.rerun,self.metal_mode,self.w34])
        

    def set_tdir(self,tdir):    
        if (tdir is not None):
            if (not os.path.isdir(tdir)):
                tools.error(f"{tdir} is not a valid directory for option <tdir>")
            else:
                self.tdir = tdir
        else:
            self.tdir = None

            
    def set_path(self,path):
        if (not os.path.isdir(path)):
            tools.error(f"{path} is not a valid value for <path>")
        self.path = path
            
            
    def set_prog_input(self,prog_input):
        self.baseinput = prog_input
        
        
    def set_ncpu(self,ncpu):
        try:
            ncpu = int(ncpu)
            if (ncpu > 0):
                self.ncpu = ncpu
            else:
                tools.error(f"{ncpu} must be > 0")
        except ValueError:
            tools.error(f"{ncpu} is not a valid integer number. <ncpu> must be integer")

            
    def set_formula(self,formula):
        if (formula in self.allowed_acm_formulas):
            self.acm_formula = formula
        else:
            tools.error(f"{formula} is not a valid input for <formula>")

                
    def set_wfunc(self,wfunc,formula):
        if (wfunc == None):
            self.wfunc = acmlib.autoset_wfunc(formula)
        else:
            if (wfunc in self.allowed_wfunc):
                self.wfunc = wfunc
            else:
                tools.error(f"{wfunc} is not a valid input for <wfunc>")


    def set_w34(self,formula):
        if (formula == "hfac24"):
            self.w34 = True
        else:
            self.w34 = False
            

    def set_metal_mode(self,mode):
        if isinstance(mode, bool):
            self.metal_mode = mode
        else:
            tools.error(f"{mode} is not a valid input for <metal_mode>; it must be boolean")


    def set_program(self,program):
        if (program in self.allowed_programs_list):
            self.program_name = program
            if (self.program_name == "turbomole"):
                self.program = turbomolelib.turbomole(self.tdir,self.ncpu,self.wfunc,self.metal_mode,self.path)
            elif (self.program_name == "crystal"):
                self.program = crystallib.crystal(self.tdir,self.ncpu,self.wfunc,self.rerun,self.baseinput,self.metal_mode,self.path)
            self.set_tdir(self.program.tdir)
        else:
            tools.error(f"{program} is not a valid option for <program>")


    def run_program(self):
        here = os.getcwd()
        os.chdir(self.path)
#       SCF and W_inf part
        print("Running SCF calculation...")
        self.program.run_scf_and_w()
#       W_3/4 part
        self.program.w34ene = 0.0
        if (self.w34):
            print("Running W_3/4 calculation...")
            self.program.run_w34()
#       MP2 part
        if (not self.metal_mode):
            print("Running MP2 calculation...")
            self.program.run_mp2()
        else:
            self.program.print_metallic_mp2_file()
        os.chdir(here)
        

    def extract_results(self):
        here = os.getcwd()
        os.chdir(self.path)
        self.program.extract_results()
        self.scfene = self.program.scfene
        self.xene = self.program.xene
        self.wene = self.program.wene
        self.w1ene = self.program.w1ene
        self.w34ene = self.program.w34ene
        self.mp2ene = self.program.mp2ene
        os.chdir(here)
        if (self.verbose):
            print()
            print(" %-20s %16.10f" %("SCF energy:",self.scfene))
            print(" %-20s %16.10f" %("HF-exchange energy:",self.xene))
            if (self.acm_formula == "hfac24"):
                print(" %-20s %16.10f" %("E_el energy:",self.wene))
                print(" %-20s %16.10f" %("W_1/2 energy:",self.w1ene))
                print(" %-20s %16.10f" %("W_{c,inf}^HF energy:",self.wene+self.xene))
            else:
                print(" %-20s %16.10f" %("W_inf energy:",self.wene))
                print(" %-20s %16.10f" %("W1_inf energy:",self.w1ene))
            if (self.w34):
                print( " %-20s %16.10f" %("W_3/4 energy:",self.w34ene))
            print(" %-20s %16.10f" %("MP2 corr. energy:",self.mp2ene))


    def compute_acm_xc_energy(self):
        self.correne = acmlib.compute_acm(self.acm_formula,self.xene,self.wene,self.w1ene,self.w34ene,self.mp2ene)
        self.xcene = self.xene + self.correne
        self.totene = self.scfene + self.correne


    def print_results(self):
        print()
        print("RESULTS:")
        print(" Correlation energy:          %16.10f" %(self.correne))
        print(" Exchange-correlation energy: %16.10f" %(self.xcene))
        print(" Total energy:                %16.10f" %(self.totene))
        print()
