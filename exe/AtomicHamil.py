#!/usr/bin/env python
#-*- coding:utf-8 -*-
# script to submit the job
# execute ./AtHamil.py
import os
import sys
import subprocess
from collections import OrderedDict
exe = 'AtHamil.exe'

def set_input(params,zeta):
    """
    Generate the input file to run AtHamil. Edit parameters below. 
    Available Hamiltonians are:
        Coulomb     Nonrelativistic Coulomb Hamiltonian with one-body and two-body interaction
        Breit_1b    Coulomb Hamiltonian with one-body leading-order relativistic corrections. This is a 
                    much cheaper calculation than the full Breit Hamiltonian. Laguerre orbitals only. 
        Breit       Breit Hamiltonian with all O(alpha^2) corrections. Currently work in progress. 

    The matrix elements of isotope shift operators can also be generated.
        NormalMassShift     p_i \cdot p_i, proportional to the leading-order normal mass shift. 
        SpecificMassShift   \sum_{i \neq j} p_i \cdot p_j, proportional to the leading-order specific
                            mass shift. Two-body portion of the mass shift. 
        MassShift           Sum of normal and specific mass shifts.
        FieldShift          Leading-order field shift operator. 
    """
    # basic parameters
    params["basis"] = "LO"
    params["zeta"] = zeta
    params['emax'] = 6
    params['e2max'] = 12
    params['opname'] = 'Coulomb'
    params["file_name"]="Coulomb_Hamil_"+params["basis"]+str(params["zeta"])+\
            "_emax"+str(params["emax"])+"_e2max"+str(params["e2max"])
    if( "lmax" in params): params["file_name_nn"]+= "_lmax"+str(params["lmax"])
    params["file_name"]+=".snt"

def gen_script(params, batch, machine):
    fbase = "aHamil_"+os.path.splitext(params["file_name"])[0]
    file_input = "Input_" + fbase + ".dat"
    file_log   = "log_" + fbase + ".dat"
    fsh = "run_" + fbase + ".sh"
    prt = ""
    if(machine=="oak"):
        prt += "#!/bin/bash \n"
        prt += "#PBS -q oak \n"
        prt += "#PBS -l mem=128gb,nodes=1:ppn=32,walltime=072:00:00 \n"
        prt += "cd $PBS_O_WORKDIR\n"
    if(machine=="cedar"):
        prt = "#!/bin/bash\n"
        prt += "#SBATCH --account="+account+"\n"
        prt += "#SBATCH --nodes=1\n"
        prt += "#SBATCH --ntasks=1\n"
        prt += "#SBATCH --cpus-per-task=1\n"
        prt += "#SBATCH --mem=125G\n"
        prt += "#SBATCH --time=3-00:00\n\n"

    prt += 'echo "run ' +fsh + '..."\n'
    prt += "cat > "+file_input + " <<EOF\n"
    prt += "&input\n"
    for key, value in params.items():
        if(isinstance(value, str)):
            prt += str(key) + '= "' + str(value) + '" \n'
            continue
        if(isinstance(value, list)):
            prt += str(key) + "= "
            for x in value[:-1]:
                prt += str(x) + ", "
            prt += str(value[-1]) + "\n"
            continue
        prt += str(key) + '=' + str(value) + '\n'
    prt += "&end\n"
    prt += "EOF\n"
    if(batch):
        prt += exe + " " + file_input + " > " + file_log + " 2>&1\n"
        prt += "rm " + file_input + "\n"
    if(not batch):
        prt += exe + " " + file_input + "\n"
        prt += "rm " + file_input + "\n"
    f = open(fsh, "w")
    f.write(prt)
    f.close()
    os.chmod(fsh, 0o755)
    return fsh

def main(machinename=None):
    if(machinename==None):
        batch = False
        machine = 'local'
    if(machinename!=None):
        batch=True
        if(machinename.lower() == "local"):
            machine = 'local'
        if(machinename.lower() == "oak"):
            machine = "oak"
        if(machinename.lower() =="cedar"):
            machine = "cedar"

    #for zeta in [1,2,4,6,8,10,12,14,16]:
    for zeta in [1]:
        params = OrderedDict()
        set_input(params,zeta)
        fsh = gen_script(params, batch, machine)
        if(machine == 'local'):
            cmd = "./" + fsh
        if(machine=="oak"):
            cmd = "qsub " + fsh
        if(machine=="cedar"):
            cmd = "srun " + fsh
        subprocess.call(cmd,shell=True)

if(__name__ == "__main__"):
    if(len(sys.argv) == 1):
        main()
    if(len(sys.argv) > 1):
        main(sys.argv[1])

