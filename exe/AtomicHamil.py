#!/usr/bin/env python
#-*- coding:utf-8 -*-
# script to submit the job
# execute ./NuHamil.py $machinename
# $machinename is optional argument. If $machinename is not given, job will be submitted as interactive.
import os
import sys
import subprocess
from collections import OrderedDict
exe = 'NuHamil.exe'
account="rrg-holt"     # for cedar
account="rrg-navratil" # for cedar

def set_input(params,zeta):
    # basic parameters
    params['rank'] = 2
    params["is_atomic"] = True
    params["basis"] = "LO"
    params["zeta"] = zeta
    params['emax'] = 10
    params['e2max'] = 20
    params['e3max'] = 6
    params["file_name_nn"]="Coulomb_Hamil_"+params["basis"]+str(params["zeta"])+\
            "_emax"+str(params["emax"])+"_e2max"+str(params["e2max"])
    if( "lmax" in params): params["file_name_nn"]+= "_lmax"+str(params["lmax"])
    params["file_name_nn"]+=".snt"

def gen_script(params, batch, machine):
    fbase = "vHamil_"+os.path.splitext(params["file_name_nn"])[0]
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
        header = "#!/bin/bash\n"
        header += "#SBATCH --account="+account+"\n"
        header += "#SBATCH --nodes=1\n"
        header += "#SBATCH --ntasks=1\n"
        header += "#SBATCH --cpus-per-task=1\n"
        header += "#SBATCH --mem=125G\n"
        header += "#SBATCH --time=3-00:00\n\n"

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

    for zeta in [1,2,4,6,8,10,12,14,16]:
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

