import os
import shutil
import numpy as np

num_turns = 1000000
num_particles = 20000
sigma = 6.5
theta_x = 0
theta_y = 0
zeta_norm = 0
pzeta_norm = np.linspace(0,8.2e-4, 10)
pzeta_norm = pzeta_norm.tolist()
repetition_period = 10000
optimization = 4
ecloud = ['no','mb','mq','all']
sey = 1.30
ppb = 1.20


parent_dir = '/afs/cern.ch/user/j/jpotdevi/public/simulations/DA/'

for i in range(len(ecloud)):
    for j in range(len(pzeta_norm)):

        directory = "DA_%secloud_sey%.2f_%.2fe11ppb_pzeta%.2fe-4"%(ecloud[i], sey, ppb, pzeta_norm[j]*1e4)
        print(directory)

        os.mkdir(directory)

        with open ('executable.sh', 'w') as ex:
            ex.write('''\
#! /bin/bash

source /usr/local/xsuite/miniforge3/bin/activate xsuite

nvidia-smi

pwd


cp /afs/cern.ch/user/k/kparasch/work/public/for_josephine/ecloud_xsuite_filemanager.py .
cp /afs/cern.ch/user/j/jpotdevi/public/simulations/DA/dynamic_aperture.py .
cp /afs/cern.ch/user/j/jpotdevi/public/simulations/line_and_particle.json .
cp /afs/cern.ch/user/j/jpotdevi/public/simulations/eclouds.json .
cp /afs/cern.ch/user/j/jpotdevi/public/simulations/collimators.py .

''')
    
            ex.write('''\
xrdcp root://eosuser.cern.ch//eos/project/e/ecloud-simulations/kparasch/refined_Pinches/MB/refined_LHC_MB_450GeV_sey%.2f_%.2fe11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 EC_MB.h5
xrdcp root://eosuser.cern.ch//eos/project/e/ecloud-simulations/kparasch/refined_Pinches/MQF/refined_LHC_MQF_450GeV_sey%.2f_%.2fe11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 EC_MQF.h5
xrdcp root://eosuser.cern.ch//eos/project/e/ecloud-simulations/kparasch/refined_Pinches/MQD/refined_LHC_MQD_450GeV_sey%.2f_%.2fe11ppb_symm2D_MTI2.0_MLI2.0_DTO1.0_DLO1.0.h5 EC_MQD.h5 

'''%(sey, ppb, sey, ppb, sey,ppb))
    
            ex.write('''\
python dynamic_aperture.py --filename %s.h5 --num_turns %d --num_particles %d --sigma %f --theta_x %f --theta_y %f --zeta_norm %f --pzeta_norm %f --optimization %d --ecloud %s --sey %f --ppb %f

'''%(directory, num_turns, num_particles, sigma, theta_x, theta_y, zeta_norm, pzeta_norm[j], optimization, ecloud[i], sey, ppb))

            ex.write('''\
cp %s.h5 /eos/project/e/ecloud-simulations/jpotdevi/simulations/DA/'''%(directory))


        shutil.move('executable.sh', directory)

        with open ('gpu_submit.sub','w') as sub:
                sub.write('''\
executable = executable.sh
arguments = $(ClusterId).$(ProcId)
output = $(ClusterId).$(ProcId).out
error = $(ClusterId).$(ProcId).err
log = $(ClusterID).$(ProcId).log
transfer_output_files   = ""
+SingularityImage = "/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/kparasch/incoherent-ecloud-docker-image:latest"
on_exit_remove = (ExitBySignal == False) && (ExitCode == 0)
max_retries = 3
requirements = Machine =!= LastRemoteHost && regexp("V100", Target.CUDADeviceName) && (regexp("9.0",TARGET.CondorVersion))
request_GPUs = 1
request_CPUs = 1
+MaxRunTime = 82800
queue
''')

        shutil.move('gpu_submit.sub', directory)


        os.system('cd '+directory+'; condor_submit gpu_submit.sub ; cd ..')
