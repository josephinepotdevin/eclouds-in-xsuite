import xtrack as xt
import xobjects as xo
import xpart as xp
import xfields as xf

import collimators

import json
import ecloud_xsuite_filemanager as exfm

import numpy as np
import time

import argparse

start_running = time.time()

parser = argparse.ArgumentParser()
parser.add_argument('--filename', nargs='?', default='delete.h5', type=str)
parser.add_argument('--num_turns', nargs='?', default=10000, type=int)
parser.add_argument('--num_particles', nargs='?', default=20000, type=int)
parser.add_argument('--sigma', nargs='?', default=6, type=float)
parser.add_argument('--theta_x', nargs='?', default=0, type=float)
parser.add_argument('--theta_y', nargs='?', default=0, type=float)
parser.add_argument('--zeta_norm', nargs='?', default=0, type=float)
parser.add_argument('--pzeta_norm', nargs='?', default=8.31e-4, type=float)
parser.add_argument('--optimization', nargs='?', default=4, type=int)
parser.add_argument('--ecloud', nargs='?', default="no", type=str)
parser.add_argument('--sey', nargs='?', default=1.30, type=float)
parser.add_argument('--ppb', nargs='?', default=1.20, type=float)
args = parser.parse_args()

output_filename = args.filename
num_turns = args.num_turns
num_particles = args.num_particles
num_r = int(np.sqrt(num_particles))
num_theta = int(np.sqrt(num_particles))
sigma_init = args.sigma
theta_x = args.theta_x
theta_y = args.theta_y
zeta_norm = args.zeta_norm
pzeta_norm = args.pzeta_norm
optimization = args.optimization
ecloud = args.ecloud
sey = args.sey
ppb = args.ppb

with open('eclouds.json', 'r') as fid:
    ecloud_info = json.load(fid)

with open('line_and_particle.json', 'r') as fid:
    input_data = json.load(fid)
line = xt.Line.from_dict(input_data["line"])
line.particle_ref = xp.Particles.from_dict(input_data['particle'])

context = xo.ContextCupy()

tau_max = 0.4
collimators.collimator_setup(line, context=context)

EC_MB = "EC_MB.h5"
EC_MQF = "EC_MQF.h5"
EC_MQD = "EC_MQD.h5"


if ecloud == "no":
    filenames = {}
if ecloud == "mb":
    filenames = {"mb":EC_MB}
if ecloud == "mqf":
   filenames = {"mqf":EC_MQD}
if ecloud == "mqd":
   filenames = {"mqd":EC_MQF}
if ecloud == "mq":
    filenames = {"mqf":EC_MQF, "mqd":EC_MQD}
if ecloud == "all":
    filenames = {"mb":EC_MB, "mqf":EC_MQF, "mqd":EC_MQD}

if optimization > 0:
    line.remove_inactive_multipoles()
    if optimization > 1:
        line.remove_zero_length_drifts()
        if optimization >2:
            line.merge_consecutive_drifts()
            if optimization >3:
                line.merge_consecutive_multipoles()
                if optimization >4:
                    line.use_simple_bends()
                    if optimization>5:
                        line.use_simple_quadrupoles()

tracker, twiss_without_ecloud, twiss_with_ecloud = xf.full_electroncloud_setup(line=line, 
ecloud_info=ecloud_info, filenames=filenames, context=context, tau_max=tau_max)

r = np.linspace(0, sigma_init, num_r+1)[1:]
theta = np.linspace(0,np.pi/2, num_theta+1, endpoint = True)[1:]
rr, tcos= np.meshgrid(r, np.cos(theta))
rr, tsin = np.meshgrid(r, np.sin(theta))
Ax_norm = (rr*tcos).flatten()
Ay_norm = (rr*tsin).flatten()
x_norm, y_norm = Ax_norm*np.cos(theta_x), Ay_norm*np.cos(theta_y)
px_norm, py_norm = 0,0

particles = xp.build_particles(_context=context, tracker=tracker, x_norm=x_norm, y_norm=y_norm, px_norm = px_norm, py_norm = py_norm,
scale_with_transverse_norm_emitt=(2.0e-6, 2.0e-6), 
zeta=zeta_norm, delta=pzeta_norm)

x_init = context.nparray_from_context_array(particles.x)
px_init = context.nparray_from_context_array(particles.px)
y_init = context.nparray_from_context_array(particles.y)
py_init = context.nparray_from_context_array(particles.py)
zeta_init = context.nparray_from_context_array(particles.zeta)
pzeta_init = context.nparray_from_context_array(particles.pzeta)

at_turn_init = context.nparray_from_context_array(particles.at_turn)
id_init = context.nparray_from_context_array(particles.particle_id)
state_init = context.nparray_from_context_array(particles.state)

print('Build tracker...')

start_tracking = time.time()
tracker.track(particles, num_turns=num_turns)
context.synchronize()
end_tracking = time.time()
print(f'Tracking time:{(end_tracking-start_tracking)/60.:.4f}mins')


x_fin = context.nparray_from_context_array(particles.x)
px_fin = context.nparray_from_context_array(particles.px)
y_fin = context.nparray_from_context_array(particles.y)
py_fin = context.nparray_from_context_array(particles.py)
zeta_fin = context.nparray_from_context_array(particles.zeta)
pzeta_fin = context.nparray_from_context_array(particles.pzeta)

at_turn_fin =context.nparray_from_context_array(particles.at_turn)
id_fin = context.nparray_from_context_array(particles.particle_id)
state_fin = context.nparray_from_context_array(particles.state)

end_running = time.time()

IC_norm = {"Ax_norm":Ax_norm, "Ay_norm":Ay_norm, "x_norm":x_norm, 
"px_norm":px_norm, "y_norm":y_norm, "py_norm":py_norm, 
"zeta_norm":zeta_norm, "pzeta_norm":pzeta_norm}

IC_phys = {"x_init":x_init, "px_init":px_init, "y_init":y_init, 
"py_init":py_init, "zeta_init":zeta_init, "pzeta_init":pzeta_init, 
"at_turn_init":at_turn_init, "id_init":id_init, "state_init":state_init}

tbt_phys = {"x_fin":x_fin, "px_fin":px_fin, "y_fin":y_fin, "py_fin":py_fin, 
"zeta_fin":zeta_fin, "pzeta_fin":pzeta_fin, "at_turn_fin":at_turn_fin, 
"id_fin":id_fin, "state_fin":state_fin}

Param = {'num_turns':num_turns, "num_r":num_r, "num_theta":num_theta, "sigma_init":sigma_init, "time_track":(end_tracking-start_tracking)/60, "time_run":(end_running-start_running)/60}

Part_ref = line.particle_ref.to_dict()

exfm.dict_to_h5(IC_norm, output_filename, group='IC_norm', readwrite_opts='w')
exfm.dict_to_h5(IC_phys, output_filename, group='IC_phys', readwrite_opts='a')
exfm.dict_to_h5(tbt_phys, output_filename, group='tbt_phys', readwrite_opts='a')
exfm.dict_to_h5(Param, output_filename, group='parameters', readwrite_opts='a')
exfm.dict_to_h5(Part_ref, output_filename, group='particle_ref', readwrite_opts='a')
