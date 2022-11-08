import numpy as np
import matplotlib.pyplot as plt
import json

import xtrack as xt
import xfields as xf
import xpart as xp
import xobjects as xo

context = xo.ContextCpu()

line_folder = 'Lines/run3_collisions_30cm_160urad_1.2e11_2.0um_62.310_60.320_15_430_0.001/'
#fname_line = 'line_bb_for_tracking_no_coupling.json'
# with open(line_folder + "line_bb_for_tracking.json", 'r') as fid:
#      input_data_tracking = json.load(fid)
with open(line_folder + "line_b1_tracking.json", 'r') as fid:
     input_data_b1 = json.load(fid)
with open(line_folder + "line_b4_tracking.json", 'r') as fid:
     input_data_b2 = json.load(fid)

# remove bb lenses
for key in input_data_b1["elements"].keys():
    if "BeamBeam" in input_data_b1["elements"][key]["__class__"]:
        input_data_b1["elements"][key] = xt.Drift(length=0.).to_dict()
for key in input_data_b2["elements"].keys():
    if "BeamBeam" in input_data_b2["elements"][key]["__class__"]:
        input_data_b2["elements"][key] = xt.Drift(length=0.).to_dict()

line_b1 = xt.Line.from_dict(input_data_b1)
line_b1.particle_ref = xp.Particles(p0c=input_data_b1["particle_on_tracker_co"]["p0c"])
line_b2 = xt.Line.from_dict(input_data_b2)
line_b2.particle_ref = xp.Particles(p0c=input_data_b2["particle_on_tracker_co"]["p0c"])

# line_b1.vars['vrf400'] = 12
# line_b2.vars['vrf400'] = 12
## line_b2.vars['lagrf400.b2'] = 0.5
## for ii, el in enumerate(line_b2.elements):
##     elcl = el.__class__
##     is_multi = isinstance(el, xt.beam_elements.elements.Multipole)
##     is_drift = isinstance(el, xt.beam_elements.elements.Drift)
##     is_cav = isinstance(el, xt.beam_elements.elements.Cavity)
##     is_diped = isinstance(el, xt.beam_elements.elements.DipoleEdge)
##     is_bb2d = isinstance(el, xf.beam_elements.beambeam2d.BeamBeamBiGaussian2D)
##     is_bb3d = isinstance(el, xf.beam_elements.beambeam3d.BeamBeamBiGaussian3D)
##     if is_cav:
##         print(el.to_dict())
##     if is_bb2d:
##         el.scale_strength=0
##     if is_bb3d:
##         el.scale_strength=0
##         line_b1.elements[ii] = xt.Drift(length=0.)

tracker_b1 = xt.Tracker(_context=context, line=line_b1)
tracker_b2 = xt.Tracker(_context=context, line=line_b2)
twiss_b1 = tracker_b1.twiss()
twiss_b2 = tracker_b2.twiss().reverse()

# parts = xp.Particles(p0c=6.8e12, x=6e-4, px=2.e-5, y=-0.003, py=-7e-5, zeta=np.linspace(-0.4, 0.4, 10))
# parts = xp.Particles(p0c=6.8e12, zeta=np.linspace(-0.4, 0.4, 10))
# tracker_b2.track(parts, num_turns=1000, turn_by_turn_monitor=True)
# plt.figure(100)
# plt.plot(tracker_b2.record_last_track.zeta.T, tracker_b2.record_last_track.delta.T, '.')
# plt.show()

n_particles = 1000000
distribution="gaussian"

std_zeta_list = []
std_delta_list = []

for rms_bunch_length in np.linspace(0.076, 0.100, 10):
    zeta, delta, matcher = xp.generate_longitudinal_coordinates(tracker=tracker_b1,
                                                     num_particles=n_particles,
                                                     sigma_z=rms_bunch_length, distribution=distribution,
                                                     engine="single-rf-harmonic", return_matcher=True)
                                                     #engine="pyheadtail", return_matcher=True)
    
    std_zeta_list.append(np.std(zeta))
    std_delta_list.append(np.std(delta))

for std_zeta, std_delta in zip(std_zeta_list, std_delta_list):
    print(f"Bunch length rms: {std_zeta:.4f}, delta rms: {std_delta:.7f}")

params = np.polyfit(std_zeta_list, std_delta_list, deg=2)
xp = np.linspace(0.075, 0.101, 10)
plt.plot(std_zeta_list, std_delta_list , 'o-')
plt.plot(xp, params[0] * xp**2 + params[1] * xp + params[2], 'r-')



plt.show()