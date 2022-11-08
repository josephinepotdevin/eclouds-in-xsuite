import numpy as np
import matplotlib.pyplot as plt
import json
import scipy.interpolate
import pickle as pkl

import xtrack as xt
import xfields as xf
import xpart as xp
import xobjects as xo


nemitt_x = 2.0e-6
nemitt_y = 2.0e-6
std_zeta = 0.090
def get_delta_std(zeta_std):
    params = [-4.64975248e-03,  1.62153983e-03, -9.61004104e-06]
    return params[0] * zeta_std**2 + params[1] * zeta_std + params[2]
std_delta =  get_delta_std(std_zeta)

def disable_beambeam(line):
    for ii, el in enumerate(line.elements):
        elcl = el.__class__
        is_multi = isinstance(el, xt.beam_elements.elements.Multipole)
        is_drift = isinstance(el, xt.beam_elements.elements.Drift)
        is_cav = isinstance(el, xt.beam_elements.elements.Cavity)
        is_diped = isinstance(el, xt.beam_elements.elements.DipoleEdge)
        is_bb2d = isinstance(el, xf.beam_elements.beambeam2d.BeamBeamBiGaussian2D)
        is_bb3d = isinstance(el, xf.beam_elements.beambeam3d.BeamBeamBiGaussian3D)
        if is_bb2d:
            el.scale_strength=0
        if is_bb3d:
            el.scale_strength=0

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


tracker_b1 = xt.Tracker(_context=context, line=line_b1)
tracker_b2 = xt.Tracker(_context=context, line=line_b2)
twiss_b1 = tracker_b1.twiss()
twiss_b2 = tracker_b2.twiss().reverse()

start_IR1 = 19725
end_IR1 = 20263
ip1_s = 19994.1624
start_IR5 = 6395
end_IR5 = 6933
ip5_s = 6664.5684327563

MQXA_3L5_s = 6614.418432756301
MQXB_B2L5_s = 6623.268432756301
MQXB_A2L5_s = 6629.768432756301
MQXA_1L5_s = 6638.418432756301

MQXA_1R5_s = 6690.7184327563
MQXB_A2R5_s = 6699.3684327563
MQXB_B2R5_s = 6705.8684327563
MQXA_3R5_s = 6714.7184327563

MQXA_L = 6.37
MQXB_L = 5.55


def plot_triplets(ax=None):
    ax.fill_between([MQXA_3L5_s - MQXA_L/2., MQXA_3L5_s + MQXA_L/2.], [-500,-500], [9000,9000], color="b", alpha=0.5)
    ax.fill_between([MQXB_B2L5_s - MQXA_L/2., MQXB_B2L5_s + MQXB_L/2.], [-500,-500], [9000,9000], color="b", alpha=0.5)
    ax.fill_between([MQXB_A2L5_s - MQXA_L/2., MQXB_A2L5_s + MQXB_L/2.], [-500,-500], [9000,9000], color="b", alpha=0.5)
    ax.fill_between([MQXA_1L5_s - MQXA_L/2., MQXA_1L5_s + MQXA_L/2.], [-500,-500], [9000,9000], color="b", alpha=0.5)
    ax.fill_between([MQXA_1R5_s - MQXA_L/2., MQXA_1R5_s + MQXA_L/2.], [-500,-500], [9000,9000], color="b", alpha=0.5)
    ax.fill_between([MQXB_A2R5_s - MQXA_L/2., MQXB_A2R5_s + MQXB_L/2.], [-500,-500], [9000,9000], color="b", alpha=0.5)
    ax.fill_between([MQXB_B2R5_s - MQXA_L/2., MQXB_B2R5_s + MQXB_L/2.], [-500,-500], [9000,9000], color="b", alpha=0.5)
    ax.fill_between([MQXA_3R5_s - MQXA_L/2., MQXA_3R5_s + MQXA_L/2.], [-500,-500], [9000,9000], color="b", alpha=0.5)


plt.close("all")

##beta function
fig1 = plt.figure(1, figsize=[16/2., 9./2.])
ax11 = fig1.add_subplot(121)
ax11.plot(twiss_b1["s"], twiss_b1["betx"], "k", label="$\\beta_x$")
ax11.plot(twiss_b1["s"], twiss_b1["bety"], "r", label="$\\beta_y$")
ax11.set_xlabel("s [m]")
ax11.set_ylabel("β [m]")
ax11.set_xlim(ip5_s - 60, ip5_s + 60)
ax11.set_ylim(-500, 9000)
ax11.legend()
ax11.set_title("Beam 1")

ax12 = fig1.add_subplot(122)
ax12.plot(twiss_b2["s"], twiss_b2["betx"], "k", label="$\\beta_x$")
ax12.plot(twiss_b2["s"], twiss_b2["bety"], "r", label="$\\beta_y$")
ax12.set_xlabel("s [m]")
ax12.set_ylabel("β [m]")
ax12.set_xlim(ip5_s - 60, ip5_s + 60)
ax12.set_ylim(-500, 9000)
ax12.legend()
ax12.set_title("Beam 2")

## closed orbit
fig2 = plt.figure(2, figsize=[16/2., 9./2.])
ax21 = fig2.add_subplot(121)
ax21.plot(twiss_b1["s"], twiss_b1["x"], "k", label="$x$")
ax21.plot(twiss_b1["s"], twiss_b1["y"], "r", label="$y$")
ax21.set_xlabel("s [m]")
ax21.set_ylabel("x,y [m]")
ax21.set_xlim(ip5_s - 60, ip5_s + 60)
ax21.set_ylim(-0.015, 0.015)
ax21.legend()
ax21.set_title("Beam 1")

ax22 = fig2.add_subplot(122)
ax22.plot(twiss_b2["s"], twiss_b2["x"], "k", label="$x$")
ax22.plot(twiss_b2["s"], twiss_b2["y"], "r", label="$y$")
ax22.set_xlabel("s [m]")
ax22.set_ylabel("x,y [m]")
ax22.set_xlim(ip5_s - 60, ip5_s + 60)
ax22.set_ylim(-0.015, 0.015)
ax22.legend()

fig3 = plt.figure(3, figsize=[16/2., 9./2.])
ax31 = fig3.add_subplot(121)
ax31.plot(twiss_b1["s"], twiss_b1["dx"], "k", label="$D_x$")
ax31.plot(twiss_b1["s"], twiss_b1["dy"], "r", label="$D_y$")
ax31.set_xlabel("s [m]")
ax31.set_ylabel("D [m]")
ax31.set_xlim(ip5_s - 60, ip5_s + 60)
ax31.set_ylim(-1.5, 1.5)
ax31.legend()

ax32 = fig3.add_subplot(122)
ax32.plot(twiss_b2["s"], twiss_b2["dx"], "k", label="$D_x$")
ax32.plot(twiss_b2["s"], twiss_b2["dy"], "r", label="$D_y$")
ax32.set_xlabel("s [m]")
ax32.set_ylabel("D [m]")
ax32.set_xlim(ip5_s - 60, ip5_s + 60)
ax32.set_ylim(-1.5, 1.5)
ax32.legend()

fig4 = plt.figure(4, figsize=[16/2., 9./2.])
ax41 = fig4.add_subplot(121)
ax41.plot(twiss_b1["s"], twiss_b1["mux"], "k", label="$\\mu_x$")
ax41.plot(twiss_b1["s"], twiss_b1["muy"], "r", label="$\\mu_y$")
ax41.set_xlabel("s [m]")
ax41.set_ylabel("$\\mu$ [$2\\pi$]")
ax41.set_xlim(ip5_s - 60, ip5_s + 60)
ax41.set_ylim(13.9, 15.4)
ax41.legend()

ax42 = fig4.add_subplot(122)
ax42.plot(twiss_b2["s"], twiss_b2["mux"], "k", label="$\\mu_x$")
ax42.plot(twiss_b2["s"], twiss_b2["muy"], "r", label="$\\mu_y$")
ax42.set_xlabel("s [m]")
ax42.set_ylabel("$\\mu$ [$2\\pi$]")
ax42.set_xlim(ip5_s - 60, ip5_s + 60)
ax42.set_ylim(13.9, 15.6)
ax42.legend()

sigx_b1 = np.sqrt(twiss_b1["betx"] * nemitt_x / (line_b1.particle_ref.beta0 * line_b1.particle_ref.gamma0) + ( twiss_b1["dx"] * std_delta ) ** 2 )
sigy_b1 = np.sqrt(twiss_b1["bety"] * nemitt_x / (line_b1.particle_ref.beta0 * line_b1.particle_ref.gamma0) + ( twiss_b1["dy"] * std_delta ) ** 2 )

sigx_b2 = np.sqrt(twiss_b2["betx"] * nemitt_x / (line_b2.particle_ref.beta0 * line_b2.particle_ref.gamma0) + ( twiss_b2["dx"] * std_delta ) ** 2 )
sigy_b2 = np.sqrt(twiss_b2["bety"] * nemitt_x / (line_b2.particle_ref.beta0 * line_b2.particle_ref.gamma0) + ( twiss_b2["dy"] * std_delta ) ** 2 )

sig_twiss_b1 = twiss_b1.get_betatron_sigmas(nemitt_x=nemitt_x, nemitt_y=nemitt_y)
sig_twiss_b2 = twiss_b2.get_betatron_sigmas(nemitt_x=nemitt_x, nemitt_y=nemitt_y)

fig5 = plt.figure(5, figsize=[16/2., 9./2.])
ax51 = fig5.add_subplot(121)
ax51.plot(twiss_b1["s"], sigx_b1, "k", label="$\\sigma_x$")
ax51.plot(twiss_b1["s"], sigy_b1, "r", label="$\\sigma_y$")
ax51.plot(sig_twiss_b1["s"], sig_twiss_b1["sigma_x"], "k--", label="bet. $\\sigma_x$")
ax51.plot(sig_twiss_b1["s"], sig_twiss_b1["sigma_y"], "r--", label="bet. $\\sigma_y$")
ax51.set_xlabel("s [m]")
ax51.set_ylabel("$\\sigma$ [m]")
ax51.set_xlim(ip5_s - 60, ip5_s + 60)
ax51.set_ylim(0., 0.002)
ax51.legend()

ax52 = fig5.add_subplot(122)
ax52.plot(twiss_b2["s"], sigx_b2, "k", label="$\\sigma_x$")
ax52.plot(twiss_b2["s"], sigy_b2, "r", label="$\\sigma_y$")
ax52.plot(sig_twiss_b2["s"], sig_twiss_b2["sigma_x"], "k--", label="bet. $\\sigma_x$")
ax52.plot(sig_twiss_b2["s"], sig_twiss_b2["sigma_y"], "r--", label="bet. $ \\sigma_y$")
ax52.set_xlabel("s [m]")
ax52.set_ylabel("$\\sigma$ [m]")
ax52.set_xlim(ip5_s - 60, ip5_s + 60)
ax52.set_ylim(0., 0.002)
ax52.legend()

fig6 = plt.figure(6, figsize=[16/2., 9./2.])
ax61 = fig6.add_subplot(121)
ax61.plot(twiss_b1["s"], twiss_b1["alfx"], "k", label="$\\alpha_x$")
ax61.plot(twiss_b1["s"], twiss_b1["alfy"], "r", label="$\\alpha_y$")
ax61.set_xlabel("s [m]")
ax61.set_ylabel("$\\alpha$")
ax61.set_xlim(ip5_s - 60, ip5_s + 60)
ax61.set_ylim(-400., 400.)
ax61.legend()

ax62 = fig6.add_subplot(122)
ax62.plot(twiss_b2["s"], twiss_b2["alfx"], "k", label="$\\alpha_x$")
ax62.plot(twiss_b2["s"], twiss_b2["alfy"], "r", label="$\\alpha_y$")
ax62.set_xlabel("s [m]")
ax62.set_ylabel("$\\alpha$")
ax62.set_xlim(ip5_s - 60, ip5_s + 60)
ax62.set_ylim(-400., 400.)
ax62.legend()

for ax in ax11, ax12, ax21, ax22, ax31, ax32, ax41, ax42, ax51, ax52, ax61, ax62:
    plot_triplets(ax=ax)

for ax in ax11, ax21, ax31, ax41, ax51, ax61:
    ax.set_title("Beam 1")
for ax in ax12, ax22, ax32, ax42, ax52, ax62:
    ax.set_title("Beam 2")

for fig in fig1, fig2, fig3, fig4, fig5, fig6:
    fig.tight_layout()

interp_dict = {
"x_b1" : scipy.interpolate.interp1d(twiss_b1['s'], twiss_b1["x"]),
"y_b1" : scipy.interpolate.interp1d(twiss_b1['s'], twiss_b1["y"]),
"sig_x_b1" : scipy.interpolate.interp1d(twiss_b1['s'], sigx_b1),
"sig_y_b1" : scipy.interpolate.interp1d(twiss_b1['s'], sigy_b1),
"x_b2" : scipy.interpolate.interp1d(twiss_b2['s'], twiss_b2["x"]),
"y_b2" : scipy.interpolate.interp1d(twiss_b2['s'], twiss_b2["y"]),
"sig_x_b2" : scipy.interpolate.interp1d(twiss_b2['s'], sigx_b2),
"sig_y_b2" : scipy.interpolate.interp1d(twiss_b2['s'], sigy_b2),
}
#t_offset = 2 * (s - s_ip) / c

pkl.dump(interp_dict, open("beam_params.pkl", "wb"))

plt.show()