import matplotlib.pyplot as plt
import matplotlib
import shapely.geometry as sg
import shapely.ops as so
import numpy as np
import pickle

def plot_rectcircle(halfwidth, halfheight, radius):
    
    p1 = sg.Point(0,0).buffer(radius)
    aperture = so.clip_by_rect(p1,-halfwidth,-halfheight, halfwidth,halfheight)
    coords = np.array(list(aperture.exterior.coords))
    return coords
    # x0 = -halfwidth
    # y0 = -halfheight
    # rect = matplotlib.patches.Rectangle((x0,y0), 2*halfwidth, 2*halfheight, linewidth=1, edgecolor='k', facecolor='none')
    # circ = matplotlib.patches.Circle((0,0), radius, linewidth=1, edgecolor='k', facecolor='none')
    # return rect, circ, aperture

n_sigma=5
beams = pickle.load(open("beam_params.pkl","rb"))
s = 6690
ip1_s = 19994.1624
ip5_s = 6664.5684327563

def plot_in_triplet(s, beams=beams):
    x_b1 = beams["x_b1"](s)
    y_b1 = beams["y_b1"](s)
    sig_x_b1 = beams["sig_x_b1"](s)
    sig_y_b1 = beams["sig_y_b1"](s)
    x_b2 = beams["x_b2"](s)
    y_b2 = beams["y_b2"](s)
    sig_x_b2 = beams["sig_x_b2"](s)
    sig_y_b2 = beams["sig_y_b2"](s)
    
    theta = np.linspace(0,2*np.pi,1000)
    circle_x = np.cos(theta)
    circle_y = np.sin(theta)
    
    bs001 = plot_rectcircle(halfwidth=0.02202, halfheight=0.01714, radius=0.02202)
    bs010 = plot_rectcircle(halfwidth=0.0192, halfheight=0.0241, radius=0.0241)
    bs012 = plot_rectcircle(halfwidth=0.0241, halfheight=0.029, radius=0.029)
    bs003 = plot_rectcircle(halfwidth=0.0241, halfheight=0.0192, radius=0.0241)
    bs005 = plot_rectcircle(halfwidth=0.029, halfheight=0.0241, radius=0.029)
    if abs(s - ip5_s) < abs(ip1_s):
        if abs(s - ip5_s) < 31.40:
            beamscreen = bs010 
        else:
            beamscreen = bs012 
    else:
        if abs(s - ip1_s) < 31.40:
            beamscreen = bs003
        else:
            beamscreen = bs005 

    
    
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(beamscreen[:,0]*1000, beamscreen[:,1]*1000,'k-')
    ax1.set_aspect("equal")
    ax1.set_xlabel("x [mm]")
    ax1.set_ylabel("y [mm]")
    ax1.plot([-30,+30], [-30,+30], 'k--')
    ax1.plot([-30,+30], [+30,-30], 'k--')
    ax1.plot((circle_x*sig_x_b1*n_sigma-x_b1)*1000, (circle_y*sig_y_b1*n_sigma-y_b1)*1000, 'b-')
    ax1.plot((circle_x*sig_x_b2*n_sigma-x_b2)*1000, (circle_y*sig_y_b2*n_sigma-y_b2)*1000, 'r-')
    ax1.set_title(f"s = {s:.1f}m")
    ax1.set_xlim(-30,30)
    ax1.set_ylim(-30,30)

for s in np.linspace(ip5_s + 23, ip5_s + 55, 10):
    plot_in_triplet(s, beams=beams)

plt.show()