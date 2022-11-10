import xtrack as xt
import numpy as np

def collimator_setup(line, context, number_of_sigmas=5.7, reference_norm_emit=3.5e-6, verbose=False):

    print("Setting up the three TCP collimators in IR7...")

    line_copy = line.copy()
    tracker = xt.Tracker(_context=context, line=line_copy)
    twiss = tracker.twiss()
    tcp1 = 'tcp.d6l7.b1'
    tcp2 = 'tcp.c6l7.b1'
    tcp3 = 'tcp.b6l7.b1'

    nc = number_of_sigmas
    epsn = reference_norm_emit
    gamma0 = line.particle_ref.gamma0
    beta0 = line.particle_ref.beta0
    epsg = epsn/beta0/gamma0


    sbx = {}
    sby = {}
    xco = {}
    yco = {}
    deg3 = 126.91 
    rot3 = deg3
    rad3 = deg3*np.pi/180.
    for i in range(len(twiss["name"])):
        elname = twiss["name"][i].split(':')[0]
        if elname in [tcp1, tcp2, tcp3]:
            sbx[elname] = np.sqrt(twiss["betx"][i]*epsg)
            sby[elname] = np.sqrt(twiss["bety"][i]*epsg)
            xco[elname] = twiss["x"][i]
            yco[elname] = twiss["y"][i]
            if verbose:
                print(twiss["name"][i])
                print(twiss["betx"][i])
                print(twiss["bety"][i])
                print(twiss["x"][i])
                print(twiss["y"][i])
                print()

    tcp_d6l7 = xt.beam_elements.LimitRect(min_y = -nc*sby[tcp1] + yco[tcp1], max_y = nc*sby[tcp1] + yco[tcp1], min_x=-1, max_x=1)
    tcp_c6l7 = xt.beam_elements.LimitRect(min_x = -nc*sbx[tcp2] + xco[tcp2], max_x = nc*sbx[tcp2] + xco[tcp2], min_y=-1, max_y=1)
    #sb3 = ((sbx[tcp3]*np.cos(rad3))**2 + (sby[tcp3]*np.sin(rad3))**2)**0.5
    phi3 = np.arctan( - np.tan(np.pi/2. - rad3) * sbx[tcp3] / sby[tcp3] )
    sb3 = sby[tcp3] * np.sin(rad3) / np.cos(phi3)
    # #sb3 = np.sqrt(1./((np.cos(rad3)/sbx[tcp3])**2 + (np.sin(rad3)/sby[tcp3])**2))
    # #sb3 = ((sbx[tcp3]*np.cos(rad3))**2 + (sby[tcp3]*np.sin(rad3))**2)**0.5
    # #co3 = (xco[tcp3]**2 + yco[tcp3]**2)**0.5
    co3=0
    shift1 = xt.beam_elements.XYShift(dx=xco[tcp3],dy=yco[tcp3])
    shift2 = xt.beam_elements.XYShift(dx=-xco[tcp3],dy=-yco[tcp3])
    rot1 = xt.beam_elements.SRotation(angle=+rot3)
    rot2 = xt.beam_elements.SRotation(angle=-rot3)
    tcp_b6l7 = xt.beam_elements.LimitRect(min_x = -nc*sb3, max_x = nc*sb3, min_y=-1, max_y=1)

    #collimators = ['tcp.d6l7.b1','tcp.c6l7.b1','tcp.b6l7.b1']
    for ii,elname in enumerate(line.element_names):
        if elname == tcp1:
            if verbose:
                print(ii)
            line.insert_element(index=ii, element=tcp_d6l7, name='tcp.d6l7.b1.coll')
            if verbose:
                print(line.element_names[ii-5:ii+5])
            break

    for ii,elname in enumerate(line.element_names):
        if elname == tcp2:
            if verbose:
                print(ii)
            line.insert_element(index=ii, element=tcp_c6l7, name='tcp.c6l7.b1.coll')
            if verbose:
                print(line.element_names[ii-5:ii+5])
            break

    for ii,elname in enumerate(line.element_names):
        if elname == tcp3:
            itcp3 = ii

    line.insert_element(index=itcp3, element=tcp_b6l7, name='tcp.b6l7.b1.coll')
    line.insert_element(itcp3+2,shift2,'tcp.b6l7.b1.shift2')
    line.insert_element(itcp3+2,rot2,'tcp.b6l7.b1.rot2')
    line.insert_element(itcp3,rot1,'tcp.b6l7.b1.rot1')
    line.insert_element(itcp3,shift1,'tcp.b6l7.b1.shift1')

    print("Collimator setup done.")
    return
