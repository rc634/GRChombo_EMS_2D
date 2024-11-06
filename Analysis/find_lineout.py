# parallel_L.py
# James Widdicombe
# Last Updated 16/12/2019
# Calculate L2M and L2H

# Load the modules
import yt
import numpy as np
import time
from scipy import interpolate

# Timings
start_time = time.time()

# Enable Parallelism
yt.enable_parallelism()

# Change this
data_location = "../2D/Collision_Bubble/hdf5/Bubble2Dp_00*"  # Data file location
save_folder = "lineout_data/2D_lambda_0.5_eps_0.025_phi0_0.025_R1.0/"
# End change this

# Loading dataset
ts = yt.load(data_location)

# Find length of box
right_edge = ts[0].domain_right_edge
L_box = float(right_edge[0])

# Define an empty storage dictionary for collecting information
# in parallel through processing
storage = {}

for sto, i in ts.piter(storage=storage):
    # Begin and end point of lineout
    ray_start = [0., 0., 0.]
    ray_end = [L_box, 0., 0.]

    # Take the lineout and sort it
    ray = i.r[ray_start:ray_end]
    ray_sort = np.argsort(ray["t"])

    # Define sorted variable arrays and turn from YTArray to float
    t = ray["t"][ray_sort]
    phi = ray["phi"][ray_sort]
    rho = ray["rho"][ray_sort]
    lapse = ray["lapse"][ray_sort]
    t_list = [float(i) for i in t]
    phi_list = [float(i) for i in phi]
    rho_list = [float(i) for i in rho]
    lapse_list = [float(i) for i in lapse]
    
    # Interpolate to make arrays from all timesteps the same length. If
    # not done, np.savetxt throws error at the end
    num_points = 200
    phi_interp = interpolate.interp1d(t_list, phi_list, bounds_error=False, fill_value='extrapolate')
    rho_interp = interpolate.interp1d(t_list, rho_list, bounds_error=False, fill_value='extrapolate')
    lapse_interp = interpolate.interp1d(t_list, lapse_list, bounds_error=False, fill_value='extrapolate')

    phi_save = [float(phi_interp(i)) for i in np.linspace(0, 1, num_points)]
    rho_save = [float(rho_interp(i)) for i in np.linspace(0, 1, num_points)]
    lapse_save = [float(lapse_interp(i)) for i in np.linspace(0, 1, num_points)]
    # Save
    array = [i.current_time, phi_save, rho_save, lapse_save]

    sto.result = array
    sto.result_id = str(i)

if yt.is_root():
    timedata = []
    phidata = []
    rhodata = []
    lapsedata = []

    for L in sorted(storage.items()):
        timedata.append(L[1][0])
        phidata.append(L[1][1])
        rhodata.append(L[1][2])
        lapsedata.append(L[1][3])
        
    np.savetxt(save_folder + "times.txt", timedata)
    np.savetxt(save_folder + "grid.txt", np.linspace(0, L_box, num_points))
    np.savetxt(save_folder + "phi.txt", np.array(np.array(phidata).T))
    np.savetxt(save_folder + "rho.txt", np.array(np.array(rhodata).T))
    np.savetxt(save_folder + "lapse.txt", np.array(np.array(lapsedata)))
