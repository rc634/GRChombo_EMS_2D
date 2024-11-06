import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from scipy import interpolate

# Plot params and read/save folders
plt.rcParams['font.size'] = '14'
plt.rcParams['lines.linewidth'] = '2.5'
figure_folder = 'Figures/'

# Change this
radius = "R1.0"
dt_refine_factor = 4
false_vacuum = 0.0236412
barrier = 0.00167417
true_vacuum = -0.026167
# End change this

figure_name = "contourplot2D_" + radius + ".png"
fig_suptitle = r"$\frac{R}{R_\mathrm{crit}}$ = " + radius[1:]
read_folder = "lineout_data/2D_lambda_0.5_eps_0.025_phi0_0.025_" + radius + "/"

# Data loading
times = np.loadtxt(read_folder + "times.txt")
grid = np.loadtxt(read_folder + "grid.txt")
phi_base = np.loadtxt(read_folder + "phi.txt", unpack=True)
lapse_base = np.loadtxt(read_folder + "lapse.txt", usecols=-1)

# Grid setup
# Interpolate to be able to refine time resolution
interp_phi = interpolate.RectBivariateSpline(times, grid, phi_base)
interp_lapse = interpolate.interp1d(times, lapse_base, bounds_error=False, fill_value='extrapolate')

dtime, dx = (times[1] - times[0]) / dt_refine_factor, grid[1] - grid[0]

times_fine = np.arange(times[0], times[-1] + dtime, dtime)
grid_fine = np.arange(grid[0], grid[-1] + dx, dx)

y, x = np.mgrid[slice(times[0], times[-1] + dtime, dtime), slice(grid[0], grid[-1] + dx, dx)]

phi = interp_phi(times_fine, grid_fine)
lapse = interp_lapse(times_fine)

# pick the desired colormap, sensible levels, and define a normalization
# instance which takes data values and translates those into levels.
levels = MaxNLocator(nbins=15).tick_values(phi.min(), phi.max())
cmap = plt.get_cmap('PiYG')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

# Field values over time at collision point
phi_collision_point = phi[:,-1]

# Plot
fig, ax = plt.subplots(3, 1, figsize=(10, 15))
fig.suptitle(fig_suptitle)
cf = ax[0].pcolormesh(x, y, phi, cmap=cmap, norm=norm, label='d')
ax[0].set_xlabel(r"$x_\mathrm{code}$")
ax[0].set_ylabel(r"$t_\mathrm{code}$")
ax[0].set_xlim(left=.0)
ax[0].ticklabel_format(axis="both", style="sci", scilimits=(0, 0))
fig.colorbar(cf, ax=ax[0], label=r"$\phi$")
ax[0].set_title(r"$\phi$ along x-axis over time")

ax[1].plot(times_fine, phi_collision_point, color='tab:olive')
ax[1].axhline(y=false_vacuum, label=r"False vacuum", linestyle=":", color='black')
ax[1].axhline(y=barrier, label=r"Saddle point", linestyle="dashdot", color='green')
ax[1].axhline(y=true_vacuum, label=r"True vacuum", linestyle="--", color='navy')
ax[1].set_xlabel(r"$t_\mathrm{code}$")
ax[1].set_ylabel(r"$\phi_\mathrm{collision point}$")
ax[1].set_xlim(left=.0)
ax[1].ticklabel_format(axis="both", style="sci", scilimits=(0, 0))
ax[1].legend(loc='lower left')
ax[1].set_title(r"$\phi$ at collision point over time")

ax[2].plot(times, lapse_base, color='tab:olive')
ax[2].axhline(y=1., label=r"Lapse = 1.0", linestyle=":", color='black')
ax[2].set_xlabel(r"$t_\mathrm{code}$")
ax[2].set_ylabel(r"$\alpha$")
ax[2].ticklabel_format(axis="both", style="sci", scilimits=(0, 0))
ax[2].legend(loc='lower left')
ax[2].set_title(r"$\alpha$ at collision point over time")

fig.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig(figure_folder + figure_name)
plt.close()

