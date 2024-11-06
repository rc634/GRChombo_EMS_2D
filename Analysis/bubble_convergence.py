from scipy.optimize import curve_fit
from scipy import interpolate
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.size'] = '18'
plt.rcParams['lines.linewidth'] = '2.5'
figure_folder = 'Figures/'


def extract_constraints(constraints_file):
    time = np.loadtxt(constraints_file, skiprows = 1, usecols = 0)
    ham = np.loadtxt(constraints_file, skiprows = 1, usecols = 1)
    mom = np.loadtxt(constraints_file, skiprows = 1, usecols = 2)

    return time, abs(ham), mom

def extract_M_ADM(M_file):
    time_constraint = np.loadtxt(M_file, skiprows = 1, usecols = 0)
    mass = np.loadtxt(M_file, skiprows = 1, usecols = 1)

    return time_constraint, mass


def main():
    time_N256, ham_N256, mom_N256 = extract_constraints("../2D/Bubble_Level4/data/constraint_norms_bubble_N256.dat")
    time_N512, ham_N512, mom_N512 = extract_constraints("../2D/Bubble_Level4/data/constraint_norms_bubble_N512.dat")
    time_N1024, ham_N1024, mom_N1024 = extract_constraints("../2D/Bubble_Level4/data/constraint_norms_bubble_N1024.dat")

    factor_3D = 4. / np.pi # accounts for difference in box shape between 2D/3D
    fig, ax = plt.subplots(2, 1, figsize=(10, 12))
    fig.suptitle(r"Bubble convergence test")
    ax[0].plot(time_N256, ham_N256, label="N256", linestyle='-')
    ax[0].plot(time_N512, ham_N512, label="N512", linestyle='-')
    ax[0].plot(time_N1024, ham_N1024, label="N1024", linestyle='-')
    ax[0].set_xlabel(r"$t_\mathrm{code}$")
    ax[0].set_ylabel(r"$\int_{\mathrm{Box}}\mathcal{H}$")
    ax[0].set_yscale('log')
    # ax[0].set_xlim([0, 1.])
    ax[0].set_ylim([1.e-8, 1.e-1])
    ax[0].legend(loc='lower left')
    ax[0].set_title(r"$\mathcal{H}$")
    
    ax[1].plot(time_N256, mom_N256, label="N256", linestyle='-')
    ax[1].plot(time_N512, mom_N512, label="N512", linestyle='-')
    ax[1].plot(time_N1024, mom_N1024, label="N1024", linestyle='-')
    ax[1].set_xlabel(r"$t_\mathrm{code}$")
    ax[1].set_ylabel(r"$\int_{\mathrm{Box}}\mathcal{M}$")
    ax[1].set_yscale('log')
    # ax[1].set_xlim([0, 1.])
    ax[1].set_ylim([1.e-2, 1.e0])
    ax[1].legend(loc='lower right')
    ax[1].set_title(r"$\mathcal{M}$")

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(figure_folder + "bubble_constraints.png")
    plt.close()

if __name__ == "__main__":
    main()
