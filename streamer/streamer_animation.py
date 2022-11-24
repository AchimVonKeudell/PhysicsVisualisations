#
# Propagation of a streamer
# Examples for plasma pyhsics lectures
# Elia Juengling
# (Original Lazarus Code) - Achim von Keudell
# Ruhr University Bochum, 2022
#
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, writers

# ---------------------------
# simulation settings 
# ---------------------------
derivation = 'Upwind'   
time_step = 'Euler'     

# ----------------------------------
# --Parameters of the simulation ---
# -----------------------------------
n0 = 1e14
znb = 1000  # number of grid points
pressure_0 = 1e5  # pressure
scaledmue = 1
scaled = 1e-1
scalerecombination = 0*1e-12
scaledionization = 1
scalee0 = 1
y_axis = 4e18
E_0 = -5e6
    
# --- simulation constants ---
gp = 1  # number of ghost points per side
zw = 0.01   # length of simulation area in m
dz = zw / (znb-1) # cell length
dt0 = 1e-12 # time step
time = 0
tend = 20e-9
dtplot = 5e-11

# --- natural constants ----
el_charge = 1.602176e-19  # elementary charge
me = 9.109383e-31  # electron mass
eps0 = 8.854187e-12
kb = 1.380649e-23  # Boltzmann constant

# --- gas / plasma parameters
temp_n = 300.  # temperature of neutral atoms
temp_0 = 1.1
temp_e = temp_0 * np.ones(znb)


# --- functions for complex parameters and coefficients ---
def flux_limiter(rf):
    return max(0, min(2 * rf, 0.33 + 0.66 * rf, 2))

def de_air(ef):
    if ef != 0:
        result = 1 / scaled * 4.3628e-3 * np.exp(0.22 * np.log(abs(ef)))
    else:
        result = 0
    return result


def mue_air(ef):
    if ef != 0:
        result = - 1 / scaledmue * 2.398*np.exp(-0.26*np.log(abs(ef)))
    else:
        result = 0
    return result


def vcollision_electron(tn, p):
    result = (p / (kb * tn)) * 5.3e-13
    return result

def vionization_electron_air(ef, tn):
    if abs(ef) > 0:
       result = abs(scaledionization * (1.1944e6 + 4.366e26 / abs(ef ** 3)) * np.exp(-2.73e7 * scalee0/abs(ef))) - 340.75
    else:
        result = 0
    return result

def vtherm(te):
    global el_charge, me
    result = np.sqrt(3 * te * el_charge/me)
    return result

def vcollisionel_ion(tn, p):
    result = (p / (kb * tn)) * 8.3e-16
    return result

# -------------------------------------
# Initial conditions 
# -------------------------------------
z = np.zeros(znb + 2 * gp)  # initialisation of z
n_i = np.zeros(znb + 2 * gp)
n_e = np.zeros(znb + 2 * gp)# initialisation of ni
for i in range(znb + 2 * gp):
    z[i] = i/(znb-1) * zw
    n_i[i] = n0 * np.exp(-((z[i] - 0.2 * zw) ** 2) / ((0.01 * 0.005) ** 2))
    n_e[i] = n_i[i]


# ---------------------------------------
# Calculation of potential and field
# ---------------------------------------
efeld = np.zeros(znb + 2 * gp)
def calculate_potential_maxwell():
    global efeld, n_e, n_i, dz
    efeld_right = np.zeros(znb + 2 * gp)
    efeld_left = np.zeros(znb + 2 * gp)

    for i in range(znb + gp, znb + 2 * gp):
        efeld_right[i] = E_0
        efeld_left[i] = efeld_right[i]
        #efeld_right[i] = potential_0 / zw
    for i in range(znb + 2 * gp -1 ,1,-1):
        efeld_left[i] = efeld_right[i] + dz * el_charge/eps0 * (n_e[i] - n_i[i])
        efeld_right[i-1] = efeld_left[i]
    efeld_left[1] = efeld_right[2]
    efeld_right[0] = efeld_left[1]
    efeld_left[0] = efeld_right[0]

    for i in range(znb + 2 * gp):
        efeld[i] = 0.5 * (efeld_right[i]+efeld_left[i])


# -------------------------------------
# Calculation of densities 
# -------------------------------------
dtadvection = 1e-13
dtrelaxation = 1e-13
dtdiffusion = 1e-13
dn_e = np.zeros(znb + 2 * gp)
dn_i = np.zeros(znb + 2 * gp)

def calculate_derivation(ne,ni):
    global dz, dn_i, dn_e, dtadvection, dtdiffusion, dtrelaxation
    Erightmax = 0
    nemax = 1e8
    demax = 0
    dn_e_advection = np.zeros(znb+ 2* gp)
    dn_e_diffusion = np.zeros(znb+ 2* gp)
    dn_e_ionization = np.zeros(znb+ 2* gp)
    dn_e_recombination = np.zeros(znb+ 2* gp)


    for i in range(gp, znb+gp):

        # electrons
        de_el = de_air(efeld[i])
        mue_el = mue_air(efeld[i])

        Eleft = 0.5 * (efeld[i-1] + efeld[i])
        Eright = 0.5 * (efeld[i] + efeld[i+1])
        mue_el_left = mue_air(Eleft)
        mue_el_right = mue_air(Eright)
        de_el_right = de_air(Eright)
        de_el_left = de_air(Eleft)

        # Check time constant
        if abs(mue_el * Eright) > Erightmax:
            Erightmax = abs(mue_el * Eright)
            dtadvection = 0.1 * (dz * 0.5) / Erightmax

        if abs(ne[i]-ni[i]) > nemax:
            nemax = abs(ne[i]-ni[i])
            dtrelaxation = 0.1 * eps0 / (nemax * el_charge)

        if de_el > demax:
            demax = de_el
            dtdiffusion = 0.1 * 0.5 * (dz ** 2) / demax


        # Advection term
        if derivation == 'FTCS':
            dn_e_advection[i] = 0.5 * (ne[i]+ne[i-1])/dz * mue_el_left * Eleft - 0.5 * (ne[i]+ne[i+1])/dz * mue_el_right * Eright

        elif derivation == 'Upwind':
            if Eleft * mue_el > 0:
                if ne[i] - ne[i-1] == 0:
                    fllimiter = 2
                else:
                    rleft = (ne[i-1]-ne[i-2])/(ne[i]-ne[i-1])
                    fllimiter = flux_limiter(rleft)
                dn_e_advection[i] += (ne[i-1] + 0.5 * fllimiter * (ne[i] - ne[i-1])) / dz * mue_el_left * Eleft
            elif Eleft * mue_el < 0:
                if ne[i-1] - ne[i] == 0:
                    fllimiter = 2
                else:
                    rleft = (ne[i]-ne[i+1])/(ne[i-1]-ne[i])
                    fllimiter = flux_limiter(rleft)
                dn_e_advection[i] += (ne[i] + 0.5 * fllimiter * (ne[i-1] - ne[i]))/ dz * mue_el_left * Eleft

            if Eright * mue_el > 0:
                if ne[i+1] - ne[i] == 0:
                    fllimiter = 2
                else:
                    rright = (ne[i]-ne[i-1])/(ne[i+1]-ne[i])
                    fllimiter = flux_limiter(rright)
                dn_e_advection[i] -= (ne[i] + 0.5 * fllimiter * (ne[i+1] - ne[i]))/ dz * mue_el_right * Eright
            elif Eright * mue_el < 0:
                if ne[i] - ne[i+1] == 0:
                    fllimiter = 2
                else:
                    rleft = (ne[i+1]-ne[i+2])/(ne[i]-ne[i+1])
                    fllimiter = flux_limiter(rleft)
                dn_e_advection[i] -= (ne[i+1] + 0.5 * fllimiter * (ne[i] - ne[i+1]))/ dz * mue_el_right * Eright

        # Diffusion term
        dn_e_diffusion[i] = 1/(dz ** 2) * (de_el_right * (ne[i+1] - ne[i]) + de_el_left * (ne[i-1] - ne[i]))

        # Ionization
        vionization_e = vionization_electron_air(efeld[i], temp_n) * abs(mue_el) * abs(efeld[i])
        dn_e_ionization[i] = vionization_e * ne[i]

        # recombination
        dn_e_recombination[i] = - scalerecombination * ne[i] * ni[i]

        # Electron total change
        dn_e[i] = dn_e_advection[i] + dn_e_ionization[i] + dn_e_diffusion[i] + dn_e_recombination[i]

        # ions
        dn_i[i] = vionization_e * ne[i] - scalerecombination * ne[i] * ni[i]

        


def propagation_euler():
    global dn_e, n_e, dt0, dn_i, n_i
    calculate_derivation(n_e, n_i)

    dt = dt0

    if dtadvection < dt:
        dt = dtadvection
    if dtdiffusion < dt:
        dt = dtdiffusion

    for i in range(gp, znb + gp):
        n_e[i] += dn_e[i] * dt
        n_i[i] += dn_i[i] * dt
        if n_e[i] < 0:
            n_e[i] = 0
        if n_i[i] < 0:
            n_i[i] = 0

    for i in range(gp):
        n_e[i] = n_e[gp]
        n_i[i] = n_i[gp]

    for i in range(znb+gp, znb + 2 * gp):
        n_e[i] = n_e[znb + gp - 1]
        n_i[i] = n_i[znb + gp - 1]


def run():
    calculate_potential_maxwell()
    propagation_euler()


# -------------------------------------
# Define Plots
# -------------------------------------
fig, ax = plt.subplots(1, 1,figsize=(8,4.5))

ax.set(xlim=(0, (znb+2 * gp)*dz/1e-3), ylim=(0, y_axis),xlabel = 'z (mm)',ylabel='density (m^-3)')
axp = ax.twinx()
axp.set(xlim=(0, (znb+2 * gp)*dz/1e-3), ylim=(0  , -6e6),ylabel='E field (V/m)')

x = np.linspace(0, znb + 2 * gp, znb + 2 * gp)
line1 = ax.plot(x*dz/1e-3, n_e, color='r', lw=1)[0]
line2 = ax.plot(x*dz/1e-3, n_i, color='b', lw=1)[0]
line4 = axp.plot(x*dz/1e-3, efeld, color='g', lw=1)[0]

line1.set_label('$n_e$')
line2.set_label('$n_i$')
line4.set_label('E')

ax.legend(loc=2)
axp.legend(loc=1)


# calculate velocity
result_space = []
result_time = []
average_span = 5
velocity_left = []
velocity_right = []
v = 0


def calculate_velocity(space,time):
    k = len(space)
    if len(time) != k:
        print('laengen stimmen nicht ueberein')
        exit()
    v = (space[k-1] - space[k-1-average_span]) / (time[k-1] - time[k-1]- average_span) * dz/dt0
    return v


def animate(i):
    global time, dt0, v

    timelocal = 0
    
    while timelocal<dtplot:
        dt = dt0
        timelocal += dt
        if dtadvection < dt:
            dt = dtadvection

        if dtdiffusion < dt:
            dt = dtdiffusion

        run()
        time += dt
        store = np.zeros(int(znb / 2 + gp))

        for i in range(int(znb / 2 + gp)):
            store[i] = n_e[i]
        result_1 = np.where(store == np.max(store))

        if result_1[0][0] < znb /2 - 20:
            result_space.append(result_1[0][0])
            result_time.append(time)
        if len(result_space) > average_span +1:
            v = calculate_velocity(result_space, result_time)
            velocity_left.append(v)

    print('time: '+ "{:.2f}".format(time/1e-9)+' ns of '+"{:.0f}".format(tend/1e-9)+' ns ')
    title = "time = " + "{:.2f}".format(time/1e-9) + 'ns  '
    fig.suptitle(title, y = 1)
    line1.set_ydata(n_e)
    line2.set_ydata(n_i)
    line4.set_ydata(efeld)

anim = FuncAnimation(fig, animate, interval=1, frames=int(tend /dtplot), repeat=False)
anim.save('streamer_air_pot30e3.gif', writer='pillow')

plt.draw()
plt.tight_layout()
plt.show()


