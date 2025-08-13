import math

# -----------------------------
# Physical constants
# -----------------------------
g = 9.80665                  # m/s^2
R = 8.31446261815324         # J/(mol K)
M_air = 0.0289647            # kg/mol
M_H2  = 0.00201588           # kg/mol
mu_air = 1.81e-5              # Pa.s
D_H2_air = 6.1e-5             # m^2/s binary diffusion coefficient

# Model coefficients
alpha = 0.10                  # entrainment coefficient
C_added = 0.5                 # added mass coeff
C_d = 0.47                     # drag coeff

def rho_from_ideal(P, T, M):
    return P * M / (R * T)

def drag_force(v, rho_air, radius):
    A = math.pi * radius**2
    return 0.5 * C_d * rho_air * A * v * abs(v)

def step_rk4(state, dt, params):
    def deriv(s):
        z, v, V, m_air, m_h2 = s
        rho_air = params['rho_air']

        radius = (3.0*V/(4.0*math.pi))**(1/3)
        A_surface = 4 * math.pi * radius**2

        # Diffusion
        rho_H2_parcel = m_h2 / V
        rho_H2_air = 0.0  # ambient hydrogen
        dm_h2_dt = -D_H2_air * A_surface * (rho_H2_parcel - rho_H2_air) / radius
        dm_air_from_diffusion = -dm_h2_dt  # same volume replaced by air

        # Entrainment
        A_cross = math.pi * radius**2
        dVdt = alpha * A_cross * v if v > 0 else 0.0
        dm_air_from_entrain = rho_air * dVdt

        # Parcel mass & density
        m_parcel = m_h2 + m_air
        rho_p = m_parcel / V

        # Buoyancy
        F_b = (rho_air - rho_p) * g * V

        #Added mass
        m_added = C_added * rho_air * V
        m_eff = m_parcel + m_added

        # Drag
        F_d = drag_force(v, rho_air, radius)

        # Motion equations
        dvdt = (F_b - F_d) / m_eff
        dzdt = v
        dVdt_total = dVdt  # only entrainment changes volume
        dm_air_dt = dm_air_from_entrain + dm_air_from_diffusion

        return (dzdt, dvdt, dVdt_total, dm_air_dt, dm_h2_dt)

    k1 = deriv(state)
    s2 = tuple(state[i] + 0.5*dt*k1[i] for i in range(5))
    k2 = deriv(s2)
    s3 = tuple(state[i] + 0.5*dt*k2[i] for i in range(5))
    k3 = deriv(s3)
    s4 = tuple(state[i] + dt*k3[i] for i in range(5))
    k4 = deriv(s4)
    return tuple(state[i] + dt*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i])/6.0 for i in range(5))

def simulate(T=293.15, P=101325.0, V0=0.01, dt=0.0001, t_max=5.0):
    rho_air = rho_from_ideal(P, T, M_air)
    rho_h2  = rho_from_ideal(P, T, M_H2)
    m_h2_0 = rho_h2 * V0

    # state: z, v, V, m_air, m_h2
    state = (0.0, 0.0, V0, 0.0, m_h2_0)
    params = {'rho_air': rho_air}
    t = 0.0
    history = [(t, *state, rho_h2)]

    while t < t_max:
        state = step_rk4(state, dt, params)
        t += dt
        z, v, V, m_air, m_h2 = state
        rho_p = (m_h2 + m_air) / V
        history.append((t, z, v, V, m_air, m_h2, rho_p))
        if rho_p >= rho_air:
            break
    return history

#Run & print
hist = simulate()

print("t[s]    z[m]    v[m/s] ")
for row in hist[:: max(1, len(hist)//20)]:
    t,z,v,V,m_air,m_h2,rho_p = row
    print(f"{t:6.3f}  {z:7.3f}  {v:7.3f} ")

