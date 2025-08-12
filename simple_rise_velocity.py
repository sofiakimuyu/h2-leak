import math

# Constants
g = 9.80665  # gravity (m/s^2)
rho_air = 1.225  # air density (kg/m^3) at STP
rho_h2 = 0.0899  # hydrogen density (kg/m^3) at STP
V = 0.01  # parcel volume in m^3 (adjust as needed)
C_d = 0.47  # drag coefficient for sphere
dt = 0.0001  # time step (s)
t_max = 5  # total simulation time (s)
print_interval = .25  # seconds between velocity outputs

# Derived values
radius = (3*V/(4*math.pi))**(1/3)
A = math.pi * radius**2  # cross-sectional area (m^2)
m_parcel = rho_h2 * V  # parcel mass (kg)

# Initial conditions
v = 0.0
t = 0.0
next_print = 0.0

print(f"{'Time (s)':>10} | {'Velocity (m/s)':>15}")
print("-" * 28)

while t <= t_max:
    # Buoyant force
    F_buoy = (rho_air - rho_h2) * g * V

    # Drag force:
    F_drag = 0.5 * C_d * rho_air * A * v * abs(v)

    # Net force and acceleration
    F_net = F_buoy - F_drag
    a = F_net / m_parcel

    # Update velocity
    v += a * dt
    t += dt

    # Print
    if t >= next_print:
        print(f"{t:10.1f} | {v:15.5f}")
        next_print += print_interval
