#!/usr/bin/env python3
"""orbital_mech - Orbital mechanics calculator (Kepler, Hohmann transfer, vis-viva)."""
import sys, math

G = 6.674e-11  # gravitational constant
M_EARTH = 5.972e24
R_EARTH = 6.371e6

def orbital_velocity(mass, radius, semi_major=None):
    if semi_major is None:
        semi_major = radius
    return math.sqrt(G * mass * (2/radius - 1/semi_major))

def orbital_period(mass, semi_major):
    return 2 * math.pi * math.sqrt(semi_major**3 / (G * mass))

def escape_velocity(mass, radius):
    return math.sqrt(2 * G * mass / radius)

def hohmann_transfer(mass, r1, r2):
    a_transfer = (r1 + r2) / 2
    v1 = orbital_velocity(mass, r1)
    v2 = orbital_velocity(mass, r2)
    vt1 = orbital_velocity(mass, r1, a_transfer)
    vt2 = orbital_velocity(mass, r2, a_transfer)
    dv1 = abs(vt1 - v1)
    dv2 = abs(v2 - vt2)
    transfer_time = orbital_period(mass, a_transfer) / 2
    return dv1, dv2, dv1 + dv2, transfer_time

def kepler_third_law(period=None, semi_major=None, mass=M_EARTH):
    if period is not None:
        return (G * mass * period**2 / (4 * math.pi**2))**(1/3)
    elif semi_major is not None:
        return 2 * math.pi * math.sqrt(semi_major**3 / (G * mass))
    raise ValueError("Provide period or semi_major")

def test():
    # LEO ~400km altitude
    r_leo = R_EARTH + 400e3
    v_leo = orbital_velocity(M_EARTH, r_leo)
    assert abs(v_leo - 7672) < 50  # ~7.67 km/s
    T_leo = orbital_period(M_EARTH, r_leo)
    assert abs(T_leo - 5560) < 100  # ~92.7 min
    v_esc = escape_velocity(M_EARTH, R_EARTH)
    assert abs(v_esc - 11186) < 50  # ~11.2 km/s
    # GEO
    r_geo = kepler_third_law(period=86400, mass=M_EARTH)
    assert abs(r_geo - 42164e3) < 100e3  # ~42,164 km
    # Hohmann: LEO to GEO
    dv1, dv2, total_dv, t_transfer = hohmann_transfer(M_EARTH, r_leo, r_geo)
    assert total_dv > 3000  # ~3.9 km/s total delta-v
    assert t_transfer > 10000  # several hours
    print("OK: orbital_mech")

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "test":
        test()
    else:
        print("Usage: orbital_mech.py test")
