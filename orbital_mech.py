#!/usr/bin/env python3
"""Orbital mechanics — Kepler's laws, orbital elements, vis-viva."""
import sys, math

G = 6.674e-11  # gravitational constant
M_EARTH = 5.972e24
R_EARTH = 6.371e6

def orbital_velocity(r, a, mu=None):
    mu = mu or G * M_EARTH
    return math.sqrt(mu * (2/r - 1/a))

def escape_velocity(r, mu=None):
    mu = mu or G * M_EARTH
    return math.sqrt(2 * mu / r)

def orbital_period(a, mu=None):
    mu = mu or G * M_EARTH
    return 2 * math.pi * math.sqrt(a**3 / mu)

def hohmann_transfer(r1, r2, mu=None):
    mu = mu or G * M_EARTH
    a_transfer = (r1 + r2) / 2
    v1 = math.sqrt(mu / r1)
    v1_transfer = math.sqrt(mu * (2/r1 - 1/a_transfer))
    v2 = math.sqrt(mu / r2)
    v2_transfer = math.sqrt(mu * (2/r2 - 1/a_transfer))
    dv1 = abs(v1_transfer - v1); dv2 = abs(v2 - v2_transfer)
    transfer_time = math.pi * math.sqrt(a_transfer**3 / mu)
    return {"dv1": dv1, "dv2": dv2, "total_dv": dv1 + dv2, "transfer_time": transfer_time}

def kepler_third(T, mu=None):
    mu = mu or G * M_EARTH
    return (mu * T**2 / (4 * math.pi**2)) ** (1/3)

def main():
    if len(sys.argv) < 2: print("Usage: orbital_mech.py <demo|test>"); return
    if sys.argv[1] == "test":
        # LEO circular orbit (~400km)
        r_leo = R_EARTH + 400e3
        v = orbital_velocity(r_leo, r_leo)
        assert 7000 < v < 8000  # ~7.67 km/s
        T = orbital_period(r_leo)
        assert 5000 < T < 6000  # ~92 min = 5520s
        v_esc = escape_velocity(r_leo)
        assert v_esc > v  # escape > orbital
        assert abs(v_esc / v - math.sqrt(2)) < 0.01
        # GEO
        T_geo = 86164  # sidereal day in seconds
        r_geo = kepler_third(T_geo)
        assert abs(r_geo / 1e6 - 42.164) < 0.5  # ~42,164 km
        # Hohmann LEO to GEO
        h = hohmann_transfer(r_leo, r_geo)
        assert h["total_dv"] > 3000  # ~3.9 km/s
        assert h["transfer_time"] > 10000  # ~5 hours
        print("All tests passed!")
    else:
        r = R_EARTH + 400e3
        print(f"LEO v={orbital_velocity(r,r)/1000:.2f} km/s, T={orbital_period(r)/60:.1f} min")
        h = hohmann_transfer(r, kepler_third(86164))
        print(f"Hohmann to GEO: dv={h['total_dv']/1000:.2f} km/s, time={h['transfer_time']/3600:.1f} h")

if __name__ == "__main__": main()
