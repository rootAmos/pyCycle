# simple_turbojet.py — Output guide

The script runs one **design** point and two **off-design** points (OD0, OD1), then prints a report for each. Units are **English** (e.g. lbf, lbm/s, ft, °R) unless noted.

---

## 1. Cycle points (what DESIGN / OD0 / OD1 are)

| Point   | Meaning |
|--------|---------|
| **DESIGN** | Design point: sea-level static (alt=0, Mach≈0). You set **thrust target** (11800 lbf), **T4** (2370 °R), **compressor PR** (13.5), efficiencies. The solver finds mass flow **W**, FAR, and turbine PR so that thrust and shaft power match. |
| **OD0** | Off-design: same flight (alt=0, Mach≈0) but **lower thrust target** (11000 lbf). Nozzle area is fixed from design. Solver finds W, FAR, Nmech. |
| **OD1** | Off-design: **different flight** (Mach 0.2, 5000 ft) and **thrust target** 8000 lbf. Again nozzle area from design; solver finds W, FAR, Nmech. |

---

## 2. PERFORMANCE CHARACTERISTICS (one row per point)

| Column | Meaning | Typical values (from your run) |
|--------|--------|----------------------------------|
| **Mach** | Flight Mach number | 0 (static), 0.2 (OD1) |
| **Alt** | Altitude | ft (0, 5000) |
| **W** | Engine inlet mass flow | lbm/s (e.g. 147.3 design, 119.1 OD1) |
| **Fn** | Net thrust (Fg − Fram) | lbf (target thrust) |
| **Fg** | Gross thrust (nozzle) | lbf |
| **Fram** | Ram drag (inlet momentum) | lbf (0 at Mach 0, ~812 at Mach 0.2) |
| **OPR** | Overall pressure ratio (Pt3/Pt2) | e.g. 13.5 design, 12.2 OD1 |
| **TSFC** | Thrust-specific fuel consumption | (lbm/h)/lbf (e.g. 0.798 design) |

- At **Mach 0**, Fn = Fg and Fram = 0. At **Mach 0.2**, Fg > Fn because ram drag is subtracted.

---

## 3. FLOW STATIONS

One row per **station** along the engine: freestream → inlet exit → compressor exit → burner exit → turbine exit → nozzle exit.

| Column | Meaning | Units |
|--------|--------|--------|
| **tot:P** | Total (stagnation) pressure | psi |
| **tot:T** | Total temperature | °R |
| **tot:h** | Total enthalpy | Btu/lbm |
| **tot:S** | Total entropy | — |
| **stat:P** | Static pressure | psi |
| **stat:W** | Mass flow (same at each station in 1-D) | lbm/s |
| **stat:MN** | Mach number | — |
| **stat:V** | Velocity | ft/s |
| **stat:area** | Cross-sectional area | in² |

- **fc.Fl_O** = freestream (upstream of inlet).  
- **nozz.Fl_O** = nozzle exit (expansion to ambient **stat:P** = 14.696 psi at sea level).  
- Total pressure **drops** through compressor (loss), burner (dPqP), turbine, nozzle. Total **temperature** rises in compressor and burner, drops in turbine and nozzle.

---

## 4. COMPRESSOR PROPERTIES

| Symbol | Meaning |
|--------|--------|
| **Wc** | Corrected mass flow (map coordinate) |
| **Pr** | Pressure ratio (OPR for this compressor) |
| **eta_a, eta_p** | Adiabatic / polytropic efficiency |
| **Nc** | Corrected speed | rpm |
| **pwr** | Power absorbed (negative = into compressor) | hp |
| **RlineMap, NcMap** | Map coordinates (operating point on map) |
| **PRmap, WcMap, effMap** | Map grid / values |
| **alphaMap** | Map interpolation index (e.g. IGV angle) |

---

## 5. BURNER PROPERTIES

| Symbol | Meaning |
|--------|--------|
| **dPqP** | Total pressure loss (fraction) |
| **TtOut** | Burner exit total temperature (T4) | °R |
| **Wfuel** | Fuel flow | lbm/s |
| **FAR** | Fuel–air ratio |

---

## 6. TURBINE PROPERTIES

| Symbol | Meaning |
|--------|--------|
| **Wp** | Corrected flow (turbine map) |
| **PR** | Turbine pressure ratio |
| **eff_a, eff_p** | Adiabatic / polytropic efficiency |
| **Np** | Corrected speed |
| **pwr** | Power output (positive) | hp |

Shaft balance: compressor **pwr** (magnitude) ≈ turbine **pwr**.

---

## 7. NOZZLE PROPERTIES

| Symbol | Meaning |
|--------|--------|
| **PR** | Nozzle total-to-static pressure ratio |
| **Cv** | Velocity coefficient (loss) |
| **Ath** | Throat area | in² (fixed from design in off-design) |
| **MNth** | Throat Mach (≈ 1 when choked) |
| **MNout** | Exit Mach |
| **V** | Exit velocity | ft/s |
| **Fg** | Gross thrust | lbf |

---

## 8. SHAFT PROPERTIES

| Symbol | Meaning |
|--------|--------|
| **Nmech** | Shaft speed | rpm |
| **trqin / trqout** | Torque in (compressor) / out (turbine) | ft·lbf |
| **pwrin / pwrout** | Power in / out | hp (equal magnitude at convergence) |

---

## Quick comparison (your run)

| Quantity | DESIGN | OD0 | OD1 |
|----------|--------|-----|-----|
| Flight | 0 ft, M≈0 | 0 ft, M≈0 | 5000 ft, M=0.2 |
| Thrust (Fn) | 11800 lbf | 11000 lbf | 8000 lbf |
| W (lbm/s) | 147.3 | 142.8 | 119.1 |
| OPR | 13.50 | 12.86 | 12.20 |
| Nmech (rpm) | 8070 | 7944 | 7700 |
| Nozzle Ath | 245.25 in² | 245.25 (same) | 245.25 (same) |

Off-design uses the **same nozzle area** as design; at lower thrust (OD0, OD1) the solver finds lower W and Nmech that still satisfy shaft power and thrust targets.
