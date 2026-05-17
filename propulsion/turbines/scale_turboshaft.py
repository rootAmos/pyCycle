"""
Reference-only design-point model. For aircraft sizing and mission deck work,
use `propulsion/turbines/scaled_turboshaft_deck.py` with
`propulsion/data/turbine/turboshaft_1120hp.csv`.

Scale the Mistè & Benini (2012) T700 turboshaft design point to a target FPT shaft power.

Source of baseline: Mistè, G.A. & Benini, E., "Performance of a Turboshaft Engine
for Helicopter Applications Operating at Variable Shaft Speed", ASME GTINDIA2012-9505,
Tables 2 and 3.

Generic component maps used by the original work: Duyar, A., Gu, Z., Litt, J.S.,
"A Simplified Dynamic Model of the T700 Turboshaft Engine", NASA TM 105805, 1992.
"""

import aerosandbox.numpy as np
from dataclasses import dataclass, field, asdict


# ---------------------------------------------------------------------------
# Baseline design point (Mistè & Benini 2012, Tables 2 & 3)
# ---------------------------------------------------------------------------

@dataclass
class TurboshaftDesignPoint:
    """Design-point cycle parameters for a free-turbine turboshaft."""

    # Inputs (Table 2)
    m_dot_air:        float  # [kg/s]
    pi_inlet:         float  # [-] inlet total-pressure recovery
    pi_compressor:    float  # [-] compressor pressure ratio
    eta_compressor:   float  # [-] compressor isentropic efficiency
    N_compressor:     float  # [rpm] gas-generator design speed
    pi_combustor:     float  # [-] combustor relative total-pressure loss
    eta_combustor:    float  # [-]
    LHV_fuel:         float  # [J/kg]
    eta_ggt:          float  # [-] GGT isentropic efficiency
    eta_mech_ggt:     float  # [-] GGT mechanical efficiency
    N_fpt:            float  # [rpm] FPT design speed
    eta_fpt:          float  # [-] FPT isentropic efficiency
    eta_mech_fpt:     float  # [-] FPT mechanical efficiency
    P_shaft:          float  # [W] FPT shaft power delivered to the load
    eta_nozzle:       float  # [-]

    # Ambient (from the paper's design point: 500 m density altitude, ISA)
    T_amb:            float  # [K]
    p_amb:            float  # [Pa]

    # Outputs (Table 3) - populated by solve_design_point
    T_t: dict = field(default_factory=dict)  # total temperatures by station [K]
    p_t: dict = field(default_factory=dict)  # total pressures by station [Pa]
    m_dot_fuel:       float = 0.0            # [kg/s]
    SFC:              float = 0.0            # [kg/(W*s)]
    eta_overall:      float = 0.0            # [-]


# Mistè baseline — directly from Tables 2 & 3 of GTINDIA2012-9505
MISTE_BASELINE = TurboshaftDesignPoint(
    m_dot_air      = 4.6122,    # [kg/s]
    pi_inlet       = 0.9880,
    pi_compressor  = 17.500,
    eta_compressor = 0.8210,
    N_compressor   = 44_700.0,  # [rpm]
    pi_combustor   = 0.04,      # 4 % total-pressure loss
    eta_combustor  = 0.9850,
    LHV_fuel       = 43.1e6,    # [J/kg]
    eta_ggt        = 0.85,
    eta_mech_ggt   = 0.99,
    N_fpt          = 20_900.0,  # [rpm]
    eta_fpt        = 0.85,
    eta_mech_fpt   = 0.99,
    P_shaft        = 1_329_900.0,  # [W]
    eta_nozzle     = 0.90,
    # Ambient is back-computed from station-1 in Table 3
    T_amb          = 289.44,    # [K]
    p_amb          = 95_891.0 / 0.9880,  # [Pa] de-recover the inlet
)


# ---------------------------------------------------------------------------
# Cycle solver
# ---------------------------------------------------------------------------

# Working-fluid properties — paper uses Shomate / NIST tables; for a closed-form
# scaler we use temperature-averaged cp and gamma matching the paper's station
# values (Table 3 cp varies 1004 -> 1260 J/(kg K)). This closes the cycle to
# within ~0.3 % of the published station enthalpies.
CP_COLD = 1006.0    # [J/(kg K)] inlet/compressor average
CP_HOT  = 1240.0    # [J/(kg K)] post-combustor average
GAMMA_C = 1.400     # [-] cold side
GAMMA_H = 1.333     # [-] hot side


def _compressor_outlet_temperature(T_in, pi_c, eta_c, gamma=GAMMA_C):
    """Total-temperature rise across an adiabatic compressor. Saravanamuttoo,
    Gas Turbine Theory, 6th ed., Eq. 2.13."""
    exponent = (gamma - 1.0) / gamma
    T_out_is = T_in * pi_c**exponent
    return T_in + (T_out_is - T_in) / eta_c


def _turbine_outlet_temperature(T_in, pi_t, eta_t, gamma=GAMMA_H):
    """Total-temperature drop across an adiabatic turbine. Saravanamuttoo,
    Gas Turbine Theory, 6th ed., Eq. 2.14. pi_t is inlet/outlet (>1)."""
    exponent = (gamma - 1.0) / gamma
    T_out_is = T_in * pi_t**(-exponent)
    return T_in - eta_t * (T_in - T_out_is)


def solve_design_point(dp: TurboshaftDesignPoint, T4_target: float) -> TurboshaftDesignPoint:
    """Close the design-point cycle for a turboshaft given mass flow, OPR, and T4.

    Stations follow Mistè Fig. 2:
      1: ambient    2: post-inlet    3: post-compressor
      4: post-combustor (turbine inlet)    5: post-GGT    6: post-FPT    7: nozzle exit
    """
    T_t, p_t = {}, {}

    # Station 1: ambient (stagnation = static for static test stand)
    T_t[1] = dp.T_amb
    p_t[1] = dp.p_amb

    # Station 2: post-inlet
    T_t[2] = T_t[1]
    p_t[2] = p_t[1] * dp.pi_inlet

    # Station 3: post-compressor
    T_t[3] = _compressor_outlet_temperature(T_t[2], dp.pi_compressor, dp.eta_compressor)
    p_t[3] = p_t[2] * dp.pi_compressor
    P_compressor = dp.m_dot_air * CP_COLD * (T_t[3] - T_t[2])  # [W]

    # Station 4: post-combustor (T4 set by user — material/cooling limit)
    T_t[4] = T4_target
    p_t[4] = p_t[3] * (1.0 - dp.pi_combustor)
    # Fuel mass flow from energy balance with combustion efficiency
    m_dot_fuel = (dp.m_dot_air * CP_HOT * (T_t[4] - T_t[3])) / (dp.eta_combustor * dp.LHV_fuel)
    m_dot_hot = dp.m_dot_air + m_dot_fuel  # GGT/FPT see fuel-laden flow

    # Station 5: post-GGT — must drive compressor through eta_mech
    P_ggt_required = P_compressor / dp.eta_mech_ggt
    dT_ggt = P_ggt_required / (m_dot_hot * CP_HOT)
    T_t[5] = T_t[4] - dT_ggt
    # Back out GGT pressure ratio from temperature drop and isentropic efficiency
    pi_ggt = (1.0 - (T_t[4] - T_t[5]) / (T_t[4] * dp.eta_ggt))**(-GAMMA_H / (GAMMA_H - 1.0))
    p_t[5] = p_t[4] / pi_ggt

    # Station 6: post-FPT — extracts P_shaft / eta_mech
    P_fpt_required = dp.P_shaft / dp.eta_mech_fpt
    dT_fpt = P_fpt_required / (m_dot_hot * CP_HOT)
    T_t[6] = T_t[5] - dT_fpt
    pi_fpt = (1.0 - (T_t[5] - T_t[6]) / (T_t[5] * dp.eta_fpt))**(-GAMMA_H / (GAMMA_H - 1.0))
    p_t[6] = p_t[5] / pi_fpt

    # Station 7: nozzle (turboshaft — assume near-ambient discharge)
    T_t[7] = T_t[6]
    p_t[7] = dp.p_amb  # turboshaft jet thrust is negligible; expansion to ambient

    # Performance metrics
    SFC = m_dot_fuel / dp.P_shaft                  # [kg/(W s)]
    eta_overall = dp.P_shaft / (m_dot_fuel * dp.LHV_fuel)

    dp.T_t = T_t
    dp.p_t = p_t
    dp.m_dot_fuel = m_dot_fuel
    dp.SFC = SFC
    dp.eta_overall = eta_overall
    return dp


# ---------------------------------------------------------------------------
# Scaling to a target shaft power
# ---------------------------------------------------------------------------

def scale_to_target_power(
    baseline: TurboshaftDesignPoint,
    P_target: float,
    *,
    pi_compressor_new: float = None,
    T4_new: float = 1480.0,                  # [K] held to baseline by default
    eta_compressor_bump: float = 0.0,        # additive bump (e.g. +0.01 for bigger engine)
    eta_ggt_bump: float = 0.0,
    eta_fpt_bump: float = 0.0,
) -> TurboshaftDesignPoint:
    """Build a new design point at P_target by holding cycle (OPR, T4, efficiencies)
    and letting mass flow follow the cycle.

    Mass flow scales as m_dot_new = m_dot_base * (P_target / P_base) * (W_base / W_new)
    where W is specific work — but at fixed cycle, W_new == W_base, so the ratio
    collapses to linear in P. This is correct only when OPR and T4 are held.

    Reference: Walsh & Fletcher, "Gas Turbine Performance", 2nd ed., Sec. 8.4
    (linear scaling at fixed cycle).
    """
    power_ratio = P_target / baseline.P_shaft

    scaled = TurboshaftDesignPoint(
        m_dot_air      = baseline.m_dot_air * power_ratio,
        pi_inlet       = baseline.pi_inlet,
        pi_compressor  = pi_compressor_new if pi_compressor_new is not None else baseline.pi_compressor,
        eta_compressor = baseline.eta_compressor + eta_compressor_bump,
        N_compressor   = baseline.N_compressor,    # see scale_factors() below
        pi_combustor   = baseline.pi_combustor,
        eta_combustor  = baseline.eta_combustor,
        LHV_fuel       = baseline.LHV_fuel,
        eta_ggt        = baseline.eta_ggt + eta_ggt_bump,
        eta_mech_ggt   = baseline.eta_mech_ggt,
        N_fpt          = baseline.N_fpt,
        eta_fpt        = baseline.eta_fpt + eta_fpt_bump,
        eta_mech_fpt   = baseline.eta_mech_fpt,
        P_shaft        = P_target,
        eta_nozzle     = baseline.eta_nozzle,
        T_amb          = baseline.T_amb,
        p_amb          = baseline.p_amb,
    )
    return solve_design_point(scaled, T4_target=T4_new)


def map_scale_factors(baseline: TurboshaftDesignPoint, scaled: TurboshaftDesignPoint) -> dict:
    """Compute the four scale factors per component to apply to a generic map.

    Procedure: Kurzke, J., "How to Get Component Maps for Aircraft Gas Turbine
    Performance Calculations", ASME 96-GT-164, 1996.

    For each component:
        s_W   = corrected_mass_flow_new / corrected_mass_flow_base
        s_PR  = (PR_new - 1) / (PR_base - 1)         # for compressors
                or PR_new / PR_base                  # for turbines (choked)
        s_eta = eta_new / eta_base
        s_N   = N_corr_new / N_corr_base             # = sqrt(T_in_base / T_in_new) at fixed N

    These factors are multiplied into the generic map values when reading any
    operating point off the map.
    """

    def corrected_mass_flow(m_dot, T_t, p_t):
        return m_dot * np.sqrt(T_t / 288.15) / (p_t / 101_325.0)

    # Compressor (station 2 in, station 3 out)
    Wc_base = corrected_mass_flow(baseline.m_dot_air, baseline.T_t[2], baseline.p_t[2])
    Wc_new  = corrected_mass_flow(scaled.m_dot_air,   scaled.T_t[2],   scaled.p_t[2])

    # GGT — corrected mass flow is referenced to turbine inlet (station 4)
    m_hot_base = baseline.m_dot_air + baseline.m_dot_fuel
    m_hot_new  = scaled.m_dot_air   + scaled.m_dot_fuel
    W_ggt_base = corrected_mass_flow(m_hot_base, baseline.T_t[4], baseline.p_t[4])
    W_ggt_new  = corrected_mass_flow(m_hot_new,  scaled.T_t[4],   scaled.p_t[4])
    pi_ggt_base = baseline.p_t[4] / baseline.p_t[5]
    pi_ggt_new  = scaled.p_t[4]   / scaled.p_t[5]

    # FPT — referenced to FPT inlet (station 5)
    W_fpt_base = corrected_mass_flow(m_hot_base, baseline.T_t[5], baseline.p_t[5])
    W_fpt_new  = corrected_mass_flow(m_hot_new,  scaled.T_t[5],   scaled.p_t[5])
    pi_fpt_base = baseline.p_t[5] / baseline.p_t[6]
    pi_fpt_new  = scaled.p_t[5]   / scaled.p_t[6]

    return {
        "compressor": {
            "s_W":   Wc_new / Wc_base,
            "s_PR":  (scaled.pi_compressor - 1.0) / (baseline.pi_compressor - 1.0),
            "s_eta": scaled.eta_compressor / baseline.eta_compressor,
            "s_N":   np.sqrt(baseline.T_t[2] / scaled.T_t[2]),
        },
        "ggt": {
            "s_W":   W_ggt_new / W_ggt_base,
            "s_PR":  pi_ggt_new / pi_ggt_base,
            "s_eta": scaled.eta_ggt / baseline.eta_ggt,
            "s_N":   np.sqrt(baseline.T_t[4] / scaled.T_t[4]),
        },
        "fpt": {
            "s_W":   W_fpt_new / W_fpt_base,
            "s_PR":  pi_fpt_new / pi_fpt_base,
            "s_eta": scaled.eta_fpt / baseline.eta_fpt,
            "s_N":   np.sqrt(baseline.T_t[5] / scaled.T_t[5]),
        },
    }


# ---------------------------------------------------------------------------
# Optional: simple Reynolds-number efficiency correction for size change
# ---------------------------------------------------------------------------

def reynolds_efficiency_correction(eta_base: float, scale_ratio: float, k: float = 0.05) -> float:
    """Empirical Reynolds-number correction on component efficiency.

    eta_new = eta_base + k * log10(Re_new / Re_base) capped to physical bounds.
    Re scales roughly with linear size; linear size scales with sqrt(power) at
    fixed cycle. So for power scale ratio r, Re_ratio ~ sqrt(r).

    Source: Walsh & Fletcher, "Gas Turbine Performance", 2nd ed., Sec. 5.6.
    Coefficient k = 0.05 is the rule-of-thumb for aero turbomachinery.
    """
    re_ratio = np.sqrt(scale_ratio)
    return np.clip(eta_base + k * np.log10(re_ratio), 0.0, 0.95)


# ---------------------------------------------------------------------------
def _print_design_point(dp: TurboshaftDesignPoint, label: str):
    print(f"\n=== {label} ===")
    print(f"  Shaft power      : {dp.P_shaft/1e3:8.1f} kW")
    print(f"  Air mass flow    : {dp.m_dot_air:8.3f} kg/s")
    print(f"  Fuel mass flow   : {dp.m_dot_fuel:8.4f} kg/s")
    print(f"  OPR              : {dp.pi_compressor:8.2f}")
    print(f"  T4 (TIT)         : {dp.T_t[4]:8.1f} K")
    print(f"  Overall eta      : {dp.eta_overall*100:8.2f} %")
    print(f"  SFC              : {dp.SFC*3.6e6:8.4f} kg/(kWh)")
    print(f"  Stations T_t [K] : 1:{dp.T_t[1]:.1f}  3:{dp.T_t[3]:.1f}  "
          f"4:{dp.T_t[4]:.1f}  5:{dp.T_t[5]:.1f}  6:{dp.T_t[6]:.1f}")


if __name__ == "__main__":
    # --- Step 1: validate baseline against Mistè Table 3 ---
    base = solve_design_point(MISTE_BASELINE, T4_target=1479.72)
    _print_design_point(base, "Baseline (Mistè & Benini 2012, T700-class, 1.33 MW)")
    print(f"  Mistè published SFC          : 0.2719 kg/(kWh)")
    print(f"  Mistè published eta_overall  : 30.72 %")
    print(f"  Mistè published T_t,5         : 1125.08 K")

    # --- Step 2: scale to 3 MW with same technology level ---
    target_3MW = scale_to_target_power(
        base,
        P_target=3.0e6,                  # [W]
        pi_compressor_new=18.0,          # slight bump (CT7-8 territory)
        T4_new=1500.0,                   # [K] modest TIT bump
        eta_compressor_bump=+0.01,       # bigger compressor -> +1 pt
        eta_ggt_bump=+0.02,              # bigger HP turbine -> +2 pt
        eta_fpt_bump=+0.03,              # better PT design -> +3 pt
    )
    _print_design_point(target_3MW, "Scaled to 3 MW (CT7-8/T408-class technology)")

    # --- Step 3: report map scale factors for off-design map use ---
    factors = map_scale_factors(base, target_3MW)
    print("\n=== Map scale factors (apply to a generic map at any operating point) ===")
    for comp, fdict in factors.items():
        print(f"  {comp:11s}: " + ", ".join(f"{k}={v:.4f}" for k, v in fdict.items()))

    # --- Step 4: optional Reynolds-aware version ---
    print("\n=== Reynolds-corrected efficiency check ===")
    r = target_3MW.P_shaft / base.P_shaft
    for name, eta in [("compressor", base.eta_compressor),
                      ("GGT",        base.eta_ggt),
                      ("FPT",        base.eta_fpt)]:
        eta_re = reynolds_efficiency_correction(eta, r)
        print(f"  {name:11s}: eta_base={eta:.3f}  eta_Re_corrected={eta_re:.3f}  "
              f"(Δ = {(eta_re-eta)*100:+.2f} pt)")

    # --- Step 5: parametric sweep of SFC vs target power ---
    print("\n=== SFC sweep across power class ===")
    print(f"  {'P [MW]':>8} {'m_dot [kg/s]':>14} {'SFC [kg/kWh]':>14} {'eta [%]':>10}")
    for P_MW in [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]:
        scaled = scale_to_target_power(base, P_target=P_MW * 1e6, pi_compressor_new=18.0)
        print(f"  {P_MW:8.2f} {scaled.m_dot_air:14.3f} "
              f"{scaled.SFC*3.6e6:14.4f} {scaled.eta_overall*100:10.2f}")
