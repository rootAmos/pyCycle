"""Shared constraint-analysis equations."""


def design_point_thrust_to_weight_from_wing_loading(
    wing_loading_N_m2,
    dynamic_pressure_Pa,
    velocity_m_s,
    installed_full_throttle_thrust_lapse=1.0,
    instantaneous_weight_fraction=1.0,
    load_factor=1.0,
    drag_polar_k1=0.0,
    drag_polar_k2=0.0,
    zero_lift_drag_coefficient=0.0,
    climb_rate_m_s=0.0,
    acceleration_m_s2=0.0,
    gravity_m_s2=9.80665,
):
    """Return design T/W for a condition at a given wing loading W/S."""
    beta = instantaneous_weight_fraction
    alpha = installed_full_throttle_thrust_lapse
    wing_loading_term = dynamic_pressure_Pa / (beta * wing_loading_N_m2)
    lift_loading_term = load_factor * beta * wing_loading_N_m2 / dynamic_pressure_Pa
    excess_power_m_s = climb_rate_m_s + velocity_m_s * acceleration_m_s2 / gravity_m_s2

    return (
        beta
        / alpha
        * (
            wing_loading_term
            * (
                drag_polar_k1 * lift_loading_term**2
                + drag_polar_k2 * lift_loading_term
                + zero_lift_drag_coefficient
            )
            + excess_power_m_s / velocity_m_s
        )
    )
