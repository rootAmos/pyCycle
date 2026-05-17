"""Aero drag build-up and engine-deck drag-point helpers."""

import aerosandbox as asb
import aerosandbox.numpy as np
import aerosandbox.tools.units as u

from aero.data.interpolators import (
    cla_cla_theory_ratio as airfoil_cla_theory_ratio,
    leading_edge_suction_factor,
)


def swept_wing_tip_le_x_m(span_m, root_chord_m, tip_chord_m, sweep_25_rad):
    """Return tip leading-edge x offset from quarter-chord sweep."""
    return (
        0.5 * span_m * np.tan(sweep_25_rad)
        + 0.25 * (root_chord_m - tip_chord_m)
    )


def drag_geometry_from_planform_area(planform_area_m2, config, number_engines):
    """Return sizing geometry needed for parasite and wave drag build-up."""
    aspect_ratio = 3.0
    taper_ratio = 0.25
    root_thickness_to_chord = 0.06
    span_m = (planform_area_m2 * aspect_ratio) ** 0.5
    root_chord_m = 2.0 * planform_area_m2 / (span_m * (1.0 + taper_ratio))
    tip_chord_m = taper_ratio * root_chord_m
    main_wing_tip_le_x_m = swept_wing_tip_le_x_m(
        span_m,
        root_chord_m,
        tip_chord_m,
        config.main_wing_quarter_chord_sweep_rad,
    )
    leading_edge_sweep_rad = np.arctan(main_wing_tip_le_x_m / (0.5 * span_m))
    half_chord_sweep_rad = np.arctan(
        (
            main_wing_tip_le_x_m
            + 0.5 * tip_chord_m
            - 0.5 * root_chord_m
        )
        / (0.5 * span_m)
    )
    fuselage_length_m = 4.0 * planform_area_m2**0.5
    fuselage_height_m = 0.12 * fuselage_length_m
    fuselage_width_m = 0.10 * fuselage_length_m
    fuselage_radius_a_m = 0.5 * fuselage_width_m
    fuselage_radius_b_m = 0.5 * fuselage_height_m
    fuselage_perimeter_m = np.pi * (
        3.0 * (fuselage_radius_a_m + fuselage_radius_b_m)
        - (
            (3.0 * fuselage_radius_a_m + fuselage_radius_b_m)
            * (fuselage_radius_a_m + 3.0 * fuselage_radius_b_m)
        )
        ** 0.5
    )
    equivalent_fuselage_diameter_m = (fuselage_height_m * fuselage_width_m) ** 0.5
    engine_diameter_m = 2.0 * u.foot
    engine_length_m = config.nacelle_length_to_diameter * engine_diameter_m
    vtail_area_m2 = 0.26 * planform_area_m2

    return {
        "reference_area_m2": planform_area_m2,
        "aspect_ratio": aspect_ratio,
        "root_thickness_to_chord": root_thickness_to_chord,
        "tip_chord_m": tip_chord_m,
        "mean_aerodynamic_chord_m": 2.0
        / 3.0
        * root_chord_m
        * (1.0 + taper_ratio + taper_ratio**2)
        / (1.0 + taper_ratio),
        "leading_edge_sweep_rad": leading_edge_sweep_rad,
        "half_chord_sweep_rad": half_chord_sweep_rad,
        "wing_wetted_area_m2": 2.0
        * planform_area_m2
        * (1.0 + 0.25 * root_thickness_to_chord),
        "tail_wetted_area_m2": 2.0
        * vtail_area_m2
        * (1.0 + 0.25 * root_thickness_to_chord),
        "tail_mean_chord_m": vtail_area_m2 / (vtail_area_m2 * 1.4) ** 0.5,
        "fuselage_length_m": fuselage_length_m,
        "fuselage_wetted_area_m2": fuselage_perimeter_m * fuselage_length_m,
        "fuselage_fineness_ratio": fuselage_length_m / equivalent_fuselage_diameter_m,
        "max_cross_section_area_m2": 0.25
        * np.pi
        * fuselage_width_m
        * fuselage_height_m,
        "nacelle_wetted_area_m2": number_engines
        * np.pi
        * engine_diameter_m
        * engine_length_m,
        "nacelle_length_m": engine_length_m,
    }


def air_dynamic_viscosity_kg_m_s(temperature_K):
    """Sutherland-law dynamic viscosity for air."""
    reference_temperature_K = 273.15
    reference_viscosity_kg_m_s = 1.716e-5
    sutherland_temperature_K = 110.4
    return (
        reference_viscosity_kg_m_s
        * (temperature_K / reference_temperature_K) ** 1.5
        * (reference_temperature_K + sutherland_temperature_K)
        / (temperature_K + sutherland_temperature_K)
    )


def turbulent_skin_friction_coefficient(reynolds_number, mach):
    """Raymer-style turbulent flat-plate skin friction coefficient."""
    reynolds_number = np.maximum(reynolds_number, 1.0e5)
    return 0.455 / (
        np.log10(reynolds_number) ** 2.58 * (1.0 + 0.144 * mach**2) ** 0.65
    )


def swept_wing_oswald_efficiency(aspect_ratio, leading_edge_sweep_rad):
    """Raymer swept-wing Oswald efficiency correlation for Lambda_LE > 30 deg."""
    return (
        4.61
        * (1.0 - 0.045 * aspect_ratio**0.68)
        * np.cos(leading_edge_sweep_rad) ** 0.15
        - 3.1
    )


def smoothstep(x):
    x = np.clip(x, 0.0, 1.0)
    return x**2.0 * (3.0 - 2.0 * x)


def airfoil_theory_lift_curve_slope(config, drag_geometry):
    """Return theoretical 2D airfoil lift curve slope in 1/rad."""
    thickness_to_chord = drag_geometry["root_thickness_to_chord"]
    return (
        2.0 * np.pi
        + 4.7
        * thickness_to_chord
        * (1.0 + 0.00375 * config.airfoil_trailing_edge_angle_deg)
    )


def subsonic_finite_wing_lift_curve_slope(
    config,
    drag_geometry,
    mach,
    reynolds_number,
):
    """Return 3D subsonic CL_alpha using the airfoil-ratio data and finite-wing relation."""
    aspect_ratio = drag_geometry["aspect_ratio"]
    beta = np.sqrt(np.maximum(1.0 - mach**2.0, 1.0e-9))
    tan_half_te_ang = np.tan(np.radians(0.5 * config.airfoil_trailing_edge_angle_deg))
    clalpha_theory = airfoil_theory_lift_curve_slope(config, drag_geometry)
    clalpha_ratio = airfoil_cla_theory_ratio(
        tan_half_te_ang_deg=tan_half_te_ang,
        reynolds_number=reynolds_number,
    )
    airfoil_kappa = 1.05 * clalpha_ratio * clalpha_theory / (2.0 * np.pi)
    return (
        2.0
        * np.pi
        * aspect_ratio
        / (
            2.0
            + np.sqrt(
                aspect_ratio**2.0
                * beta**2.0
                / airfoil_kappa**2.0
                * (
                    1.0
                    + np.tan(drag_geometry["half_chord_sweep_rad"]) ** 2.0
                    / beta**2.0
                )
                + 4.0
            )
        )
    )


def supersonic_ackeret_lift_curve_slope(mach):
    """Return Ackeret 2D supersonic lift curve slope in 1/rad."""
    return 4.0 / np.sqrt(np.maximum(mach**2.0 - 1.0, 1.0e-9))


def blended_lift_curve_slope(config, drag_geometry, mach, reynolds_number):
    """Smoothly blend subsonic finite-wing CL_alpha to Ackeret CL_alpha."""
    subsonic_clalpha = subsonic_finite_wing_lift_curve_slope(
        config=config,
        drag_geometry=drag_geometry,
        mach=mach,
        reynolds_number=reynolds_number,
    )
    m_start = 1.0
    m_end = leading_edge_sonic_mach(drag_geometry)
    m_span = np.maximum(m_end - m_start, 1.0e-6)
    transonic_target_clalpha = supersonic_ackeret_lift_curve_slope(m_end)
    supersonic_clalpha = supersonic_ackeret_lift_curve_slope(np.maximum(mach, m_end))
    blend = smoothstep((mach - m_start) / m_span)
    transonic_clalpha = (
        (1.0 - blend) * subsonic_clalpha
        + blend * transonic_target_clalpha
    )
    return np.where(mach < m_end, transonic_clalpha, supersonic_clalpha)


def leading_edge_sonic_mach(drag_geometry):
    """Return Mach where the leading-edge normal component becomes sonic."""
    return 1.0 / np.maximum(np.cos(drag_geometry["leading_edge_sweep_rad"]), 1.0e-9)


def lift_dependent_drag_factor(
    config,
    drag_geometry,
    mach,
    reynolds_number,
    lift_coefficient,
):
    """Return K from leading-edge suction split between K100 and K0."""
    subsonic_clalpha = subsonic_finite_wing_lift_curve_slope(
        config=config,
        drag_geometry=drag_geometry,
        mach=mach,
        reynolds_number=reynolds_number,
    )
    m_start = 1.0
    m_end = leading_edge_sonic_mach(drag_geometry)
    m_span = np.maximum(m_end - m_start, 1.0e-6)
    supersonic_clalpha_at_transition = supersonic_ackeret_lift_curve_slope(m_end)
    supersonic_clalpha = supersonic_ackeret_lift_curve_slope(np.maximum(mach, m_end))
    design_cl = drag_geometry.get("design_lift_coefficient", 0.3)
    suction = leading_edge_suction_factor(
        cl=np.maximum(lift_coefficient, 0.0),
        cl_design=design_cl,
    )
    suction = np.clip(suction, 0.0, 1.0)
    aspect_ratio = drag_geometry["aspect_ratio"]
    k100 = 1.0 / (np.pi * aspect_ratio)
    subsonic_k = suction * k100 + (1.0 - suction) / subsonic_clalpha
    transition_supersonic_k = (
        suction * k100 + (1.0 - suction) / supersonic_clalpha_at_transition
    )
    supersonic_k = suction * k100 + (1.0 - suction) / supersonic_clalpha
    blend = smoothstep((mach - m_start) / m_span)
    transonic_k = (
        (1.0 - blend) * subsonic_k
        + blend * transition_supersonic_k
    )
    lift_dependent_k = np.where(mach < m_end, transonic_k, supersonic_k)
    clalpha = blended_lift_curve_slope(
        config=config,
        drag_geometry=drag_geometry,
        mach=mach,
        reynolds_number=reynolds_number,
    )
    return lift_dependent_k, clalpha, suction


def supersonic_wave_drag_coefficient(config, drag_geometry, mach):
    """Return the Ma >= 1.2 wave drag estimate from the supplied paper."""
    leading_edge_sweep_deg = np.degrees(drag_geometry["leading_edge_sweep_rad"])
    wave_factor = (
        1.0
        - 0.386
        * np.maximum(mach - config.supersonic_wave_drag_start_mach, 0.0) ** 0.57
        * (1.0 - np.pi * leading_edge_sweep_deg / 100.0) ** 2.0
    )
    return (
        1.5
        * np.maximum(wave_factor, 0.0)
        * 9.0
        * np.pi
        / 2.0
        * (drag_geometry["max_cross_section_area_m2"] / drag_geometry["fuselage_length_m"]) ** 2.0
        / drag_geometry["reference_area_m2"]
    )


def transonic_wave_drag_coefficient(config, drag_geometry, mach):
    """Bezier drag-rise estimate between Mcrit and the Ma 1.2 wave-drag model."""
    mdd = config.drag_divergence_mach
    mcrit = mdd - config.critical_mach_offset_from_mdd
    msup = config.supersonic_wave_drag_start_mach
    cdw_mdd = config.drag_divergence_wave_cd
    cdw_msup = supersonic_wave_drag_coefficient(config, drag_geometry, msup)

    t = np.clip((mach - mcrit) / (msup - mcrit), 0.0, 1.0)
    t_mdd = (mdd - mcrit) / (msup - mcrit)

    # Cubic Bezier ordinate. P0 is zero at Mcrit, P2 has the same CDw as P3
    # so the curve reaches the Ma 1.2 value with a flat tangent as in points B/A.
    p0 = 0.0
    p2 = cdw_msup
    p3 = cdw_msup
    p1_denominator = 3.0 * (1.0 - t_mdd) ** 2.0 * t_mdd
    p1_numerator = cdw_mdd - (
        3.0 * (1.0 - t_mdd) * t_mdd**2.0 * p2 + t_mdd**3.0 * p3
    )
    p1 = p1_numerator / p1_denominator
    wave_cd = (
        (1.0 - t) ** 3.0 * p0
        + 3.0 * (1.0 - t) ** 2.0 * t * p1
        + 3.0 * (1.0 - t) * t**2.0 * p2
        + t**3.0 * p3
    )
    return np.maximum(wave_cd, 0.0)


def drag_build_up_coefficients(
    config,
    drag_geometry,
    altitude_m,
    velocity_m_s,
    lift_coefficient=None,
):
    """Return condition-dependent CD0 and lift-dependent K."""
    atmosphere = asb.Atmosphere(altitude=altitude_m)
    density_kg_m3 = atmosphere.density()
    speed_of_sound_m_s = atmosphere.speed_of_sound()
    mach = velocity_m_s / speed_of_sound_m_s
    viscosity_kg_m_s = air_dynamic_viscosity_kg_m_s(atmosphere.temperature())

    reference_area_m2 = drag_geometry["reference_area_m2"]
    wing_re = density_kg_m3 * velocity_m_s * drag_geometry["mean_aerodynamic_chord_m"] / viscosity_kg_m_s
    tail_re = density_kg_m3 * velocity_m_s * drag_geometry["tail_mean_chord_m"] / viscosity_kg_m_s
    fuselage_re = density_kg_m3 * velocity_m_s * drag_geometry["fuselage_length_m"] / viscosity_kg_m_s
    nacelle_re = density_kg_m3 * velocity_m_s * drag_geometry["nacelle_length_m"] / viscosity_kg_m_s

    wing_cd0 = (
        turbulent_skin_friction_coefficient(wing_re, mach)
        * config.wing_form_factor
        * drag_geometry["wing_wetted_area_m2"]
        / reference_area_m2
    )
    tail_cd0 = (
        turbulent_skin_friction_coefficient(tail_re, mach)
        * config.tail_form_factor
        * drag_geometry["tail_wetted_area_m2"]
        / reference_area_m2
    )
    fuselage_form_factor = (
        1.0
        + 60.0 / drag_geometry["fuselage_fineness_ratio"] ** 3.0
        + drag_geometry["fuselage_fineness_ratio"] / 400.0
    )
    fuselage_cd0 = (
        turbulent_skin_friction_coefficient(fuselage_re, mach)
        * fuselage_form_factor
        * drag_geometry["fuselage_wetted_area_m2"]
        / reference_area_m2
    )
    nacelle_cd0 = (
        turbulent_skin_friction_coefficient(nacelle_re, mach)
        * config.nacelle_form_factor
        * drag_geometry["nacelle_wetted_area_m2"]
        / reference_area_m2
    )
    parasite_cd0 = wing_cd0 + tail_cd0 + fuselage_cd0 + nacelle_cd0

    transonic_wave_cd0 = transonic_wave_drag_coefficient(
        config,
        drag_geometry,
        mach,
    )
    supersonic_wave_cd0 = supersonic_wave_drag_coefficient(config, drag_geometry, mach)
    wave_cd0 = np.where(
        mach < config.drag_divergence_mach - config.critical_mach_offset_from_mdd,
        0.0,
        np.where(
            mach < config.supersonic_wave_drag_start_mach,
            transonic_wave_cd0,
            supersonic_wave_cd0,
        ),
    )

    if lift_coefficient is None:
        lift_coefficient = drag_geometry.get("design_lift_coefficient", 0.3)
    lift_dependent_k, lift_curve_slope, leading_edge_suction = lift_dependent_drag_factor(
        config=config,
        drag_geometry=drag_geometry,
        mach=mach,
        reynolds_number=wing_re,
        lift_coefficient=lift_coefficient,
    )
    oswald_efficiency = swept_wing_oswald_efficiency(
        drag_geometry["aspect_ratio"],
        drag_geometry["leading_edge_sweep_rad"],
    )

    return {
        "mach": mach,
        "parasite_cd0": parasite_cd0,
        "wave_cd0": wave_cd0,
        "zero_lift_drag_coefficient": parasite_cd0 + wave_cd0,
        "oswald_efficiency": oswald_efficiency,
        "lift_curve_slope": lift_curve_slope,
        "leading_edge_suction": leading_edge_suction,
        "lift_dependent_k": lift_dependent_k,
    }


def engine_deck_drag_point(aircraft_sizing, mach, altitude_m):
    """Return atmospheric and drag-derived thrust metrics for one deck condition."""
    config = aircraft_sizing["config"]
    atmosphere = asb.Atmosphere(altitude=altitude_m)
    velocity_m_s = mach * atmosphere.speed_of_sound()
    dynamic_pressure_Pa = 0.5 * atmosphere.density() * velocity_m_s**2
    lift_coefficient = (
        aircraft_sizing["takeoff_weight_N"]
        / (dynamic_pressure_Pa * aircraft_sizing["planform_area_m2"])
    )
    drag = drag_build_up_coefficients(
        config=config,
        drag_geometry=aircraft_sizing["drag_geometry"],
        altitude_m=altitude_m,
        velocity_m_s=velocity_m_s,
        lift_coefficient=lift_coefficient,
    )
    induced_drag_coefficient = drag["lift_dependent_k"] * lift_coefficient**2.0
    total_drag_coefficient = (
        drag["zero_lift_drag_coefficient"] + induced_drag_coefficient
    )
    drag_N = (
        dynamic_pressure_Pa
        * aircraft_sizing["planform_area_m2"]
        * total_drag_coefficient
    )
    return {
        "mach": mach,
        "altitude_m": altitude_m,
        "temperature_K": atmosphere.temperature(),
        "pressure_Pa": atmosphere.pressure(),
        "density_kg_m3": atmosphere.density(),
        "speed_of_sound_m_s": atmosphere.speed_of_sound(),
        "velocity_m_s": velocity_m_s,
        "dynamic_pressure_Pa": dynamic_pressure_Pa,
        "lift_coefficient": lift_coefficient,
        "zero_lift_drag_coefficient": drag["zero_lift_drag_coefficient"],
        "parasite_cd0": drag["parasite_cd0"],
        "wave_cd0": drag["wave_cd0"],
        "lift_dependent_k": drag["lift_dependent_k"],
        "induced_drag_coefficient": induced_drag_coefficient,
        "total_drag_coefficient": total_drag_coefficient,
        "drag_N": drag_N,
        "required_thrust_N": drag_N,
        "required_thrust_to_weight": drag_N / aircraft_sizing["takeoff_weight_N"],
        "is_stall_limited": lift_coefficient > config.drag_plot_cl_max,
    }
