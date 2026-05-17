"""Parametric motor/generator efficiency model for the Duality electric fans.

The loss model follows the AIAA motor-map equations supplied by the user:

    P_L = C0 + C1 omega + C2 omega^3 + C3 Q^2

The coefficients are parameterized by the peak-efficiency speed, torque, and
efficiency plus a parasite-loss ratio k0. The same loss model is used for
motors and generators; generators are solved with a small fixed-point loop
because the required shaft power appears on both sides of the efficiency
relationship.
"""

from dataclasses import dataclass
import csv
import json
from pathlib import Path

import aerosandbox.numpy as np


hp_to_W = 745.6998715822702
W_to_hp = 1.0 / hp_to_W
gamma_air = 1.4
gas_constant_ft_lbf_slug_R = 1716.59
standard_pressure_psf = 2116.21662367394
standard_temperature_R = 518.67
lbm_per_slug = 32.1740485564


@dataclass(frozen=True)
class ElectricMachineMap:
    peak_speed_rad_s: object
    peak_torque_N_m: object
    peak_efficiency: object = 0.96
    parasite_loss_ratio: object = 0.35
    rated_torque_ratio: object = 1.5
    rated_power_ratio: object = 1.5
    speed_limit_ratio: object = 1.5


def electric_machine_loss_coefficients(machine: ElectricMachineMap):
    """Return C0, C1, C2, C3 for the power-loss map."""
    omega_hat = machine.peak_speed_rad_s
    torque_hat = machine.peak_torque_N_m
    eta_hat = machine.peak_efficiency

    # Source: supplied motor-map literature, Eqs. 7 and 9.
    c0 = machine.parasite_loss_ratio * omega_hat * torque_hat / 6.0 * (1.0 - eta_hat) / eta_hat
    c1 = -3.0 * c0 / (2.0 * omega_hat) + torque_hat * (1.0 - eta_hat) / (4.0 * eta_hat)
    c2 = c0 / (2.0 * omega_hat**3) + torque_hat * (1.0 - eta_hat) / (4.0 * eta_hat * omega_hat**2)
    c3 = omega_hat * (1.0 - eta_hat) / (2.0 * torque_hat * eta_hat)
    return c0, c1, c2, c3


def electric_machine_loss_power_W(machine: ElectricMachineMap, speed_rad_s, torque_N_m):
    """Return loss power P_L in W."""
    c0, c1, c2, c3 = electric_machine_loss_coefficients(machine)
    # Source: supplied motor-map literature, Eq. 1.
    return c0 + c1 * speed_rad_s + c2 * speed_rad_s**3 + c3 * torque_N_m**2


def electric_machine_efficiency(machine: ElectricMachineMap, speed_rad_s, torque_N_m):
    """Return machine efficiency at a speed/torque operating point."""
    shaft_power_W = speed_rad_s * torque_N_m
    loss_power_W = electric_machine_loss_power_W(machine, speed_rad_s, torque_N_m)
    # Source: supplied motor-map literature, Eq. 2.
    return np.where(
        shaft_power_W > 0.0,
        shaft_power_W / (shaft_power_W + loss_power_W),
        0.0,
    )


def motor_electrical_power_W(machine: ElectricMachineMap, shaft_power_W, speed_rad_s):
    """Return electrical input power for a motor driving a fan."""
    if shaft_power_W <= 0.0:
        return 0.0
    torque_N_m = shaft_power_W / speed_rad_s
    loss_power_W = electric_machine_loss_power_W(machine, speed_rad_s, torque_N_m)
    return shaft_power_W + loss_power_W


def generator_shaft_power_W(machine: ElectricMachineMap, electrical_output_W, speed_rad_s, iterations=12):
    """Return shaft power required for a generator to supply electrical output."""
    if electrical_output_W <= 0.0:
        return 0.0
    shaft_power_W = electrical_output_W / machine.peak_efficiency
    for _ in range(iterations):
        torque_N_m = shaft_power_W / speed_rad_s
        loss_power_W = electric_machine_loss_power_W(machine, speed_rad_s, torque_N_m)
        shaft_power_W = electrical_output_W + loss_power_W
    return shaft_power_W


def electric_machine_limits(machine: ElectricMachineMap):
    """Return rated torque, rated power, speed limit, and rated speed."""
    # Source: supplied motor-map literature, Eqs. 10-13.
    rated_torque_N_m = machine.rated_torque_ratio * machine.peak_torque_N_m
    rated_power_W = machine.rated_power_ratio * machine.peak_speed_rad_s * machine.peak_torque_N_m
    speed_limit_rad_s = machine.speed_limit_ratio * machine.peak_speed_rad_s
    rated_speed_rad_s = machine.rated_power_ratio / machine.rated_torque_ratio * machine.peak_speed_rad_s
    return {
        "rated_torque_N_m": rated_torque_N_m,
        "rated_power_W": rated_power_W,
        "speed_limit_rad_s": speed_limit_rad_s,
        "rated_speed_rad_s": rated_speed_rad_s,
    }


def rpm_to_rad_s(speed_rpm):
    return speed_rpm * 2.0 * np.pi / 60.0


def electric_machine_for_required_power(
    required_power_W,
    peak_speed_rpm,
    peak_efficiency=0.96,
    parasite_loss_ratio=0.35,
    rated_power_ratio=1.5,
    rated_torque_ratio=1.5,
    speed_limit_ratio=1.5,
    sizing_margin=1.0,
):
    """Scale an electric-machine map so rated power clears the required power."""
    peak_speed_rad_s = rpm_to_rad_s(peak_speed_rpm)
    peak_torque_N_m = (
        sizing_margin
        * required_power_W
        / (rated_power_ratio * peak_speed_rad_s)
    )
    return ElectricMachineMap(
        peak_speed_rad_s=peak_speed_rad_s,
        peak_torque_N_m=max(peak_torque_N_m, 1.0e-9),
        peak_efficiency=peak_efficiency,
        parasite_loss_ratio=parasite_loss_ratio,
        rated_torque_ratio=rated_torque_ratio,
        rated_power_ratio=rated_power_ratio,
        speed_limit_ratio=speed_limit_ratio,
    )


def _row_float(row, name, default=0.0):
    value = row.get(name)
    if value in (None, ""):
        return default
    return float(value)


def _machine_summary(machine, prefix):
    limits = electric_machine_limits(machine)
    return {
        f"{prefix}_peak_speed_rpm": float(machine.peak_speed_rad_s * 60.0 / (2.0 * np.pi)),
        f"{prefix}_peak_torque_N_m": float(machine.peak_torque_N_m),
        f"{prefix}_peak_efficiency": float(machine.peak_efficiency),
        f"{prefix}_rated_power_hp": float(limits["rated_power_W"] * W_to_hp),
        f"{prefix}_rated_power_kw": float(limits["rated_power_W"] / 1000.0),
        f"{prefix}_rated_torque_N_m": float(limits["rated_torque_N_m"]),
    }


def _equivalent_diameter_in(area_in2):
    if area_in2 in ("", None) or area_in2 <= 0.0:
        return ""
    return float((4.0 * area_in2 / np.pi) ** 0.5)


def _sizing_flowpath_summary(row):
    summary = {}
    for prefix, column in (
        ("inlet", "inlet_area_in2"),
        ("fan1", "fan1_area_in2"),
        ("fan2", "fan2_area_in2"),
        ("nozzle_throat", "nozzle_throat_area_in2"),
    ):
        area_in2 = "" if column not in row or row.get(column) in (None, "") else _row_float(row, column)
        summary[f"sizing_{prefix}_area_in2"] = float(area_in2) if area_in2 != "" else ""
        summary[f"sizing_{prefix}_equivalent_diameter_in"] = _equivalent_diameter_in(area_in2)
    return summary


def standard_atmosphere_imperial(altitude_ft):
    altitude_m = float(altitude_ft) * 0.3048
    g0 = 9.80665
    gas_constant = 287.05287
    temperature_K = 288.15
    pressure_Pa = 101325.0
    base_altitude_m = 0.0
    layers = ((11000.0, -0.0065), (20000.0, 0.0), (32000.0, 0.0010))

    for top_altitude_m, lapse_rate in layers:
        next_altitude_m = min(altitude_m, top_altitude_m)
        delta_h = next_altitude_m - base_altitude_m
        if delta_h > 0.0:
            if abs(lapse_rate) < 1.0e-12:
                pressure_Pa *= np.exp(-g0 * delta_h / (gas_constant * temperature_K))
            else:
                next_temperature_K = temperature_K + lapse_rate * delta_h
                pressure_Pa *= (next_temperature_K / temperature_K) ** (-g0 / (lapse_rate * gas_constant))
                temperature_K = next_temperature_K
        base_altitude_m = next_altitude_m
        if altitude_m <= top_altitude_m:
            break

    return {
        "temperature_R": float(temperature_K * 1.8),
        "pressure_psf": float(pressure_Pa * 0.020885434273039),
        "density_slug_ft3": float(pressure_Pa / (gas_constant * temperature_K) * 0.001940320331954),
    }


def generator_turbine_pressure_recovery(
    mach,
    subsonic_pressure_recovery=0.98,
    min_pressure_recovery=0.35,
    supersonic_recovery_coefficient=0.075,
    supersonic_recovery_exponent=1.35,
):
    if mach <= 1.0:
        return float(subsonic_pressure_recovery)
    recovery = subsonic_pressure_recovery - supersonic_recovery_coefficient * (mach - 1.0) ** supersonic_recovery_exponent
    return float(np.clip(recovery, min_pressure_recovery, subsonic_pressure_recovery))


def generator_turbine_inlet_state(
    mach,
    altitude_ft,
    engine_face_mach=0.35,
    subsonic_pressure_recovery=0.98,
    min_pressure_recovery=0.35,
    supersonic_recovery_coefficient=0.075,
    supersonic_recovery_exponent=1.35,
):
    atmosphere = standard_atmosphere_imperial(altitude_ft)
    total_temperature_R = atmosphere["temperature_R"] * (1.0 + 0.5 * (gamma_air - 1.0) * mach**2)
    freestream_total_pressure_psf = atmosphere["pressure_psf"] * (
        1.0 + 0.5 * (gamma_air - 1.0) * mach**2
    ) ** (gamma_air / (gamma_air - 1.0))
    pressure_recovery = generator_turbine_pressure_recovery(
        mach,
        subsonic_pressure_recovery=subsonic_pressure_recovery,
        min_pressure_recovery=min_pressure_recovery,
        supersonic_recovery_coefficient=supersonic_recovery_coefficient,
        supersonic_recovery_exponent=supersonic_recovery_exponent,
    )
    total_pressure_psf = pressure_recovery * freestream_total_pressure_psf
    static_temperature_R = total_temperature_R / (1.0 + 0.5 * (gamma_air - 1.0) * engine_face_mach**2)
    static_pressure_psf = total_pressure_psf / (
        1.0 + 0.5 * (gamma_air - 1.0) * engine_face_mach**2
    ) ** (gamma_air / (gamma_air - 1.0))
    density_slug_ft3 = static_pressure_psf / (gas_constant_ft_lbf_slug_R * static_temperature_R)
    velocity_ft_s = engine_face_mach * np.sqrt(gamma_air * gas_constant_ft_lbf_slug_R * static_temperature_R)
    return {
        "freestream_total_temperature_R": float(total_temperature_R),
        "freestream_total_pressure_psf": float(freestream_total_pressure_psf),
        "engine_face_mach": float(engine_face_mach),
        "pressure_recovery": float(pressure_recovery),
        "total_temperature_R": float(total_temperature_R),
        "total_pressure_psf": float(total_pressure_psf),
        "static_temperature_R": float(static_temperature_R),
        "static_pressure_psf": float(static_pressure_psf),
        "density_slug_ft3": float(density_slug_ft3),
        "velocity_ft_s": float(velocity_ft_s),
    }


def corrected_airflow_lbm_s(turboshaft_deck):
    metadata = turboshaft_deck.metadata or {}
    base_corrected_airflow_lbm_s = float(metadata.get("sls_corrected_airflow", 0.0) or 0.0)
    if base_corrected_airflow_lbm_s <= 0.0:
        return 0.0
    return base_corrected_airflow_lbm_s * float(turboshaft_deck.scale_factor)


def generator_turbine_inlet_sizing(inlet_state, corrected_airflow_lbm_s):
    theta = inlet_state["total_temperature_R"] / standard_temperature_R
    delta = inlet_state["total_pressure_psf"] / standard_pressure_psf
    airflow_lbm_s = corrected_airflow_lbm_s * delta / np.sqrt(theta)
    area_ft2 = (
        airflow_lbm_s
        / lbm_per_slug
        / max(inlet_state["density_slug_ft3"] * inlet_state["velocity_ft_s"], 1.0e-12)
    )
    return {
        "corrected_airflow_lbm_s": float(corrected_airflow_lbm_s),
        "airflow_lbm_s": float(airflow_lbm_s),
        "engine_face_area_ft2": float(area_ft2),
        "engine_face_area_in2": float(area_ft2 * 144.0),
    }


def write_sizing_summary(summary, output_json, output_csv=None):
    output_json = Path(output_json)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    with output_json.open("w", encoding="utf-8") as stream:
        json.dump(summary, stream, indent=2)
    if output_csv is not None:
        output_csv = Path(output_csv)
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        with output_csv.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=("quantity", "value"))
            writer.writeheader()
            writer.writerows(
                {"quantity": key, "value": value}
                for key, value in summary.items()
            )


def write_motor_efficiency_map(
    motor_map,
    data_csv,
    plot_png=None,
    rpm_min=500.0,
    rpm_max=None,
    torque_min_N_m=1.0,
    torque_max_N_m=None,
    n_rpm=180,
    n_torque=160,
):
    limits = electric_machine_limits(motor_map)
    rpm_max = (
        float(limits["speed_limit_rad_s"] * 60.0 / (2.0 * np.pi))
        if rpm_max is None
        else float(rpm_max)
    )
    torque_max_N_m = (
        float(limits["rated_torque_N_m"])
        if torque_max_N_m is None
        else float(torque_max_N_m)
    )
    rpm = np.linspace(rpm_min, rpm_max, n_rpm)
    torque = np.linspace(torque_min_N_m, torque_max_N_m, n_torque)
    rpm_grid, torque_grid = np.meshgrid(rpm, torque)
    efficiency_grid = electric_machine_efficiency(
        motor_map,
        speed_rad_s=rpm_to_rad_s(rpm_grid),
        torque_N_m=torque_grid,
    )

    data_csv = Path(data_csv)
    data_csv.parent.mkdir(parents=True, exist_ok=True)
    with data_csv.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=("rpm", "torque_N_m", "efficiency"),
        )
        writer.writeheader()
        for rpm_value, torque_value, efficiency in zip(
            np.ravel(rpm_grid),
            np.ravel(torque_grid),
            np.ravel(efficiency_grid),
        ):
            writer.writerow(
                {
                    "rpm": f"{float(rpm_value):.12g}",
                    "torque_N_m": f"{float(torque_value):.12g}",
                    "efficiency": f"{float(efficiency):.12g}",
                }
            )

    if plot_png is not None:
        import matplotlib.pyplot as plt

        plot_png = Path(plot_png)
        plot_png.parent.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(10.0, 6.5), constrained_layout=True)
        filled = ax.contourf(
            rpm_grid,
            torque_grid,
            efficiency_grid,
            levels=np.linspace(0.70, 0.97, 28),
            cmap="viridis",
        )
        contours = ax.contour(
            rpm_grid,
            torque_grid,
            efficiency_grid,
            levels=(0.80, 0.85, 0.90, 0.92, 0.94, 0.95, 0.96),
            colors="white",
            linewidths=0.9,
        )
        ax.clabel(contours, inline=True, fmt=lambda value: f"{100 * value:.0f}%", fontsize=8)
        ax.scatter(
            [float(motor_map.peak_speed_rad_s * 60.0 / (2.0 * np.pi))],
            [float(motor_map.peak_torque_N_m)],
            color="tab:red",
            label="Peak efficiency point",
            zorder=5,
        )
        ax.axvline(
            float(limits["speed_limit_rad_s"] * 60.0 / (2.0 * np.pi)),
            color="tab:orange",
            linestyle="--",
            linewidth=1.4,
            label="Speed limit",
        )
        ax.axhline(
            float(limits["rated_torque_N_m"]),
            color="tab:blue",
            linestyle="--",
            linewidth=1.4,
            label="Rated torque",
        )
        colorbar = fig.colorbar(filled, ax=ax)
        colorbar.set_label("Efficiency")
        ax.set_title("Sized motor efficiency map")
        ax.set_xlabel("Speed, rpm")
        ax.set_ylabel("Torque, N-m")
        ax.grid(True, alpha=0.25)
        ax.legend(loc="best")
        fig.savefig(plot_png, dpi=180)
        plt.close(fig)

    return {
        "motor_efficiency_map_csv": str(data_csv),
        "motor_efficiency_map_png": str(plot_png) if plot_png is not None else "",
    }


powertrain_output_columns = {
    "electric_power_hp",
    "fan1_motor_terminal_hp",
    "fan2_motor_terminal_hp",
    "fan1_generator_electric_hp",
    "fan2_generator_electric_hp",
    "fan1_generator_shaft_hp",
    "fan2_generator_shaft_hp",
    "per_engine_motor_terminal_hp",
    "per_engine_generator_electric_hp",
    "per_engine_generator_shaft_hp",
    "aircraft_motor_terminal_hp",
    "aircraft_generator_electric_hp",
    "aircraft_generator_shaft_hp",
    "generator_shaft_power_hp",
    "fan1_motor_efficiency",
    "fan2_motor_efficiency",
    "fan1_generator_efficiency",
    "fan2_generator_efficiency",
    "distribution_efficiency",
    "inverter_efficiency",
    "cable_efficiency",
    "per_turbine_shaft_power_hp",
    "per_turbine_shaft_power_kw",
    "turboshaft_required_throttle",
    "turboshaft_power_available_hp",
    "turboshaft_fuel_flow_lb_hr_per_turbine",
    "aircraft_turboshaft_fuel_flow_lb_hr",
    "aircraft_turboshaft_fuel_flow_lbm_s",
    "generator_turbine_freestream_mach",
    "generator_turbine_deck_mach",
    "generator_turbine_engine_face_mach",
    "generator_turbine_pressure_recovery",
    "generator_turbine_freestream_total_temperature_R",
    "generator_turbine_freestream_total_pressure_psf",
    "generator_turbine_total_temperature_R",
    "generator_turbine_total_pressure_psf",
    "generator_turbine_static_temperature_R",
    "generator_turbine_static_pressure_psf",
    "generator_turbine_density_slug_ft3",
    "generator_turbine_velocity_ft_s",
    "generator_turbine_corrected_airflow_lbm_s_per_turbine",
    "generator_turbine_airflow_lbm_s_per_turbine",
    "generator_turbine_engine_face_area_in2_per_turbine",
    "aircraft_generator_turbine_engine_face_area_in2",
    "turboshaft_status",
}


def build_duality_powertrain_deck(
    deck_csv,
    output_csv=None,
    *,
    sizing_summary_json=None,
    sizing_summary_csv=None,
    motor_efficiency_map_csv=None,
    motor_efficiency_map_png=None,
    number_duality_engines=2,
    motors_per_duality_engine=2,
    generators_per_motor=1,
    number_turbines=2,
    turboshaft_deck=None,
    generator_turbine_engine_face_mach=0.35,
    generator_turbine_subsonic_pressure_recovery=0.98,
    generator_turbine_min_pressure_recovery=0.35,
    generator_turbine_supersonic_recovery_coefficient=0.075,
    generator_turbine_supersonic_recovery_exponent=1.35,
    motor_peak_speed_rpm=6000.0,
    generator_peak_speed_rpm=12000.0,
    motor_peak_efficiency=0.96,
    generator_peak_efficiency=0.965,
    motor_parasite_loss_ratio=0.35,
    generator_parasite_loss_ratio=0.30,
    inverter_efficiency=0.97,
    cable_efficiency=0.99,
    sizing_margin=1.0,
):
    """Size electric machines from a Duality deck and write a powertrain deck.

    The deck is expected to use the imperial CSV convention generated by
    `write_pycycle_engine_deck()`. Fan shaft powers from pyCycle size the motor;
    the recomputed motor terminal power plus distribution losses then size the
    generator map.
    """
    deck_csv = Path(deck_csv)
    output_csv = (
        deck_csv.with_name(f"{deck_csv.stem}_powertrain.csv")
        if output_csv is None
        else Path(output_csv)
    )
    with deck_csv.open(newline="") as stream:
        reader = csv.DictReader(stream)
        rows = list(reader)
        fieldnames = [
            name
            for name in list(reader.fieldnames or [])
            if name not in powertrain_output_columns
        ]

    if not rows:
        raise ValueError(f"No rows found in {deck_csv}.")
    number_motors = int(number_duality_engines * motors_per_duality_engine)
    number_generators = int(number_motors * generators_per_motor)
    if number_motors <= 0 or number_generators <= 0 or number_turbines <= 0:
        raise ValueError("number_motors, number_generators, and number_turbines must be positive.")

    max_motor_shaft_hp = max(
        max(
            _row_float(row, "fan1_shaft_power_hp"),
            _row_float(row, "fan2_shaft_power_hp"),
        )
        for row in rows
    )
    motor_map = electric_machine_for_required_power(
        required_power_W=max_motor_shaft_hp * hp_to_W,
        peak_speed_rpm=motor_peak_speed_rpm,
        peak_efficiency=motor_peak_efficiency,
        parasite_loss_ratio=motor_parasite_loss_ratio,
        sizing_margin=sizing_margin,
    )

    generator_electric_requests_W = []
    for row in rows:
        fan1_speed_rpm = _row_float(row, "fan1_speed_rpm", motor_peak_speed_rpm)
        fan2_speed_rpm = _row_float(row, "fan2_speed_rpm", motor_peak_speed_rpm)
        if fan1_speed_rpm <= 0.0:
            fan1_speed_rpm = motor_peak_speed_rpm
        if fan2_speed_rpm <= 0.0:
            fan2_speed_rpm = motor_peak_speed_rpm
        fan1_motor_terminal_W = motor_electrical_power_W(
            motor_map,
            _row_float(row, "fan1_shaft_power_hp") * hp_to_W,
            rpm_to_rad_s(fan1_speed_rpm),
        )
        fan2_motor_terminal_W = motor_electrical_power_W(
            motor_map,
            _row_float(row, "fan2_shaft_power_hp") * hp_to_W,
            rpm_to_rad_s(fan2_speed_rpm),
        )
        distribution_efficiency = inverter_efficiency * cable_efficiency
        generator_electric_requests_W.extend(
            [
                fan1_motor_terminal_W / distribution_efficiency,
                fan2_motor_terminal_W / distribution_efficiency,
            ]
        )

    max_generator_electric_W = max(generator_electric_requests_W)
    generator_map = electric_machine_for_required_power(
        required_power_W=max_generator_electric_W,
        peak_speed_rpm=generator_peak_speed_rpm,
        peak_efficiency=generator_peak_efficiency,
        parasite_loss_ratio=generator_parasite_loss_ratio,
        sizing_margin=sizing_margin,
    )

    output_rows = []
    max_aircraft_generator_shaft_hp = 0.0
    max_aircraft_generator_turbine_inlet_area_in2 = 0.0
    max_unit_generator_shaft_hp = 0.0
    max_generator_row = None
    for row in rows:
        fan1_speed_rpm = _row_float(row, "fan1_speed_rpm", motor_peak_speed_rpm)
        fan2_speed_rpm = _row_float(row, "fan2_speed_rpm", motor_peak_speed_rpm)
        if fan1_speed_rpm <= 0.0:
            fan1_speed_rpm = motor_peak_speed_rpm
        if fan2_speed_rpm <= 0.0:
            fan2_speed_rpm = motor_peak_speed_rpm

        powertrain = duality_electric_powertrain(
            fan1_shaft_power_W=_row_float(row, "fan1_shaft_power_hp") * hp_to_W,
            fan2_shaft_power_W=_row_float(row, "fan2_shaft_power_hp") * hp_to_W,
            fan1_speed_rpm=fan1_speed_rpm,
            fan2_speed_rpm=fan2_speed_rpm,
            motor_map=motor_map,
            generator_map=generator_map,
            generator_speed_rpm=generator_peak_speed_rpm,
            number_engines=number_duality_engines,
            inverter_efficiency=inverter_efficiency,
            cable_efficiency=cable_efficiency,
        )
        fan1_shaft_W = _row_float(row, "fan1_shaft_power_hp") * hp_to_W
        fan2_shaft_W = _row_float(row, "fan2_shaft_power_hp") * hp_to_W
        fan1_motor_efficiency = (
            fan1_shaft_W / powertrain["fan1_motor_terminal_W"]
            if powertrain["fan1_motor_terminal_W"] > 0.0
            else 0.0
        )
        fan2_motor_efficiency = (
            fan2_shaft_W / powertrain["fan2_motor_terminal_W"]
            if powertrain["fan2_motor_terminal_W"] > 0.0
            else 0.0
        )
        fan1_generator_efficiency = (
            powertrain["fan1_generator_electric_W"] / powertrain["fan1_generator_shaft_W"]
            if powertrain["fan1_generator_shaft_W"] > 0.0
            else 0.0
        )
        fan2_generator_efficiency = (
            powertrain["fan2_generator_electric_W"] / powertrain["fan2_generator_shaft_W"]
            if powertrain["fan2_generator_shaft_W"] > 0.0
            else 0.0
        )

        updates = {
            "fan1_motor_terminal_hp": powertrain["fan1_motor_terminal_W"] * W_to_hp,
            "fan2_motor_terminal_hp": powertrain["fan2_motor_terminal_W"] * W_to_hp,
            "fan1_generator_electric_hp": powertrain["fan1_generator_electric_W"] * W_to_hp,
            "fan2_generator_electric_hp": powertrain["fan2_generator_electric_W"] * W_to_hp,
            "fan1_generator_shaft_hp": powertrain["fan1_generator_shaft_W"] * W_to_hp,
            "fan2_generator_shaft_hp": powertrain["fan2_generator_shaft_W"] * W_to_hp,
            "per_engine_motor_terminal_hp": powertrain["per_engine_motor_terminal_W"] * W_to_hp,
            "per_engine_generator_electric_hp": powertrain["per_engine_generator_electric_W"] * W_to_hp,
            "per_engine_generator_shaft_hp": powertrain["per_engine_generator_shaft_W"] * W_to_hp,
            "aircraft_motor_terminal_hp": powertrain["aircraft_motor_terminal_W"] * W_to_hp,
            "aircraft_generator_electric_hp": powertrain["aircraft_generator_electric_W"] * W_to_hp,
            "aircraft_generator_shaft_hp": powertrain["aircraft_generator_shaft_W"] * W_to_hp,
            "generator_shaft_power_hp": powertrain["per_engine_generator_shaft_W"] * W_to_hp,
            "fan1_motor_efficiency": fan1_motor_efficiency,
            "fan2_motor_efficiency": fan2_motor_efficiency,
            "fan1_generator_efficiency": fan1_generator_efficiency,
            "fan2_generator_efficiency": fan2_generator_efficiency,
            "distribution_efficiency": powertrain["distribution_efficiency"],
            "inverter_efficiency": powertrain["inverter_efficiency"],
            "cable_efficiency": powertrain["cable_efficiency"],
        }
        per_turbine_shaft_hp = updates["aircraft_generator_shaft_hp"] / number_turbines
        inlet_state = generator_turbine_inlet_state(
            mach=_row_float(row, "mach"),
            altitude_ft=_row_float(row, "altitude_ft"),
            engine_face_mach=generator_turbine_engine_face_mach,
            subsonic_pressure_recovery=generator_turbine_subsonic_pressure_recovery,
            min_pressure_recovery=generator_turbine_min_pressure_recovery,
            supersonic_recovery_coefficient=generator_turbine_supersonic_recovery_coefficient,
            supersonic_recovery_exponent=generator_turbine_supersonic_recovery_exponent,
        )
        turboshaft_corrected_airflow_lbm_s = (
            corrected_airflow_lbm_s(turboshaft_deck)
            if turboshaft_deck is not None
            else 0.0
        )
        inlet_sizing = generator_turbine_inlet_sizing(
            inlet_state,
            turboshaft_corrected_airflow_lbm_s,
        )
        turbine_updates = {
            "per_turbine_shaft_power_hp": per_turbine_shaft_hp,
            "per_turbine_shaft_power_kw": per_turbine_shaft_hp * hp_to_W / 1000.0,
            "turboshaft_required_throttle": "",
            "turboshaft_power_available_hp": "",
            "turboshaft_fuel_flow_lb_hr_per_turbine": "",
            "aircraft_turboshaft_fuel_flow_lb_hr": "",
            "aircraft_turboshaft_fuel_flow_lbm_s": "",
            "generator_turbine_freestream_mach": _row_float(row, "mach"),
            "generator_turbine_deck_mach": min(
                generator_turbine_engine_face_mach,
                turboshaft_deck.interpolator.bounds.mach_max if turboshaft_deck is not None else generator_turbine_engine_face_mach,
            ),
            "generator_turbine_engine_face_mach": inlet_state["engine_face_mach"],
            "generator_turbine_pressure_recovery": inlet_state["pressure_recovery"],
            "generator_turbine_freestream_total_temperature_R": inlet_state["freestream_total_temperature_R"],
            "generator_turbine_freestream_total_pressure_psf": inlet_state["freestream_total_pressure_psf"],
            "generator_turbine_total_temperature_R": inlet_state["total_temperature_R"],
            "generator_turbine_total_pressure_psf": inlet_state["total_pressure_psf"],
            "generator_turbine_static_temperature_R": inlet_state["static_temperature_R"],
            "generator_turbine_static_pressure_psf": inlet_state["static_pressure_psf"],
            "generator_turbine_density_slug_ft3": inlet_state["density_slug_ft3"],
            "generator_turbine_velocity_ft_s": inlet_state["velocity_ft_s"],
            "generator_turbine_corrected_airflow_lbm_s_per_turbine": inlet_sizing["corrected_airflow_lbm_s"],
            "generator_turbine_airflow_lbm_s_per_turbine": inlet_sizing["airflow_lbm_s"],
            "generator_turbine_engine_face_area_in2_per_turbine": inlet_sizing["engine_face_area_in2"],
            "aircraft_generator_turbine_engine_face_area_in2": inlet_sizing["engine_face_area_in2"] * number_turbines,
            "turboshaft_status": "not_evaluated",
        }
        if turboshaft_deck is not None and per_turbine_shaft_hp > 0.0:
            try:
                turbine_result = turboshaft_deck.required_throttle_for_power(
                    mach=turbine_updates["generator_turbine_deck_mach"],
                    altitude_ft=_row_float(row, "altitude_ft"),
                    required_power_kw=per_turbine_shaft_hp
                    * hp_to_W
                    / (1000.0 * max(inlet_state["pressure_recovery"], 1.0e-6)),
                )
                fuel_flow_lb_hr_per_turbine = turbine_result["fuel_flow_lb_hr"]
                turbine_updates.update(
                    {
                        "turboshaft_required_throttle": turbine_result["required_throttle"],
                        "turboshaft_power_available_hp": turbine_result["shaft_power_hp"],
                        "turboshaft_fuel_flow_lb_hr_per_turbine": fuel_flow_lb_hr_per_turbine,
                        "aircraft_turboshaft_fuel_flow_lb_hr": fuel_flow_lb_hr_per_turbine * number_turbines,
                        "aircraft_turboshaft_fuel_flow_lbm_s": fuel_flow_lb_hr_per_turbine * number_turbines / 3600.0,
                        "turboshaft_status": "ok",
                    }
                )
            except ValueError as exc:
                turbine_updates["turboshaft_status"] = str(exc)
        elif per_turbine_shaft_hp <= 0.0:
            turbine_updates["turboshaft_status"] = "no_generator_load"

        output_row = {
            key: value
            for key, value in row.items()
            if key not in powertrain_output_columns
        }
        output_row.update({key: f"{float(value):.12g}" for key, value in updates.items()})
        for key, value in turbine_updates.items():
            output_row[key] = f"{float(value):.12g}" if isinstance(value, (int, float)) else value
        output_rows.append(output_row)
        aircraft_generator_shaft_hp = updates["aircraft_generator_shaft_hp"]
        max_unit_generator_shaft_hp = max(
            max_unit_generator_shaft_hp,
            updates["fan1_generator_shaft_hp"],
            updates["fan2_generator_shaft_hp"],
        )
        max_aircraft_generator_turbine_inlet_area_in2 = max(
            max_aircraft_generator_turbine_inlet_area_in2,
            turbine_updates["aircraft_generator_turbine_engine_face_area_in2"],
        )
        if aircraft_generator_shaft_hp > max_aircraft_generator_shaft_hp:
            max_aircraft_generator_shaft_hp = aircraft_generator_shaft_hp
            max_generator_row = row

    for row in output_rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(output_rows)

    motor_summary = _machine_summary(motor_map, "motor")
    generator_summary = _machine_summary(generator_map, "generator")
    flowpath_summary = _sizing_flowpath_summary(max_generator_row or {})
    summary = {
        "deck_csv": str(deck_csv),
        "output_csv": str(output_csv),
        "number_duality_engines": int(number_duality_engines),
        "motors_per_duality_engine": int(motors_per_duality_engine),
        "number_motors": number_motors,
        "generators_per_motor": int(generators_per_motor),
        "number_generators": number_generators,
        "number_turbines": int(number_turbines),
        "unit_motor_shaft_power_hp": float(max_motor_shaft_hp),
        "unit_motor_shaft_power_kw": float(max_motor_shaft_hp * hp_to_W / 1000.0),
        "unit_motor_rated_power_hp": float(motor_summary["motor_rated_power_hp"]),
        "unit_motor_rated_power_kw": float(motor_summary["motor_rated_power_kw"]),
        "unit_generator_electric_power_hp": float(max_generator_electric_W * W_to_hp),
        "unit_generator_electric_power_kw": float(max_generator_electric_W / 1000.0),
        "unit_generator_shaft_power_hp": float(max_unit_generator_shaft_hp),
        "unit_generator_shaft_power_kw": float(max_unit_generator_shaft_hp * hp_to_W / 1000.0),
        "unit_turbine_shaft_power_hp": float(max_aircraft_generator_shaft_hp / number_turbines),
        "unit_turbine_shaft_power_kw": float(max_aircraft_generator_shaft_hp * hp_to_W / (1000.0 * number_turbines)),
        "aircraft_generator_shaft_power_hp": float(max_aircraft_generator_shaft_hp),
        "aircraft_generator_shaft_power_kw": float(max_aircraft_generator_shaft_hp * hp_to_W / 1000.0),
        "aircraft_generator_turbine_inlet_area_in2": float(max_aircraft_generator_turbine_inlet_area_in2),
        "sizing_mach": float(_row_float(max_generator_row or {}, "mach")),
        "sizing_altitude_ft": float(_row_float(max_generator_row or {}, "altitude_ft")),
        **flowpath_summary,
        "motor_peak_speed_rpm": motor_summary["motor_peak_speed_rpm"],
        "motor_peak_torque_N_m": motor_summary["motor_peak_torque_N_m"],
        "motor_peak_efficiency": motor_summary["motor_peak_efficiency"],
        "motor_rated_torque_N_m": motor_summary["motor_rated_torque_N_m"],
        "generator_peak_speed_rpm": generator_summary["generator_peak_speed_rpm"],
        "generator_peak_torque_N_m": generator_summary["generator_peak_torque_N_m"],
        "generator_peak_efficiency": generator_summary["generator_peak_efficiency"],
        "generator_rated_torque_N_m": generator_summary["generator_rated_torque_N_m"],
    }
    if sizing_summary_json is not None:
        write_sizing_summary(summary, sizing_summary_json, sizing_summary_csv)
    if motor_efficiency_map_csv is not None:
        summary.update(
            write_motor_efficiency_map(
                motor_map,
                motor_efficiency_map_csv,
                plot_png=motor_efficiency_map_png,
            )
        )
        if sizing_summary_json is not None:
            write_sizing_summary(summary, sizing_summary_json, sizing_summary_csv)
    return summary


def update_duality_deck_electric_sizing(deck_csv, output_csv=None, **kwargs):
    return build_duality_powertrain_deck(deck_csv, output_csv, **kwargs)


def duality_electric_powertrain(
    fan1_shaft_power_W,
    fan2_shaft_power_W,
    fan1_speed_rpm,
    fan2_speed_rpm,
    motor_map,
    generator_map,
    generator_speed_rpm=12000.0,
    number_engines=2,
    inverter_efficiency=0.97,
    cable_efficiency=0.99,
):
    """Return power flows for the Duality two-engine electric fan system.

    Inputs `fan*_shaft_power_W` are per engine. Each engine has two fan motors
    and two generators. The total aircraft has `number_engines` engines.
    """
    fan1_speed_rad_s = fan1_speed_rpm * 2.0 * np.pi / 60.0
    fan2_speed_rad_s = fan2_speed_rpm * 2.0 * np.pi / 60.0
    generator_speed_rad_s = generator_speed_rpm * 2.0 * np.pi / 60.0

    fan1_motor_terminal_W = motor_electrical_power_W(motor_map, fan1_shaft_power_W, fan1_speed_rad_s)
    fan2_motor_terminal_W = motor_electrical_power_W(motor_map, fan2_shaft_power_W, fan2_speed_rad_s)
    distribution_efficiency = inverter_efficiency * cable_efficiency
    fan1_generator_electric_W = fan1_motor_terminal_W / distribution_efficiency
    fan2_generator_electric_W = fan2_motor_terminal_W / distribution_efficiency
    fan1_generator_shaft_W = generator_shaft_power_W(generator_map, fan1_generator_electric_W, generator_speed_rad_s)
    fan2_generator_shaft_W = generator_shaft_power_W(generator_map, fan2_generator_electric_W, generator_speed_rad_s)

    per_engine_motor_terminal_W = fan1_motor_terminal_W + fan2_motor_terminal_W
    per_engine_generator_electric_W = fan1_generator_electric_W + fan2_generator_electric_W
    per_engine_generator_shaft_W = fan1_generator_shaft_W + fan2_generator_shaft_W
    return {
        "fan1_motor_terminal_W": fan1_motor_terminal_W,
        "fan2_motor_terminal_W": fan2_motor_terminal_W,
        "fan1_generator_electric_W": fan1_generator_electric_W,
        "fan2_generator_electric_W": fan2_generator_electric_W,
        "fan1_generator_shaft_W": fan1_generator_shaft_W,
        "fan2_generator_shaft_W": fan2_generator_shaft_W,
        "per_engine_motor_terminal_W": per_engine_motor_terminal_W,
        "per_engine_generator_electric_W": per_engine_generator_electric_W,
        "per_engine_generator_shaft_W": per_engine_generator_shaft_W,
        "aircraft_motor_terminal_W": number_engines * per_engine_motor_terminal_W,
        "aircraft_generator_electric_W": number_engines * per_engine_generator_electric_W,
        "aircraft_generator_shaft_W": number_engines * per_engine_generator_shaft_W,
        "distribution_efficiency": distribution_efficiency,
        "inverter_efficiency": inverter_efficiency,
        "cable_efficiency": cable_efficiency,
        "number_motors": 2 * number_engines,
        "number_generators": 2 * number_engines,
    }


def main():
    # Edit run options here.
    save_plot = "electric_machine_efficiency_map.png"
    show_plot = False
    rpm_min = 500.0
    rpm_max = 14000.0
    torque_min_N_m = 1.0
    torque_max_N_m = 500.0
    n_rpm = 180
    n_torque = 160

    motor_map = ElectricMachineMap(
        peak_speed_rad_s=6000.0 * 2.0 * np.pi / 60.0,
        peak_torque_N_m=250.0,
        peak_efficiency=0.96,
        parasite_loss_ratio=0.35,
    )
    generator_map = ElectricMachineMap(
        peak_speed_rad_s=12000.0 * 2.0 * np.pi / 60.0,
        peak_torque_N_m=150.0,
        peak_efficiency=0.965,
        parasite_loss_ratio=0.30,
    )
    result = duality_electric_powertrain(
        fan1_shaft_power_W=120000.0,
        fan2_shaft_power_W=90000.0,
        fan1_speed_rpm=6000.0,
        fan2_speed_rpm=6000.0,
        motor_map=motor_map,
        generator_map=generator_map,
        generator_speed_rpm=12000.0,
        number_engines=2,
    )

    print(" electric machine powertrain")
    print(f"Per-engine motor terminal power: {result['per_engine_motor_terminal_W']:.3f} W")
    print(f"Per-engine generator electric power: {result['per_engine_generator_electric_W']:.3f} W")
    print(f"Per-engine generator shaft power: {result['per_engine_generator_shaft_W']:.3f} W")
    print(f"Aircraft motor terminal power: {result['aircraft_motor_terminal_W']:.3f} W")
    print(f"Aircraft generator electric power: {result['aircraft_generator_electric_W']:.3f} W")
    print(f"Aircraft gas-turbine generator shaft power: {result['aircraft_generator_shaft_W']:.3f} W")

    import matplotlib.pyplot as plt

    rpm = np.linspace(rpm_min, rpm_max, n_rpm)
    torque = np.linspace(torque_min_N_m, torque_max_N_m, n_torque)
    rpm_grid, torque_grid = np.meshgrid(rpm, torque)
    speed_grid_rad_s = rpm_grid * 2.0 * np.pi / 60.0
    efficiency_grid = electric_machine_efficiency(
        motor_map,
        speed_rad_s=speed_grid_rad_s,
        torque_N_m=torque_grid,
    )

    fig, ax = plt.subplots(figsize=(10.0, 6.5), constrained_layout=True)
    filled = ax.contourf(
        rpm_grid,
        torque_grid,
        efficiency_grid,
        levels=np.linspace(0.70, 0.97, 28),
        cmap="viridis",
    )
    contours = ax.contour(
        rpm_grid,
        torque_grid,
        efficiency_grid,
        levels=(0.80, 0.85, 0.90, 0.92, 0.94, 0.95, 0.96),
        colors="white",
        linewidths=0.9,
    )
    ax.clabel(contours, inline=True, fmt=lambda value: f"{100 * value:.0f}%", fontsize=8)
    ax.scatter(
        [motor_map.peak_speed_rad_s * 60.0 / (2.0 * np.pi)],
        [motor_map.peak_torque_N_m],
        color="tab:red",
        label="Peak efficiency point",
        zorder=5,
    )
    limits = electric_machine_limits(motor_map)
    ax.axvline(
        limits["speed_limit_rad_s"] * 60.0 / (2.0 * np.pi),
        color="tab:orange",
        linestyle="--",
        linewidth=1.4,
        label="Speed limit",
    )
    ax.axhline(
        limits["rated_torque_N_m"],
        color="tab:blue",
        linestyle="--",
        linewidth=1.4,
        label="Rated torque",
    )
    colorbar = fig.colorbar(filled, ax=ax)
    colorbar.set_label("Efficiency")
    ax.set_title("Electric machine efficiency map")
    ax.set_xlabel("Speed, rpm")
    ax.set_ylabel("Torque, N-m")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")

    if save_plot is not None:
        fig.savefig(save_plot, dpi=180)
        print(f"Efficiency map plot: {save_plot}")
    if show_plot:
        plt.show()
    else:
        plt.close(fig)


if __name__ == "__main__":
    main()
