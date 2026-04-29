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

import aerosandbox.numpy as np


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
    return shaft_power_W / (shaft_power_W + loss_power_W)


def motor_electrical_power_W(machine: ElectricMachineMap, shaft_power_W, speed_rad_s):
    """Return electrical input power for a motor driving a fan."""
    torque_N_m = shaft_power_W / speed_rad_s
    loss_power_W = electric_machine_loss_power_W(machine, speed_rad_s, torque_N_m)
    return shaft_power_W + loss_power_W


def generator_shaft_power_W(machine: ElectricMachineMap, electrical_output_W, speed_rad_s, iterations=12):
    """Return shaft power required for a generator to supply electrical output."""
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


def duality_electric_powertrain(
    fan1_shaft_power_W,
    fan2_shaft_power_W,
    fan1_speed_rpm,
    fan2_speed_rpm,
    motor_map,
    generator_map,
    generator_speed_rpm=12000.0,
    number_engines=2,
):
    """Return power flows for the Duality two-engine electric fan system.

    Inputs `fan*_shaft_power_W` are per engine. Each engine has two fan motors
    and two generators. The total aircraft has `number_engines` engines.
    """
    fan1_speed_rad_s = fan1_speed_rpm * 2.0 * np.pi / 60.0
    fan2_speed_rad_s = fan2_speed_rpm * 2.0 * np.pi / 60.0
    generator_speed_rad_s = generator_speed_rpm * 2.0 * np.pi / 60.0

    fan1_motor_electric_W = motor_electrical_power_W(motor_map, fan1_shaft_power_W, fan1_speed_rad_s)
    fan2_motor_electric_W = motor_electrical_power_W(motor_map, fan2_shaft_power_W, fan2_speed_rad_s)
    fan1_generator_shaft_W = generator_shaft_power_W(generator_map, fan1_motor_electric_W, generator_speed_rad_s)
    fan2_generator_shaft_W = generator_shaft_power_W(generator_map, fan2_motor_electric_W, generator_speed_rad_s)

    per_engine_motor_electric_W = fan1_motor_electric_W + fan2_motor_electric_W
    per_engine_generator_shaft_W = fan1_generator_shaft_W + fan2_generator_shaft_W
    return {
        "fan1_motor_electric_W": fan1_motor_electric_W,
        "fan2_motor_electric_W": fan2_motor_electric_W,
        "per_engine_motor_electric_W": per_engine_motor_electric_W,
        "per_engine_generator_shaft_W": per_engine_generator_shaft_W,
        "aircraft_motor_electric_W": number_engines * per_engine_motor_electric_W,
        "aircraft_generator_shaft_W": number_engines * per_engine_generator_shaft_W,
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
    print(f"Per-engine motor electric power: {result['per_engine_motor_electric_W']:.3f} W")
    print(f"Per-engine generator shaft power: {result['per_engine_generator_shaft_W']:.3f} W")
    print(f"Aircraft motor electric power: {result['aircraft_motor_electric_W']:.3f} W")
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
