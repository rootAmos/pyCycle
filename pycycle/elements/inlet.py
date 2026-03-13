"""
Inlet element for gas turbine cycle analysis.

Computes inlet total pressure rise (ram recovery), ram drag, and exit flow state.
Flow enters at freestream conditions (Fl_I) and exits at engine face (Fl_O) with
higher total pressure and known static state (design: specified MN; off-design: fixed area).
"""

import openmdao.api as om

from pycycle.constants import g_c  # gravitational constant for lbm<->lbf: 32.174 lbm*ft/(lbf*s^2)

from pycycle.thermo.cea import species_data
from pycycle.thermo.thermo import Thermo

from pycycle.flow_in import FlowIn
from pycycle.passthrough import PassThrough
from pycycle.element_base import Element


# ---------------------------------------------------------------------------
# MilSpecRecovery: Ram pressure recovery vs. Mach number (MIL-E-5007B style)
# ---------------------------------------------------------------------------
class MilSpecRecovery(om.ExplicitComponent):
    """
    Computes total-pressure recovery factor (eta_ram) from flight Mach number.

    Ram recovery = Pt_capture / Pt_freestream. Subsonic inlets can achieve
    ~1.0; supersonic shocks cause losses; this implements an empirical curve.
    """

    def setup(self):
        # MN: flight Mach number (dimensionless). Determines which recovery formula is used.
        self.add_input('MN', val=0.5, units=None, desc='Flight Mach Number')
        # ram_recovery_base: ideal or design recovery (e.g. 0.995). Used as-is for subsonic.
        self.add_input('ram_recovery_base', units=None, desc='Base Inlet Ram Recovery')

        # ram_recovery: effective total-pressure recovery (0 to 1). Output used in Pt_out = Pt_in * ram_recovery.
        self.add_output('ram_recovery', val=1.0, units=None, desc='Mil Spec Ram Recovery')

        self.declare_partials('ram_recovery', ['ram_recovery_base', 'MN'])

    def compute(self, inputs, outputs):
        MN = inputs['MN']
        ram_recovery_base = inputs['ram_recovery_base']

        # Subsonic (MN < 1): no shock losses; recovery = base value.
        if MN < 1.0:
            outputs['ram_recovery'] = ram_recovery_base

        # Supersonic (1 <= MN < 5): empirical loss factor from MIL-E-5007B-type correlation.
        # Formula: eta = eta_base * (1 - 0.075*(MN-1)^1.35). Loss increases with Mach.
        elif MN >= 1.0 and MN < 5.0:
            outputs['ram_recovery'] = ram_recovery_base * (1 - (0.075 * ((MN - 1) ** 1.35)))

        # Hypersonic (MN >= 5): recovery dominated by strong shocks; curve 800/(MN^4 + 935).
        elif MN >= 5.0:
            outputs['ram_recovery'] = 800 / (MN ** 4 + 935)

    def compute_partials(self, inputs, J):
        """Jacobian of ram_recovery w.r.t. MN and ram_recovery_base for Newton/optimization."""
        MN = inputs['MN']
        ram_recovery_base = inputs['ram_recovery_base']

        if MN < 1.0:
            J['ram_recovery', 'ram_recovery_base'] = 1
            J['ram_recovery', 'MN'] = 0

        elif MN >= 1.0 and MN < 5.0:
            # d(eta)/d(eta_base) = (1 - 0.075*(MN-1)^1.35)
            J['ram_recovery', 'ram_recovery_base'] = 1 - (0.075 * ((MN - 1) ** 1.35))
            # d(eta)/d(MN) = eta_base * (-0.075*1.35*(MN-1)^0.35) = -0.10125*eta_base*(MN-1)^0.35
            J['ram_recovery', 'MN'] = -0.10125 * ram_recovery_base * ((MN - 1) ** 0.35)

        elif MN >= 5.0:
            J['ram_recovery', 'ram_recovery_base'] = 0
            # d/d(MN) of 800/(MN^4+935) = -800*4*MN^3/(MN^4+935)^2
            J['ram_recovery', 'MN'] = -(3200 * MN ** 3) / ((MN ** 4 + 935) ** 2)

# ---------------------------------------------------------------------------
# Calcs: Inlet total pressure and ram drag
# ---------------------------------------------------------------------------
class Calcs(om.ExplicitComponent):
    """
    Inlet engineering calculations: exit total pressure and ram drag.

    Equations:
      Pt_out = Pt_in * ram_recovery   (total pressure after inlet)
      F_ram  = W * V / g_c            (ram drag force in lbf; momentum of captured air)
    """

    def setup(self):
        # Pt_in: total pressure at inlet entrance (freestream total), lbf/in^2 (psi).
        self.add_input('Pt_in', val=5.0, units='lbf/inch**2', desc='Entrance total pressure')
        # ram_recovery: total-pressure recovery factor (0--1) from MilSpecRecovery or user.
        self.add_input('ram_recovery', val=1.0, desc='Ram recovery')
        # V_in: freestream velocity (ft/s). Used for ram drag F_ram = W*V/g_c.
        self.add_input('V_in', val=0.0, units='ft/s', desc='Entrance velocity')
        # W_in: mass flow rate (lbm/s) at inlet. From cycle balance (e.g. balance.W).

        self.add_input('W_in', val=100.0, units='lbm/s', desc='Entrance flow rate')

        # Pt_out: total pressure at engine face (psi). Pt_out = Pt_in * ram_recovery.
        self.add_output('Pt_out', val=14.696, units='lbf/inch**2', desc='Exit total pressure')
        # F_ram: ram drag (lbf). Force required to decelerate captured air; opposes thrust.
        self.add_output('F_ram', val=1.0, units='lbf', desc='Ram drag')

        self.declare_partials('Pt_out', ['Pt_in', 'ram_recovery'])
        self.declare_partials('F_ram', ['V_in', 'W_in'])

    def compute(self, inputs, outputs):
        # Total pressure at exit: entrance total pressure times recovery factor.
        outputs['Pt_out'] = inputs['Pt_in'] * inputs['ram_recovery']
        # Ram drag: momentum flux of incoming air. F = m_dot * V; g_c converts lbm*ft/s^2 to lbf.
        outputs['F_ram'] = inputs['W_in'] * inputs['V_in'] / g_c

    def compute_partials(self, inputs, J):
        # d(Pt_out)/d(Pt_in)=ram_recovery, d(Pt_out)/d(ram_recovery)=Pt_in
        J['Pt_out', 'Pt_in'] = inputs['ram_recovery']
        J['Pt_out', 'ram_recovery'] = inputs['Pt_in']
        # d(F_ram)/d(V_in)=W/g_c, d(F_ram)/d(W_in)=V/g_c
        J['F_ram', 'V_in'] = inputs['W_in'] / g_c
        J['F_ram', 'W_in'] = inputs['V_in'] / g_c


# ---------------------------------------------------------------------------
# Inlet: Main element — Fl_I (freestream) -> Fl_O (engine face)
# ---------------------------------------------------------------------------
class Inlet(Element):
    """
    Inlet element: ram recovery, ram drag, and exit flow state.

    Flow stations:
      Fl_I  = flow at entrance (from FlowStart: Pt, Tt, composition, W, static V, etc.)
      Fl_O  = flow at exit (engine face): same composition and Tt; Pt from ram recovery; statics from MN or area

    Design: exit static state is determined by specified MN (Mach number at engine face).
    Off-design: exit static state is determined by fixed area (from design); mass flow varies.
    """

    def initialize(self):
        # If True, compute Fl_O static properties (Ps, Ts, area, etc.). If False, only pass W through.
        self.options.declare('statics', default=True,
                             desc='If True, calculate static properties.')

        # Design-to-off-design: design exit throat area becomes the fixed 'area' input for off-design.
        self.default_des_od_conns = [
            ('Fl_O:stat:area', 'area'),
        ]

        super().initialize()

    def pyc_setup_output_ports(self):
        # Declare that Fl_O has the same flow structure as Fl_I; actual values set in setup.
        self.copy_flow('Fl_I', 'Fl_O')

    def setup(self):
        # Thermo package and design/statics flags (set by Cycle or problem).
        thermo_method = self.options['thermo_method']
        thermo_data = self.options['thermo_data']
        statics = self.options['statics']
        design = self.options['design']

        # Composition (mole fractions) for Fl_I; inherited from upstream and used for Fl_O thermo.
        composition = self.Fl_I_data['Fl_I']

        # ---- Flow input port: declares all Fl_I total and static flow variables as inputs ----
        # FlowIn(Fl_I) provides: Fl_I:tot:P, Fl_I:tot:T, Fl_I:tot:composition, Fl_I:stat:W, Fl_I:stat:V, etc.
        # These are connected from upstream (e.g. FlowStart) via pyc_connect_flow.
        flow_in = FlowIn(fl_name='Fl_I')
        self.add_subsystem('flow_in', flow_in, promotes=['Fl_I:tot:*', 'Fl_I:stat:*'])

        # ---- Inlet calcs: Pt_out and F_ram ----
        # Promotes: Pt_in <- Fl_I:tot:P, V_in <- Fl_I:stat:V, W_in <- Fl_I:stat:W; ram_recovery from options/IVC.
        # Outputs: Pt_out (engine face total pressure), F_ram (ram drag, lbf).
        self.add_subsystem('calcs_inlet', Calcs(),
                           promotes_inputs=['ram_recovery', ('Pt_in', 'Fl_I:tot:P'),
                                            ('V_in', 'Fl_I:stat:V'), ('W_in', 'Fl_I:stat:W')],
                           promotes_outputs=['F_ram'])

        # ---- Exit total state (Fl_O:tot): T and composition from Fl_I; P from calcs_inlet.Pt_out ----
        # Thermo in 'total_TP' mode: given T, P, composition -> h, S, gamma, etc. for Fl_O.
        real_flow = Thermo(mode='total_TP', fl_name='Fl_O:tot',
                           method=thermo_method,
                           thermo_kwargs={'composition': composition,
                                          'spec': thermo_data})
        self.add_subsystem('real_flow', real_flow,
                           promotes_inputs=[('T', 'Fl_I:tot:T'), ('composition', 'Fl_I:tot:composition')],
                           promotes_outputs=['Fl_O:*'])

        # Exit total pressure = inlet calcs output (Pt_in * ram_recovery).
        self.connect("calcs_inlet.Pt_out", "real_flow.P")

        # ---- Exit static state (Fl_O:stat): Ps, Ts, area, V, MN, etc. ----
        if statics:
            if design:
                # Design: exit Mach number MN is specified (e.g. 0.6 at engine face).
                # Thermo static_MN: given total state (S, h, Pt, gamma), W, and MN -> static P, T, area, V, rho.
                out_stat = Thermo(mode='static_MN', fl_name='Fl_O:stat',
                                  method=thermo_method,
                                  thermo_kwargs={'composition': composition,
                                                 'spec': thermo_data})
                self.add_subsystem('out_stat', out_stat,
                                   promotes_inputs=[('composition', 'Fl_I:tot:composition'), ('W', 'Fl_I:stat:W'), 'MN'],
                                   promotes_outputs=['Fl_O:stat:*'])
                # Total state from real_flow: entropy S, enthalpy h; guess for Pt and gamma for iterative static solve.
                self.connect('Fl_O:tot:S', 'out_stat.S')
                self.connect('Fl_O:tot:h', 'out_stat.ht')
                self.connect('Fl_O:tot:P', 'out_stat.guess:Pt')
                self.connect('Fl_O:tot:gamma', 'out_stat.guess:gamt')

            else:
                # Off-design: exit area is fixed (from design via pyc_connect_des_od).
                # Thermo static_A: given total state, W, and area -> static P, T, MN, V, etc.
                out_stat = Thermo(mode='static_A', fl_name='Fl_O:stat',
                                  method=thermo_method,
                                  thermo_kwargs={'composition': composition,
                                                 'spec': thermo_data})
                prom_in = [('composition', 'Fl_I:tot:composition'),
                           ('W', 'Fl_I:stat:W'),
                           'area']  # area = design Fl_O:stat:area, passed from cycle
                prom_out = ['Fl_O:stat:*']
                self.add_subsystem('out_stat', out_stat, promotes_inputs=prom_in,
                                   promotes_outputs=prom_out)
                self.connect('Fl_O:tot:S', 'out_stat.S')
                self.connect('Fl_O:tot:h', 'out_stat.ht')
                self.connect('Fl_O:tot:P', 'out_stat.guess:Pt')
                self.connect('Fl_O:tot:gamma', 'out_stat.guess:gamt')

        else:
            # No statics: only pass mass flow through; Fl_O:stat:W = Fl_I:stat:W.
            self.add_subsystem('W_passthru', PassThrough('Fl_I:stat:W', 'Fl_O:stat:W', 0.0, units="lbm/s"),
                               promotes=['*'])

        super().setup()

    def pyc_setup_thermo(self, upstream):
        # Tell framework that Fl_O has same composition/species structure as upstream Fl_I.
        elements = self.options['elements']

        self.Fl_O_data = {
            'Fl_O': upstream['Fl_I']
        }


# ---------------------------------------------------------------------------
# Standalone test: FlowStart -> Inlet with IVC-provided freestream and options
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    """
    Run Inlet in isolation: FlowStart computes Fl_O from P, T, MN, W; that flow
    is connected to Inlet.Fl_I. Inlet then computes ram recovery, F_ram, and Fl_O.
    """
    import openmdao.api as om
    from pycycle.mp_cycle import Cycle
    from pycycle.thermo.cea.species_data import janaf
    from pycycle.elements.flow_start import FlowStart

    prob = om.Problem()
    model = prob.model = om.Group()

    # IndepVarComp: freestream and inlet inputs (no balance; direct values for testing).
    ivc = model.add_subsystem('ivc', om.IndepVarComp(), promotes_outputs=['*'])
    ivc.add_output('P', 17.0, units='psi')           # freestream static pressure
    ivc.add_output('T', 500.0, units='degR')          # freestream static temperature
    ivc.add_output('MN', 0.5)                        # freestream Mach number
    ivc.add_output('W', 100.0, units='lbm/s')         # mass flow (would come from balance in full cycle)
    ivc.add_output('inlet_MN', 0.6)                  # engine-face Mach (design exit MN for inlet)
    ivc.add_output('ram_recovery', 0.995)             # total-pressure recovery factor

    # Cycle group: FlowStart (P,T,MN,W -> Fl_O) then Inlet (Fl_I -> Pt_out, F_ram, Fl_O).
    cycle = model.add_subsystem('cycle', Cycle())
    cycle.options['thermo_method'] = 'CEA'
    cycle.options['thermo_data'] = janaf
    cycle.options['design'] = True
    cycle.add_subsystem('flow_start', FlowStart())
    cycle.add_subsystem('inlet', Inlet())
    cycle.pyc_connect_flow('flow_start.Fl_O', 'inlet.Fl_I')  # connects all flow vars; includes Fl_I:stat:V

    # Wire IVC to flow_start and inlet (do not connect Fl_I:stat:V again — it comes from flow_start).
    model.connect('P', 'cycle.flow_start.P')
    model.connect('T', 'cycle.flow_start.T')
    model.connect('MN', 'cycle.flow_start.MN')
    model.connect('W', 'cycle.flow_start.W')
    model.connect('inlet_MN', 'cycle.inlet.MN')
    model.connect('ram_recovery', 'cycle.inlet.ram_recovery')

    prob.setup(check=False)
    prob.run_model()

    print("--- Inlet standalone test ---")
    print("Fl_O: Pt (psi) =", prob.get_val('cycle.inlet.Fl_O:tot:P')[0])
    print("Fl_O: Tt (degR) =", prob.get_val('cycle.inlet.Fl_O:tot:T')[0])
    print("Fl_O: stat:W (lbm/s) =", prob.get_val('cycle.inlet.Fl_O:stat:W')[0])
    print("F_ram (lbf) =", prob.get_val('cycle.inlet.F_ram')[0])

