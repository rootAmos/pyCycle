"""
Tabular ThermoAdd: mix inflow with reactant or flow and output composition vector and mass-averaged enthalpy.

Composition is a vector of ratios (e.g. FAR and optionally others). In reactant mode,
  one reactant name (e.g. from TAB composition keys) and ratio(s) per mix port;
  in flow mode, composition vectors and W per mix port. Outputs composition_out (vector),
  Wout, mass_avg_h. Partial derivatives use complex step.
"""

import numpy as np

import openmdao.api as om

from pycycle.constants import AIR_JETA_TAB_SPEC, TAB_AIR_FUEL_COMPOSITION


class ThermoAdd(om.ExplicitComponent):
    """
    Tabular ThermoAdd: new composition_out and mass_avg_h from inflow + mix(es).

    In reactant mode all mix ports use the same reactant type (one key in composition).
    For multiple reactants use separate ThermoAdd instances.
    """

    def initialize(self):
        self.options.declare('spec', default=AIR_JETA_TAB_SPEC, recordable=False)
        self.options.declare('inflow_composition', default=None, 
                             desc='composition present in the inflow')

        self.options.declare('mix_mode', values=['reactant', 'flow'], default='reactant')

        self.options.declare('mix_composition', default=None, 
                             desc='name of the mixing reactant; must match one of the keys from the inflow composition dictionary', 
                             types=(dict, str, list, tuple), allow_none=True)
        self.options.declare('mix_names', default='mix', types=(str, list, tuple))

    def output_port_data(self):

        spec = self.options['spec']
        
        inflow_composition = self.options['inflow_composition']
        if inflow_composition is None: 
            inflow_composition = TAB_AIR_FUEL_COMPOSITION

        mix_mode = self.options['mix_mode']

        self.sorted_compo = sorted(inflow_composition.keys())

        if mix_mode == "reactant":
            reactant = self.options['mix_composition']
            self.idx_compo = self.sorted_compo.index(reactant)

    def setup(self):

        spec = self.options['spec']
        mix_mode = self.options['mix_mode']
        mix_names = self.options['mix_names']

        if isinstance(mix_names, str): # cast it to tuple 
            mix_names = (mix_names,)    
        self.mix_names = mix_names

        inflow_composition = self.options['inflow_composition']
        if inflow_composition is None: 
            inflow_composition = TAB_AIR_FUEL_COMPOSITION

        inflow_composition_vec = list(inflow_composition.values())

        self.output_port_data()

        # inputs
        self.add_input('Fl_I:stat:W', val=0.0, desc='weight flow', units='lbm/s')
        self.add_input('Fl_I:tot:h', val=0.0, desc='total enthalpy', units='Btu/lbm')
        self.add_input('Fl_I:tot:composition', val=inflow_composition_vec, 
                       desc='incoming flow composition')
        
        for name in mix_names: 
            self.add_input(f'{name}:h', val=0.0, units='Btu/lbm', desc="reactant enthalpy")

            if mix_mode == 'reactant': 
                self.add_input(f'{name}:ratio', val=0.0, desc='reactant to air mass ratio')
                self.add_output(f'{name}:W', shape=1, units="lbm/s", desc="mix input massflow")

            else: 
                self.add_input(f'{name}:composition', val=inflow_composition_vec, desc='mix flow composition' )
                self.add_input(f'{name}:W', shape=1, units="lbm/s", desc="mix input massflow")

        # outputs
        self.add_output('mass_avg_h', shape=1, units='Btu/lbm',
                        desc="mass flow rate averaged specific enthalpy")
        self.add_output('Wout', shape=1, units="lbm/s", desc="total massflow out")
        self.add_output('composition_out', val=inflow_composition_vec)


        self.declare_partials('*', '*', method='cs')

    def compute(self, inputs, outputs): 

        mix_mode = self.options['mix_mode']

        compo_in = inputs['Fl_I:tot:composition']
        n_compo = len(compo_in)


        W_in = inputs['Fl_I:stat:W']
        # composition vector is always given as vector of <something>-to-air ratios
        W_air_in = W_in/(1+np.sum(compo_in))
        W_other_in = W_air_in * compo_in

        W_out = 0 
        W_out += W_in

        W_other_out = np.zeros_like(compo_in)
        W_other_out += W_other_in

        W_air_out = W_air_in

        W_times_h = W_in*inputs['Fl_I:tot:h']

        if mix_mode == "reactant": 

            for mix_name in self.mix_names:  
                ratio = inputs[f'{mix_name}:ratio'] # scalar for reactant mode

                W_air_mix = W_air_in # for reactant mode, we reference from the incoming air
                W_other_mix = W_air_mix * ratio 
                outputs[f'{mix_name}:W'] = W_other_mix
                W_other_out[self.idx_compo] += W_other_mix
                W_out += W_other_mix
                W_times_h += W_other_mix*inputs[f'{mix_name}:h']

            outputs['composition_out'] = W_other_out/W_air_in
            outputs['Wout'] = W_out
            outputs['mass_avg_h'] = W_times_h/W_out

        else: 

            for mix_name in self.mix_names: 
                compo_mix = inputs[f'{mix_name}:composition'] # potentially a vector

                W_mix = inputs[f'{mix_name}:W']
                W_air_mix = W_mix/(1+np.sum(compo_mix))
                W_other_out += W_air_mix*compo_mix
                W_out += W_mix
                W_air_out += W_air_mix
                W_times_h += W_mix*inputs[f'{mix_name}:h']


            outputs['composition_out'] = W_other_out/W_air_out
            outputs['Wout'] = W_out
            outputs['mass_avg_h'] = W_times_h/W_out


if __name__ == "__main__":
    # Test tabular ThermoAdd in reactant mode (FAR).
    prob = om.Problem()
    prob.model.add_subsystem(
        "mix",
        ThermoAdd(spec=AIR_JETA_TAB_SPEC, mix_mode="reactant", mix_composition="FAR", mix_names="fuel"),
        promotes=["*"],
    )
    prob.model.set_input_defaults("Fl_I:stat:W", 100.0, units="lbm/s")
    prob.model.set_input_defaults("Fl_I:tot:h", 120.0, units="Btu/lbm")
    prob.model.set_input_defaults("Fl_I:tot:composition", val=[0.0])
    prob.model.set_input_defaults("fuel:h", 200.0, units="Btu/lbm")
    prob.model.set_input_defaults("fuel:ratio", 0.03)
    prob.setup()
    prob.run_model()
    print("Tabular ThermoAdd test: Wout=", prob.get_val("Wout")[0], "composition_out=", prob.get_val("composition_out"))
    print("OK")

        
