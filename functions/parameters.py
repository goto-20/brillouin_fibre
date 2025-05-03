#parameters.py
from dolfin import *
import sys

version = sys.version_info

if version.major == 3:

    class CellFunctionCases(UserExpression):
        def __init__(self, user_shape, user_mapping, user_marker_map, **kwargs):
            self.user_shape = user_shape
            self.user_mapping = user_mapping
            self.user_marker_map = user_marker_map
            super().__init__(**kwargs)
        def eval_cell(self, values, x, ufc_cell):
            mapping = self.user_mapping
            k = self.user_marker_map[ufc_cell.index]
            vals = mapping[k].values()
            for i in range(len(vals)):
                values[i] = vals[i]
        def value_shape(self):
            return self.user_shape

    def Create_CellFunctionCases(markers, mapping, **kwargs):

        shape = next(iter(mapping.values())).ufl_shape
        return CellFunctionCases(
                user_shape=shape,
                user_mapping=mapping,
                user_marker_map=markers,
                degree=1, **kwargs)

    def get_param_func(markers, constants):
        # wrapper around CellFunctionCases
        mapping = {k: Constant(v) for k,v in constants.items()}
        return Create_CellFunctionCases(markers, mapping)
else:
    def Create_CellFunctionCases(markers, mapping, **kwargs):
        shape = next(iter(mapping.values())).ufl_shape
        class CellFunctionCases(Expression):
            '''Given `cell_function` and `mapping`, this object is a dolfin Expression
            which takes the value mapping[cell_function[current_cell]] at every point. It
            effectively plays the same role as \\cases does in mathematics.'''
            def __init__(self, **kwargs):
                self.cell_function = markers
                self.mapping = mapping
                # we're on purpose not calling superclass __init__.
                # dolfin's expression.py says so, because it's doing some black magic

            def eval_cell(self, values, x, cell):
                mapping = self.mapping
                k = self.cell_function[cell.index]
                if k not in mapping:
                    k = None
                mapping[k].eval_cell(values, x, cell)

            def value_shape(self):
                return shape

        return CellFunctionCases(degree=1,**kwargs)

    def get_param_func(markers, constants):
        # wrapper around CellFunctionCases
        mapping = {k: Constant(v) for k,v in constants.items()}
        return Create_CellFunctionCases(markers, mapping)

def get_params(markers, material):

    if material == 'silica':
        Dielectric    = get_param_func(markers, {1: (1.4492**2), 0: (1.444**2)})
        Density       = get_param_func(markers, {1: (2.254), 0: (2.203)})
        Stress_Tensor = get_param_func(markers, {1: (76, 16.15, 29.9),
                                                 0: (78, 16, 31)})
        Photo_Tensor  = get_param_func(markers, {1: (0.20, 0.27, -0.073),
                                                 0: (0.20, 0.27, -0.073)})
        Qf_product    = get_param_func(markers, {1: (2*pi*6e3), 0: (2*pi*1e0)})

    if material == 'chalco':
        Dielectric    = get_param_func(markers, {1: (2.674**2), 0: (1.481**2)})
        Density       = get_param_func(markers, {1: (4.640), 0: (1.187)})
        Stress_Tensor = get_param_func(markers, {1: (23.5, 9.5, 6.99),
                                                 0: (6.63, 4.34, 1.145)})
        Photo_Tensor  = get_param_func(markers, {1: (0.314, 0.266, 0.024),
                                                 0: (0.30, 0.297,  0.0015)})
        Qf_product    = get_param_func(markers, {1: (2*pi*2e3), 0: (2*pi*7e1)})

    return Dielectric, Density, Stress_Tensor, Photo_Tensor, Qf_product
