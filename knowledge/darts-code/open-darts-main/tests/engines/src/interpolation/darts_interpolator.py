import importlib
import numpy as np

# import interpolator from engines
#from darts.engines import operator_set_evaluator_iface, timer_node, value_vector, index_vector

# import interpolator from the separate shared library only containing interpolator with all combinations exposed
from engines_interpolator import operator_set_evaluator_iface, timer_node, value_vector, index_vector

"""
        Create interpolator object according to specified parameters

        Parameters
        ----------
        func : Callable
            Function to interpolate with signature: func(*vars)->float. Single operator is therefore assumed 
        axes_n_points: Iterable of integers
            The number of supporting points for each axis.
        axes_min : Iterable of floats
            The minimum value for each axis.
        axes_max : Iterable of floats
            The maximum value for each axis.
        type : string
            interpolator type:
            'multilinear' (default) - piecewise multilinear generalization of piecewise bilinear interpolation on
                                      rectangles
            'linear' - a piecewise linear generalization of piecewise linear interpolation on triangles
        type : string
            interpolator mode:
            'adaptive' (default) - only supporting points required to perform interpolation are evaluated on-the-fly
            'static' - all supporting points are evaluated during itor object construction
        platform : string
            platform used for interpolation calculations :
            'cpu' (default) - interpolation happens on CPU
            'gpu' - interpolation happens on GPU
        precision : string
            precision used in interpolation calculations:
            'd' (default) - supporting points are stored and interpolation is performed using double precision
            's' - supporting points are stored and interpolation is performed using single precision
"""


class DartsInterpolator:
    def __init__(self, func, axes_points, axes_min, axes_max, amount_of_int, algorithm: str = 'multilinear', mode: str = 'adaptive',
                 platform: str = 'cpu', precision: str = 'd'):
        # create a darts wrapper for function
        class custom_evaluator(operator_set_evaluator_iface):
            def __init__(self, func):
                super().__init__()
                self.func = func

            def evaluate(self, state, values):
                self.func(state,values)
                return 0

        # verify then inputs are valid
        assert len(axes_points) == len(axes_min)
        assert len(axes_min) == len(axes_max)
        for n_p in axes_points:
            assert n_p > 1

        self.amount_of_interpolators = amount_of_int
        # create the instance of the wrapper
        self.evaluator = custom_evaluator(func)
        self.n_dim = len(axes_points)

        # create darts timer to measure time
        self.timer = timer_node()
        self.timer.node['init'] = timer_node()

        # Use importlib mechanism to import interpolator class:
        #m = importlib.import_module('darts.engines')

        m = importlib.import_module('engines_interpolator')

        # calculate object name using 32 bit index type (i)
        # for every dimensionality a separate class exists in darts - here we evaluate its name
        # and get the class from engines, will raise AttributeError if class cannot be found

        prefix = algorithm  + '_' + mode  + '_' + platform
        suffix = '_i_' + precision + '_' + str(self.n_dim) + '_' + str(self.amount_of_interpolators)
        itor_string = prefix + '_interpolator' + suffix

        # use exposed templatized version if it exists
        if hasattr(m, itor_string):
            itor_name = getattr(m, itor_string)
            self.interpolator = itor_name(self.evaluator, index_vector(axes_points), value_vector(axes_min),
                                          value_vector(axes_max))
        else:
            print("Info: using a generic non-templatized interpolator")
            suffix = '_dynamic_i_d'
            itor_string = prefix + '_interpolator' + suffix
            itor_name = getattr(m, itor_string)
            self.interpolator = itor_name(self.evaluator, index_vector(axes_points), value_vector(axes_min),
                                          value_vector(axes_max), self.n_dim, self.amount_of_interpolators)

        # create the instance of interpolator
        # point calculation happens here, so measure it
        self.timer.node['init'].start()
        self.interpolator.init()
        self.timer.node['init'].stop()
        # create a value for a interpolate_point output
        self.amount_vector = []
        for i in range(self.amount_of_interpolators):
            self.amount_vector.append(0)
        self.values = value_vector(self.amount_vector)
        self.derivatives = value_vector(self.amount_vector * self.n_dim)

        self.interpolator.init_timer_node(self.timer)

    def interpolate_point(self, point):
        assert (len(point) == self.n_dim)
        self.interpolator.evaluate(value_vector(point), self.values)
        return self.values

    def interpolate_point_with_derivatives(self, point):
        assert (len(point) == self.n_dim)

        self.interpolator.evaluate_with_derivatives(value_vector(point), index_vector(self.amount_vector), self.values,
                                                    self.derivatives)
        return self.values, self.derivatives

    def interpolate_array(self, grid):
        if len(grid.shape) > 1:
            assert (len(grid.shape) == self.n_dim + 1)
            assert (grid.shape[0] == self.n_dim)

            result_shape = grid.shape[1:]
            states = np.moveaxis(grid, 0, -1)
            # states = states.reshape((-1, len(grid)))
        else:
            states = grid.reshape((-1, self.n_dim))
            result_shape = (len(states))

        gradient_shape = states.shape

        # flatten all axes
        states = states.flatten()
        points = value_vector(states)

        block_idx = index_vector(np.arange(int(len(states) / self.n_dim), dtype=np.int32))

        # values should fit single value per point
        values = value_vector(self.amount_vector * int((len(states) / self.n_dim)))

        # derivatives should fit self.n_dim values per point
        derivatives = value_vector(self.amount_vector * len(states))

        # interpolate and shape the result
        self.interpolator.evaluate_with_derivatives(points, block_idx, values, derivatives)

        values = np.array(values, copy=False).reshape(result_shape)
        derivatives = np.moveaxis(np.array(derivatives, copy=False).reshape(gradient_shape), -1, 0)
        return values, derivatives

    def get_init_time(self):
        return self.timer.node['init'].get_timer()

    def get_interpolation_time(self):
        return self.timer.get_timer()
