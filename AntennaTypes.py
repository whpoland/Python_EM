import numpy as np
import sys
from matplotlib import pyplot

# sys.path.append("~/williampoland/Python/EM")
sys.path.append('/Users/williampoland/Documents/GitHub/Python_EM/EM')
import Antenna

class IsotropicAntenna(Antenna.Antenna):

    # necessary to assign these for instantiation
    element_factor = None
    pattern_factor = None
    # TODO: add expression for (far-field) E

    def __init__(self, **kwargs):

        freq = -1
        current = 1
        num = 720
        self.IsotropicAntenna_inputs = {"freq": freq,"num": num, "current": current}

        for key in kwargs:
            if key in self.IsotropicAntenna_inputs:
                self.IsotropicAntenna_inputs[key] = kwargs[key]
            else:
                print(f"Error: {key} is not a valid parameter")

        super().__init__(freq=self.IsotropicAntenna_inputs["freq"], num=self.IsotropicAntenna_inputs["num"])

        # note th,phi are arrays; operations should be element-wise
        my_ef = lambda th,phi: [1] * max(len(th),len(phi))
        my_pf = lambda th,phi: [1] * max(len(th),len(phi))
        super().set("element_factor", my_ef)
        super().set("pattern_factor", my_pf)


class IdealDipole(Antenna.Antenna):

    # necessary to assign these for instantiation
    element_factor = None
    pattern_factor = None

    def __init__(self, **kwargs):

        freq = -1
        current = 1
        num = 720
        self.IdealDipole_inputs = {"freq": freq,"num": num, "current": current}

        for key in kwargs:
            if key in self.IdealDipole_inputs:
                self.IdealDipole_inputs[key] = kwargs[key]
            else:
                print(f"Error: {key} is not a valid parameter")

        super().__init__(freq=self.IdealDipole_inputs["freq"], num=self.IdealDipole_inputs["num"])

        # note th,phi are arrays; operations should be element-wise
        my_ef = lambda th,phi: np.sin(th)
        my_pf = lambda th,phi: [1] * max(len(th),len(phi))
        super().set("element_factor", my_ef)
        super().set("pattern_factor", my_pf)

class UniformLineSource(Antenna.Antenna):

    # necessary to assign these for instantiation
    element_factor = None
    pattern_factor = None

    def __init__(self, **kwargs):

        freq = -1
        current = 1
        wavelengths = 1 # length in terms of wavelength
        num = 720
        self.UniformLineSource_inputs = {"freq": freq,"num": num, "current": current, "wavelengths": wavelengths}

        for key in kwargs:
            if key in self.UniformLineSource_inputs:
                self.UniformLineSource_inputs[key] = kwargs[key]
            else:
                print(f"Error: {key} is not a valid parameter")

        super().__init__(freq=self.UniformLineSource_inputs["freq"])

        length = wavelengths * self.get("wavelength")
        beta = self.get("beta")
        BL2 = length * beta / 2

        # NOTE: th,phi are arrays; operations should be element-wise
        # NOTE: if th OR phi is a vector, the output should be a vector
        my_ef = lambda th,phi: np.sin(th)
        my_pf = lambda th,phi: (np.subtract(phi,phi)+1) * np.sin(BL2 * np.cos(th)) / (BL2 * np.cos(th)) # (np.subtract(phi,phi)+1) term preserves vector output
        super().set("element_factor", my_ef)
        super().set("pattern_factor", my_pf)

class PatchAntenna(Antenna.Antenna):

    # necessary to assign these for instantiation
    element_factor = None
    pattern_factor = None

    def __init__(self, **kwargs):

        freq = -1
        current = 1
        num = 720
        self.PatchAntenna_inputs = {"freq": freq,"num": num, "current": current}

        for key in kwargs:
            if key in self.PatchAntenna_inputs:
                self.PatchAntenna_inputs[key] = kwargs[key]
            else:
                print(f"Error: {key} is not a valid parameter")

        super().__init__(freq=self.PatchAntenna_inputs["freq"], num=self.PatchAntenna_inputs["num"])

        # note th,phi are arrays; operations should be element-wise
        my_ef = lambda th,phi: 0 # TODO: implement
        my_pf = lambda th,phi: 0 # TODO: implement
        super().set("element_factor", my_ef)
        super().set("pattern_factor", my_pf)

class NotchAntenna(Antenna.Antenna):

    # necessary to assign these for instantiation
    element_factor = None
    pattern_factor = None

    def __init__(self, **kwargs):

        freq = -1
        current = 1
        num = 720
        self.NotchAntenna_inputs = {"freq": freq,"num": num, "current": current}

        for key in kwargs:
            if key in self.NotchAntenna_inputs:
                self.NotchAntenna_inputs[key] = kwargs[key]
            else:
                print(f"Error: {key} is not a valid parameter")

        super().__init__(freq=self.NotchAntenna_inputs["freq"], num=self.NotchAntenna_inputs["num"])

        # note th,phi are arrays; operations should be element-wise
        my_ef = lambda th,phi: 0 # TODO: implement
        my_pf = lambda th,phi: 0 # TODO: implement
        super().set("element_factor", my_ef)
        super().set("pattern_factor", my_pf)


# test_ant = IsotropicAntenna(freq=3e6)
# test_ant.polar_plot()
# ant2 = IdealDipole(freq=3e6)
# ant2.polar_plot(view="both")
# ant3 = UniformLineSource(freq=3e6,num=720)
# ant3.polar_plot(view="both")