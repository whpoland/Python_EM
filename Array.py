import numpy as np
import sys
from matplotlib import pyplot as plt

sys.path.append("~/williampoland/Python/EM")
import AntennaTypes

# class Array(Antenna.Antenna):
class Array():

    theta_deg = []
    theta_rad = []
    phi_deg = []
    phi_rad = []

    def __init__(self, **kwargs):

        # TODO: implement way to change eps & mu
        # TODO: implement spacing in 2 dimensions
    
        freq = -1           # must be specified by init function
        points = 720        # number of points to plot for array patterns
        num = 0             # defaults to single element if not specified
        type = ""
        magnitude = []      # will default to 1's if not specified
        phase = []          # will default to 0's if not specified
        position = []       # raises an Exception if not specified for non-uniform arrays
        uniform = True      # uniform spacing
        spacing = 0         # raises an Exception if not specified for uniform arrays
        alpha = -1          # defaults to 0 if not specified for uniform arrays
        
        
        self.Array_inputs = {"freq": freq, "points": points,                           # basic inputs
                                "num": num, "type": type,                             # kwargs for all arrays
                                "magnitude": magnitude,                               # (optional) kwargs for all arrays
                                "phase": phase, "position": position,                 # kwargs for nonuniform arrays
                                "uniform": uniform, "spacing": spacing, "alpha": alpha} # kwargs for uniform arrays
        # read kwargs
        for key in kwargs:
            if key in self.Array_inputs:
                self.Array_inputs[key] = kwargs[key]
                # print(f"set key {key} to {kwargs[key]}")
            else:
                print(f"Error: {key} is not a valid parameter") 

        Array.theta_deg = np.linspace(0, 359, num=self.Array_inputs["points"], endpoint=True)
        Array.theta_rad = np.linspace(0, 359.0/180*np.pi, num=self.Array_inputs["points"], endpoint=True)
        Array.phi_deg = Array.theta_deg
        Array.phi_rad = Array.theta_rad

        ## -------------------------------------- kwarg validation ------------------------------------------------
        vaild_types = ['isotropic', 'idealdipole', 'uniformlinesource' 'patch', 'notch']
        # check if necessary kwargs are input
        if self.Array_inputs["freq"] <= 0: 
            print("Error: frequency is invalid or not specified")
            raise Exception("Frequency is invalid")
        if self.Array_inputs["num"] < 1 or not isinstance(self.Array_inputs["num"], int):
            print("Warning: invalid number of elements...defaulting to 1")
            self.Array_inputs["num"] == 1
        if self.Array_inputs["type"] not in vaild_types:
            print("Info: array element type not specified...defaulting to 'isotropic")
            self.Array_inputs["type"] = 'isotropic'
        
        # check if uniform arrays have required inputs
        if self.Array_inputs["uniform"]:
            if self.Array_inputs["alpha"] == -1:
                print("Warning: alpha value not specified for uniform array...setting to 0")
                self.Array_inputs["alpha"] = 0
            if self.Array_inputs["spacing"] == 0:
                print("Error: spacing not specified for uniform array")
                raise Exception("Uniform array spacing not specified")
        
        # check if kwargs inputs agree with array size
        inputs_to_check = ["magnitude", "phase", "position"]
        for var in inputs_to_check:
            if self.Array_inputs["uniform"] and self.Array_inputs[var]: # if uniform & non-uniform inputs are not empty
                print("Warning: non-uniform list inputs specified for uniform array...inputs will be ignored")
            elif not self.Array_inputs["uniform"] and not self.Array_inputs[var]:
                if var == "position":
                    print(f"Error: non-uniform array was specified without necessary input '{var}'")
                    raise Exception("Necessary kwarg inputs not specified")
                else:
                    print(f"Warning: non-uniform array was specified without input '{var}'")
                    # raise Exception("Necessary kwarg inputs not specified")
            if self.Array_inputs[var] and len(self.Array_inputs[var]) != self.Array_inputs[var]:
                print(f"Error: size of input list '{var}' does not match specified number of array elements")
                raise Exception("Size mismatch with input list and number of elements")
        ## --------------------------------------------------------------------------------------------------------


        # variable to store the Antenna object for the type of antenna
        self.Array_obj = None
        if self.Array_inputs["type"] == 'isotropic':
                self.Array_obj = AntennaTypes.IsotropicAntenna(freq=self.Array_inputs["freq"])
        elif self.Array_inputs["type"] == 'idealdipole':
            self.Array_obj = AntennaTypes.IdealDipole(freq=self.Array_inputs["freq"])
        elif self.Array_inputs["type"] == 'uniformlinesource':
            self.Array_obj = AntennaTypes.UniformLineSource(freq=self.Array_inputs["freq"])
        elif self.Array_inputs["type"] == 'patch':
            self.Array_obj = AntennaTypes.PatchAntenna(freq=self.Array_inputs["freq"])
        elif self.Array_inputs["type"] == 'notch':
            self.Array_obj = AntennaTypes.NotchAntenna(freq=self.Array_inputs["freq"])
        else:
            print("Error: unrecognized antenna type...instantiating as isotropic")
            self.Array_obj = AntennaTypes.IsotropicAntenna(freq=self.Array_inputs["freq"])

        # update array data
        # for i in range(num):
        #     if type == 'isotropic':
        #         self.Array_data[i]["obj"] = AntennaTypes.IsotropicAntenna(freq=self.Array_inputs["freq"])
        #     elif type =='idealdipole':
        #         self.Array_data[i]["obj"] = AntennaTypes.IdealDipole(freq=self.Array_inputs["freq"])
        #     elif type == 'patch':
        #         self.Array_data[i]["obj"] = AntennaTypes.PatchAntenna(freq=self.Array_inputs["freq"])
        #     elif type == 'notch':
        #         self.Array_data[i]["obj"] = AntennaTypes.NotchAntenna(freq=self.Array_inputs["freq"])
        #     else:
        #         print("Error: unrecognized antenna type...setting to isotropic")
        #         self.Array_data[i]["obj"] = AntennaTypes.IsotropicAntenna(freq=self.Array_inputs["freq"])


        # now construct a list with all of the array data
        # each array element is a dictionary: current 'magnitude', current 'phase', 'position'
        # initialize as empty
        self.Array_data = [{}] * self.Array_inputs["num"]

        if self.Array_inputs["uniform"]:
            for i in range(self.Array_inputs["num"]):
                if self.Array_inputs["magnitude"]:
                    self.Array_data[i]["magnitude"] = self.Array_inputs["magnitude"][i]
                else:
                    self.Array_data[i]["magnitude"] = 1
                    
                self.Array_data[i]["phase"] =  i * self.Array_inputs["alpha"]
                self.Array_data[i]["position"] = i * self.Array_inputs["spacing"]
        else:
            for i in range(self.Array_inputs["num"]):
                if self.Array_inputs["magnitude"]:
                    self.Array_data[i]["magnitude"] = self.Array_inputs["magnitude"][i]
                else:
                    self.Array_data[i]["magnitude"] = 1
                if self.Array_inputs["phase"]:
                    self.Array_data[i]["phase"] =  self.Array_inputs["phase"][i]
                else:
                    self.Array_data[i]["phase"] =  0

                self.Array_data[i]["position"] = self.Array_inputs["position"][i]

        # array factor (AF) & element factor (EF)
        # NOTE: array EF = antenna field pattern = antenna ef*pf
        self.Array_AF = self.comp_array_factor # function (theta)
        self.Array_EF = self.Array_obj.get("field_pattern") # lambda function (theta,phi)

        self.Array_fp = lambda th,phi: list(abs(np.multiply(self.Array_AF(th),self.Array_EF(th,phi))))
        self.Array_pp = lambda th,phi: list(abs(np.square(self.Array_fp(th,phi))))


        # check data
        print("HELLO")
        for i in range(self.Array_inputs["num"]):
            vars = ["magnitude", "phase", "position"]
            for var in vars:
                print(f"{var} of {i}: {self.Array_data[i][var]}")


    def comp_array_factor(self, input_array):
        af = np.zeros(len(input_array))
        arr = self.Array_data  # shorthand

        beta = self.get("beta")
        theta = input_array
        for i in range(self.Array_inputs["num"]):
            curr_mag = self.Array_data[i]["magnitude"]
            curr_phase = self.Array_data[i]["phase"]
            d = self.Array_data[i]["position"] * self.get("wavelength") # NOTE: position is in terms of wavelength; do we need to multiply by wavelength???

            elem = curr_mag * np.exp(1j * (curr_phase + beta * d * np.cos(theta))) # calculate element contribution to AF
            # elem = np.real(elem) # TODO: how to handle complex numbers in plotting
            af = np.add(af, elem) # sum element contribution to total AF

        # now normalize and take the real part
        af = np.divide(af, max(af))
        af = np.real(af)
        
        return af


    def get(self, key, element=-1):
        if key == 'num':
            return self.Array_inputs["num"]
        elif element > 0 and element < self.Array_inputs["num"]:
            if key in self.Array_data[element]:
                return self.Array_data[element][key]
        else:
            return self.Array_obj.get(key)


    def polar_plot(self, pattern_type="power", dB=True, view="both"):

        if pattern_type == "power":
            pattern_func = self.Array_pp
            plot_text = "Power Pattern"
            db_factor = 10
        elif pattern_type == "field":
            pattern_func = self.Array_fp
            plot_text = "Field Pattern"
            db_factor = 20
        else:
            print("Invalid pattern type...plotting power pattern")
            pattern_func = self.Array_pp
            plot_text = "Power Pattern"
            db_factor = 10

        # define some shorthand
        theta = Array.theta_rad
        phi = Array.phi_rad

        # find theta,phi angles that maximize the pattern
        # NOTE: we must pass theta,phi arguments as 1 element lists for the lambda functions to work properly
        phi_max = [0] # start by assuming max value at phi=0
        if max(abs(np.array(pattern_func(theta, phi_max)))) == 0:
            phi_max = [90] # if pattern is 0 for phi=0, adjust by 90 degrees

        # find max output -> index of input that gives max output -> use index to get theta/phi values
        # first test for theta_max
        test_max_val = max(abs(np.array(pattern_func(theta, phi_max))))
        index_of_max = pattern_func(theta, phi_max).index(test_max_val,0) # find theta, since vals shouldn't be 0
        theta_max = [theta[index_of_max]]
        # second test for phi_max
        test_max_val = max(abs(np.array(pattern_func(theta_max, phi))))
        index_of_max = pattern_func(theta_max, phi).index(test_max_val,0)
        phi_max = [phi[index_of_max]]
        if max(abs(np.array(pattern_func(theta_max, phi_max)))) == 0:
            print("Warning: unable to determine max theta,phi angles...plot may be incorrect")
        
        # now compute pattern
        if view == "elevation":
            pattern = pattern_func(theta, phi_max)
            x_var = theta # identify the variable we are plotting against
            plot_text += " vs Elevation"
        elif view == "azimuth":
            pattern = pattern_func(theta_max, phi)
            x_var = phi
            plot_text += " vs Azimuth"
        elif view == "both":
            pattern_1 = pattern_func(theta, phi_max)
            x_var_1 = theta
            pattern_2 = pattern_func(theta_max, phi)
            x_var_2 = phi
            plot_text_1 = "Elevation"
            plot_text_2 = "Azimuth"
        else:
            print("Invalid view...plotting elevation")
            pattern = pattern_func(theta, phi_max)
            x_var = theta
            plot_text += " vs Elevation"
        print(f"identified theta_max={theta_max}, phi_max={phi_max}")

        # setup a figure and axes for polar plot
        if view == "both":
            # convert to log scale
            if dB:
                np.seterr(divide = 'ignore') # otherwise log10 will through 'RuntimeWarning: divide by zero encountered'
                # np.seterr(divide = 'warn') # turn warning back on
                pattern_1 = db_factor * np.log10(pattern_1, dtype=np.float32)
                pattern_2 = db_factor * np.log10(pattern_2, dtype=np.float32)
                plot_text += " (dB)"

            fig, (ax1, ax2) = plt.subplots(1, 2, subplot_kw={'projection': 'polar'})
            fig.suptitle(plot_text)
            ax1.set_title(plot_text_1)
            ax2.set_title(plot_text_2)
            ax1.plot(x_var_1, pattern_1,
                    linewidth=3)
            ax2.plot(x_var_2, pattern_2,
                    linewidth=3)
            if dB:
                pass
                ax1.set_rmax(0)
                ax2.set_rmax(0)
                ax1.set_rmin(-20)
                ax2.set_rmin(-20)
                ax1.set_rticks([-20,-15,-10,-5,0])  # radial ticks
                ax2.set_rticks([-20,-15,-10,-5,0])  # radial ticks
            else:
                ax1.set_rmax(1)
                ax2.set_rmax(2)
                ax1.set_rticks([0.25,0.5,0.75,1])  # radial ticks
                ax2.set_rticks([0.25,0.5,0.75,1])  # radial ticks
            # modify plot visuals
            ax1.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
            ax2.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
            ax1.grid(True)
            ax2.grid(True)
        # single view plot
        else:
            # convert to log scale
            if dB:
                np.seterr(divide = 'ignore') # otherwise log10 will through 'RuntimeWarning: divide by zero encountered'
                # np.seterr(divide = 'warn') # turn warning back on
                pattern = db_factor * np.log10(pattern, dtype=np.float32)
                plot_text += " (dB)"

            fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
            ax.set_title(plot_text, va='bottom')
            # plot the antenna pattern
            ax.plot(x_var, pattern,
                    linewidth=1,markersize=1)
            if dB:
                # pass
                ax.set_rmax(0)
                ax.set_rmin(-20)
                # ax.set_rticks([-1,-0.5,0])  # radial ticks
            else:
                ax.set_rmax(1)
                ax.set_rticks([0.25,0.5,0.75,1])  # radial ticks
            # modify plot visuals
            ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
            ax.grid(True)

        fig.tight_layout()
        plt.show() # display the plot in a window
        return


    def linear_plot(self):
        return

# arr = Array(freq=3e6, num=1, type='isotropic')
arr = Array(freq=3e8, num=2, points=720, type='isotropic', uniform=True, spacing=0.5, alpha=0)
arr.polar_plot(pattern_type="power",view="both")