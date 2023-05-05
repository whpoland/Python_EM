import numpy as np
import sys
from matplotlib import pyplot as plt

# sys.path.append("~/williampoland/Python/EM")
sys.path.append('/Users/williampoland/Documents/GitHub/Python_EM/EM')
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
        vaild_types = ['isotropic', 'idealdipole', 'uniformlinesource', 'patch', 'notch']
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
            if self.Array_inputs[var] and len(self.Array_inputs[var]) != self.Array_inputs["num"]:
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



        self.update_array_data()
        

        # array factor (AF) & element factor (EF)
        # NOTE: array EF = antenna field pattern = antenna ef*pf
        self.Array_AF = self.comp_array_factor # function (theta,phi-unused)
        self.Array_EF = self.Array_obj.get("field_pattern") # lambda function (theta,phi)

        self.Array_fp = lambda th,phi: list(abs(np.multiply(self.Array_AF(th,phi),self.Array_EF(th,phi))))
        self.Array_pp = lambda th,phi: list(abs(np.square(self.Array_fp(th,phi))))


        # plt.close('all') # doesn't work?


        # check data
        print("Check Instantiated Array Data:")
        print("(magnitude, phase, position)")
        for i in range(self.Array_inputs["num"]):
            print(f"Element {i}: (" ,end=" ")
            vars = ["magnitude", "phase", "position"]
            for var in vars:
                print(f"{self.Array_data[i][var]}, ", end="")
                # print(f"{var} of {i}: {self.Array_data[i][var]}")
            print(")")


    def update_array_data(self):
        # now construct a list with all of the array data
        # each array element is a dictionary: current 'magnitude', current 'phase', 'position'
        # initialize as empty

        # self.Array_data = [{'magnitude': None, 'phase': None, 'position': None, }] * self.Array_inputs["num"] # don't do this! list multiplication doesn't work the same as list comprehension
        self.Array_data = [{'magnitude': None, 'phase': None, 'position': None, } for xx in range(self.Array_inputs["num"])]
        
        
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
        return


    def comp_array_factor(self, input_array, unused_phi_input, scalar=False, normalize=True):
        """
        Computes the Array Factor (AF) for the array object that calls it
        param: input array - input theta values; can be an array of a particular value
        param: phi_input - not currently used in the function; but is added so the function's parameters match other ef/af function calls
        param: scalar (optional) [defaults to False]- must be set to true if input_array is a single value
        param: normalize (optional) [defaults to True] - normalized output array by max value when True
        """
        # note: phi input argument is declared so we can interchange this function with the EF, FP, PP calls
        if scalar and len(input_array) != 1:
            print("Error: calling scalar version of comp_array_factor() with array input argument")
            return -1

        af = np.zeros(len(input_array))
        arr = self.Array_data  # shorthand

        beta = self.get("beta")
        theta = input_array
        for i in range(self.Array_inputs["num"]):
            curr_mag = self.Array_data[i]["magnitude"]
            curr_phase = np.pi / 180 * self.Array_data[i]["phase"]
            d = self.Array_data[i]["position"] * self.get("wavelength") # NOTE: position is in terms of wavelength so multiply by wavelength
            
            # print(self.Array_data[i]["position"])
            # print(self.get("frequency"))
            # print(self.get("wavelength"))
            # print(f"#{i}| d:{d},phase:{curr_phase},theta:{theta}")

            elem = curr_mag * np.exp(1j * (curr_phase + beta * d * np.cos(theta))) # calculate element contribution to AF (for each angle in [theta])
            # elem = np.real(elem) # TODO: how to handle complex numbers in plotting
            af = np.add(af, elem) # sum element contribution to total AF
        
        # allow this as optional when we call this function to sweep AF for different input parameters
        if (normalize):
            af = np.divide(af, max(af)) # normalize

        af = np.real(af) # take the real part
        af = np.abs(af) # 
        
        # we want to extend the size of AF array if we are only passed 1 theta input
        if (not scalar) and len(input_array) == 1:
            af = np.full(len(unused_phi_input),af[0])

        af = list(af)  # convert to a list  

        if scalar:
            af = af[0]

        # print(f"AF: {af}")

        return af


    def get(self, key, element=-1):
        if key == 'num':
            return self.Array_inputs["num"]
        elif element > 0 and element < self.Array_inputs["num"]:
            if key in self.Array_data[element]:
                return self.Array_data[element][key]
        else:
            return self.Array_obj.get(key)


    def set(self, key, value):
        if key == 'freq':
            print("Warning: changing frequency after initialization could cause unexpected behavior")
        if key in self.Array_inputs:
            self.Array_inputs[key] = value
            self.update_array_data()
            print("updated array data")
            return True
        else:
            self.Array_obj.set(key,value)


    def polar_plot(self, pattern_type="power", dB=True, view="both"):
        plot_text = "{} Element ".format(self.Array_inputs["num"])

        if pattern_type == "power":
            pattern_func = self.Array_pp
            plot_text += "Power Pattern"
            db_factor = 10
        elif pattern_type == "field":
            pattern_func = self.Array_fp
            plot_text += "Field Pattern"
            db_factor = 20
        elif pattern_type == "af":
            pattern_func = self.Array_AF
            plot_text += "Array Factor"
            db_factor = 20
        elif pattern_type == "ef":
            pattern_func = self.Array_EF
            plot_text += "Element Factor"
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
        if np.max(np.abs(np.array(pattern_func(theta, phi_max)))) == 0:
            phi_max = [90] # if pattern is 0 for phi=0, adjust by 90 degrees

        # find max output -> index of input that gives max output -> use index to get theta/phi values
        # first test for theta_max
        test_max_val = np.max(abs(np.array(pattern_func(theta, phi_max))))
        index_of_max = pattern_func(theta, phi_max).index(test_max_val,0) # find theta, since vals shouldn't be 0
        theta_max = [theta[index_of_max]]
        # second test for phi_max
        test_max_val = np.max(np.abs(np.array(pattern_func(theta_max, phi))))
        # print(f"test val: {test_max_val}")
        # print(pattern_func(theta_max, phi))
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

        if(self.Array_inputs["uniform"]):
            plot_text += "\n(alpha = {}, spacing = {}{})".format(self.Array_inputs["alpha"],self.Array_inputs["spacing"],chr(0x03BB))

        # setup a figure and axes for polar plot
        if view == "both":
            # convert to log scale
            if dB:
                np.seterr(divide = 'ignore') # otherwise log10 will through 'RuntimeWarning: divide by zero encountered'
                # np.seterr(divide = 'warn') # turn warning back on
                pattern_1 = db_factor * np.log10(pattern_1, dtype=np.float32)
                pattern_2 = db_factor * np.log10(pattern_2, dtype=np.float32)
                plot_text += " [dB]"

            fig, (ax1, ax2) = plt.subplots(1, 2, subplot_kw={'projection': 'polar'}, num="Polar Plot")
            fig.suptitle(plot_text)
            ax1.set_title(plot_text_1)
            ax2.set_title(plot_text_2)
            ax1.plot(x_var_1, pattern_1,
                    linewidth=3)
            ax2.plot(x_var_2, pattern_2,
                    linewidth=3)
            if dB:
                # pass
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
                plot_text += " [dB]"

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
        plt.draw()
        plt.show() # display the plot in a window
        return

    
    def sweep_af(self, input_angle, min_spacing=0.01,max_spacing=1,min_alpha=0,max_alpha=359):
        """
        sweep AF at a particular angle for differnt spacing/current phases (normalizes output)
        param: input angle (deg)
        param: (optional): min_spacing, max_spacing, min_alpha, max_alpha
        return: void (only displays output plots)
        """
        
        if not self.Array_inputs["uniform"]:
            print("Error: sweep_af() function does not work on nonuniform arrays")
            return -1

        # convert angle to radians
        angle_rad = input_angle * np.pi / 180

        # orig_freq = self.Array_inputs["freq"]
        # orig_num = self.Array_inputs["num"]
        # orig_points = self.Array_inputs["points"]
        # orig_type = self.Array_inputs["type"]
        orig_spacing = self.Array_inputs["spacing"]
        orig_alpha = self.Array_inputs["alpha"]


        # AF vs spacing 
        x_spacing = np.linspace(min_spacing, max_spacing, num=self.Array_inputs["points"], endpoint=True)

        af_spacing = np.zeros(len(x_spacing))
        for i in range(len(x_spacing)):
            self.set("spacing", x_spacing[i])
            af_spacing[i] = self.comp_array_factor([angle_rad],None,scalar=True, normalize=False)
        af_spacing = np.divide(af_spacing, max(af_spacing)) # normalize


        # AF vs alpha
        x_alpha = np.linspace(min_alpha, max_alpha, num=self.Array_inputs["points"], endpoint=True)

        af_alpha = np.zeros(len(x_alpha))
        for i in range(len(x_alpha)):
            self.set("alpha", x_alpha[i])
            af_alpha[i] = self.comp_array_factor([angle_rad],None,scalar=True, normalize=False)
        af_alpha = np.divide(af_alpha, max(af_alpha)) # normalize

        # convert to dB
        af_spacing = 20 * np.log10(af_spacing)
        af_alpha = 20 * np.log10(af_alpha)

        fig, (ax1, ax2) = plt.subplots(1, 2, num="AF Sweep")
        num_elems = self.Array_inputs["num"]
        fig.suptitle(f"{num_elems} Element Array Factor at {chr(920)}={input_angle}º  [dB]")
        ax1.set_title(f"AF vs Element Spacing (alpha = {orig_alpha})")
        ax2.set_title(f"AF vs alpha (spacing = {orig_spacing}{chr(0x03BB)})")
        ax1.plot(x_spacing, af_spacing,
                linewidth=5)
        ax2.plot(x_alpha, af_alpha,
                linewidth=5)


        # ax1.set_ylim(bottom=-1,top=0)
        # ax2.set_ylim(bottom=-1,top=0)

        # ax1.set_xlim(left=min_spacing,right=max_spacing)
        # ax2.set_xlim(left=min_alpha,right=max_alpha)

        ax1.grid(True)
        ax2.grid(True)

        fig.tight_layout()
        plt.draw()
        plt.show() # display the plot in a window


        # return values to their originals
        self.set("spacing",orig_spacing)
        self.set("alpha",orig_alpha)

        return


    def linear_plot(self):
        return

# arr = Array(freq=3e6, num=1, type='isotropic')
# arr = Array(freq=3e8, num=5, points=720, type='isotropic', uniform=True, spacing=1, alpha=0) 
# arr.polar_plot(pattern_type="af",view="both")
# arr.sweep_af(0)

# arr2 = Array(freq=3e8, num=2, points=720, type='idealdipole', uniform=True, spacing=0.5, alpha=0)
# arr2.polar_plot(pattern_type="power", view="both")


arr = Array(freq=5.775e9, num=4, points=720, type='idealdipole', uniform=True, spacing=1, alpha=0) 
# arr = Array(freq=5.775e9, num=4, points=720, type='patch', uniform=False, magnitude=[1,1,1,1], position=[0,0.5,1,1.5],phase=[270,270,270,270])
# arr.sweep_af(90)
arr.polar_plot(pattern_type="af",view="elevation")