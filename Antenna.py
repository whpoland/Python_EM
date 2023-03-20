import numpy as np
import sys
from matplotlib import pyplot as plt
from abc import ABC, abstractmethod

# sys.path.append("~/williampoland/Python/EM")
sys.path.append('/Users/williampoland/Documents/GitHub/Python_EM/EM')
import Variables

class Antenna(Variables.Variables, ABC):

    theta_deg = []
    theta_rad = []
    phi_deg = []
    phi_rad = []

    def __init__(self, **kwargs):

        freq = -1
        current = 1 # shift implementation to Array class
        num = 720
        self.Antenna_inputs = {"freq": freq, "num": num, "current": current}

        for key in kwargs:
            if key in self.Antenna_inputs:
                self.Antenna_inputs[key] = kwargs[key]
            else:
                print(f"Error: {key} is not a valid parameter")

        # self.vars_obj = Variables.Variables(frequency=freq)
        Variables.Variables.__init__(self)

        num_points = self.Antenna_inputs["num"]
        Antenna.theta_deg = np.linspace(0, 359, num=num_points, endpoint=True)
        Antenna.theta_rad = np.linspace(0, 359.0/180*np.pi, num=num_points, endpoint=True)
        Antenna.phi_deg = Antenna.theta_deg
        Antenna.phi_rad = Antenna.theta_rad

        self.element_factor
        self.pattern_factor
        # TODO: allow computation of E field at a point (r,theta,phi)
        # self.E_far_field
    
        # dict will pass by reference
        self.Antenna_data = {"theta_deg": Antenna.theta_deg, "theta_rad": Antenna.theta_rad, 
                        "element_factor": self.element_factor, "pattern_factor": self.pattern_factor,
                        "field_pattern": None, "power_pattern": None}

        plt.close('all') # resest all plot windows

    @property
    @abstractmethod
    def element_factor(self):
        pass

    @property
    @abstractmethod
    def pattern_factor(self):
        pass

    # @property
    # @abstractmethod
    # def E_far_field(self):
    #     pass

    def get(self, key):
        if key in self.Antenna_data:
            return self.Antenna_data[key]
        else:
            return super().get(key)

    def set(self, key, value):
        if key in self.Antenna_data:
            self.Antenna_data[key] = value
            self.update_data()
            return True
        else:
            super().set(key, value)

    def update_data(self):
        self.Antenna_data["field_pattern"] = lambda th,phi: list(abs(np.multiply(self.get("element_factor")(th,phi), self.get("pattern_factor")(th,phi))))
        # self.Antenna_data["power_pattern"] = lambda th,phi: list(np.multiply(self.get("field_pattern")(th,phi), self.get("field_pattern")(th,phi),dtype=np.float64))
        self.Antenna_data["power_pattern"] = lambda th,phi: list(np.square(self.get("field_pattern")(th,phi)))
        # self.Antenna_data["power_pattern"] = lambda th,phi: list(np.multiply(np.square(self.get("element_factor")(th,phi)),np.square(self.get("pattern_factor")(th,phi))))

    # will output a polar plot of field or powe pattern
    def polar_plot(self, pattern_type="power", dB=True, view="both"):
        if pattern_type == "power":
            pattern_key = "power_pattern"
            plot_text = "Power Pattern"
            db_factor = 10
        elif pattern_type == "field":
            pattern_key = "field_pattern"
            plot_text = "Field Pattern"
            db_factor = 20
        else:
            print("Invalid pattern type...plotting power pattern")
            pattern_key = "power_pattern"
            plot_text = "Power Pattern"
            db_factor = 10

        # define some shorthand
        pat = self.get(pattern_key)
        theta = Antenna.theta_rad
        phi = Antenna.phi_rad

        

        # find theta,phi angles that maximize the pattern
        # NOTE: we must pass theta,phi arguments as 1 element lists for the lambda functions to work properly
        phi_max = [0] # start by assuming max value at phi=0
        if max(abs(np.array(pat(Antenna.theta_rad, phi_max)))) == 0:
            phi_max = [90] # if pattern is 0 for phi=0, adjust by 90 degrees

        # find max output -> index of input that gives max output -> use index to get theta/phi values
        # first test for theta_max
        test_max_val = max(abs(np.array(pat(theta, phi_max))))
        index_of_max = pat(theta, phi_max).index(test_max_val,0) # find theta, since vals shouldn't be 0
        theta_max = [theta[index_of_max]]
        # second test for phi_max
        test_max_val = max(abs(np.array(pat(theta_max, phi))))
        index_of_max = pat(theta_max, phi).index(test_max_val,0)
        phi_max = [phi[index_of_max]]
        if max(abs(np.array(pat(theta_max, phi_max)))) == 0:
            print("Warning: unable to determine max theta,phi angles...plot may be incorrect")
        
        # now compute pattern
        if view == "elevation":
            pattern = pat(theta, phi_max)
            x_var = theta # identify the variable we are plotting against
            plot_text += " vs Elevation"
        elif view == "azimuth":
            pattern = pat(theta_max, phi)
            x_var = phi
            plot_text += " vs Azimuth"
        elif view == "both":
            pattern_1 = pat(theta, phi_max)
            x_var_1 = theta
            pattern_2 = pat(theta_max, phi)
            x_var_2 = phi
            plot_text_1 = "Elevation"
            plot_text_2 = "Azimuth"
        else:
            print("Invalid view...plotting elevation")
            pattern = pat(theta, phi_max)
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


    def linear_plot(self, pattern_type="power"):
        return

