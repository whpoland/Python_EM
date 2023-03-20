import numpy as np
import sys

# sys.path.append("~/williampoland/Python/EM")
sys.path.append('/Users/williampoland/Documents/GitHub/Python_EM/EM')
import Constants

class Variables(Constants.Constants):


    def __init__(self, **kwargs):

        Constants.Constants.__init__(self)

        freq = -1
        epsr = 1
        # self.loss // dielectric loss tangent
        mur = 1

        self.Variables_vars = {"frequency": freq,  "epsr": epsr,  "mur": mur}

        for key in kwargs:
            if key in self.Variables_vars:
                self.Variables_vars[key] = kwargs[key]
            else:
                print(f"Error: {key} is not a valid parameter")

        omega = 2 * np.pi * self.Variables_vars["frequency"]
        eps = super().get("eps0") * self.Variables_vars["epsr"]
        mu = super().get("mu0") * self.Variables_vars["mur"]
        imp = np.sqrt(mu/eps)
        beta = omega * np.sqrt(mu * eps)
        wavelen = 2 * np.pi / beta
        vel = self.Variables_vars["frequency"] * wavelen

        self.Variables_data = {"eps": eps, "mu": mu, "imp": imp, "beta": beta,
                     "wavelength": wavelen, "velocity": vel}
        
    def update_data(self):
        omega = 2 * np.pi * self.Variables_vars["frequency"]
        eps = super().get("eps0") * self.Variables_vars["epsr"]
        mu = super().get("mu0") * self.Variables_vars["mur"]
        imp = np.sqrt(mu / eps)
        beta = omega * np.sqrt(mu * eps)
        wavelen = 2 * np.pi / beta
        vel = self.Variables_vars["frequency"] * wavelen

        self.Variables_data = {"eps": eps, "mu": mu, "imp": imp, "beta": beta,
                     "wavelength": wavelen, "velocity": vel}

    def get(self, key):
        # print("I GOT HERE")
        # print(f"vars: {self.Variables_vars}")
        # print(f"data: {self.Variables_data}")
        if key in self.Variables_vars:
            return self.Variables_vars[key]
        if key in self.Variables_data:
            return self.Variables_data[key]
        else:
            return super().get(key)

    def set(self, key, value):
        if key in self.Variables_vars:
            self.Variables_vars[key] = value
            self.update_data()
            return True
        else:
            print(f"Error: {key} is not a mutable parameter")
            return False

        
# TEST
# my_vars = Variables(frequency=200e6)
# print("---STARTING TESTS---")
# my_freq = my_vars.get("frequency")
# print(f"frequency: {my_freq}")
# my_eps = my_vars.get("eps")
# print(f"eps: {my_eps}")
# my_mu = my_vars.get("mu")
# print(f"mu: {my_mu}")
# my_beta = my_vars.get("beta")
# print(f"beta: {my_beta}")
# my_wavelength = my_vars.get("wavelength")
# print(f"wavelength: {my_wavelength}")
# free_imp = my_vars.get("imp0")
# print(f"freespace impedance: {free_imp}")
# free_eps = my_vars.get("eps0")
# print(f"freespace epsilon: {free_eps}")
# my_imp = my_vars.get("imp")
# print(f"impedance: {my_imp}")
