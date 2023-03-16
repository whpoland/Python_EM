import numpy as np


class Constants:

    constants = {}

    def __init__(self):
    
        eps0 = 8.854e-12
        mu0 = np.pi * 4e-7
        imp0 = np.sqrt(mu0/eps0)
        vel0 = 3e8

        Constants.constants = {"eps0": eps0, "mu0": mu0, "imp0": imp0, "vel0": vel0}

    def list_all(self):
        for key in Constants.constants:
            print(f"{key}, ")
        print("\n")
        return

    def get(self, key):
        if key in Constants.constants:
            return Constants.constants[key]
        else:
            print(f"Error: constant {key} not found")
            return None


