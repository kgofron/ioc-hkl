import math
import numpy as np
import pandas as pd
import gi
from gi.repository import GLib
gi.require_version('Hkl', '5.0')
from gi.repository import Hkl
#from epics import PV
import inspect
from enum import IntEnum

class hklCalculator_E4CV():
    def __init__(self, num_axes_solns=15):
        # initials
        self.wavelength = 0
        self.geom_name = 'E4CV'
        self.geometry = np.nan # hkl object placeholder
        self.detector = np.nan # hkl object placeholder
        self.sample = np.nan # hkl object placeholder
        self.engines = np.nan # hkl object placeholder
        self.latt = [0, 0, 0, 0, 0, 0] 
        # ^ [a1, a2, a3, alpha, beta, gamma]
        
        # sample orientation
        self.refl1_input = [0, 0, 0, 0, 0, 0, 0] 
        self.refl2_input = [0, 0, 0, 0, 0, 0, 0] 
        self.refl1 = np.nan # hkl object placeholder
        self.refl2 = np.nan # hkl object placeholder
        # ^ [h, k l, omega, chi, phi, tth]
        self.UB_matrix = np.zeros((3,3))

        # axes
        self.num_axes_solns = num_axes_solns
        self.axes_omega = 0
        self.axes_chi = 0
        self.axes_phi = 0
        self.axes_tth = 0

        # pseduoaxes 
        self.pseudoaxes_h = 0
        self.pseudoaxes_k = 0
        self.pseudoaxes_l = 0
        self.pseudoaxes_psi = 0
        self.pseudoaxes_q = 0
        self.pseudoaxes_incidence = 0
        self.pseudoaxes_azimuth1 = 0
        self.pseudoaxes_emergence = 0
        self.pseudoaxes_azimuth2 = 0
       
        # axes solutions 
        self.axes_solns_omega = []
        self.axes_solns_chi = []
        self.axes_solns_phi = []
        self.axes_solns_tth = []
        for _ in range(self.num_axes_solns):
            self.axes_solns_omega.append(0)
            self.axes_solns_chi.append(0)
            self.axes_solns_phi.append(0)
            self.axes_solns_tth.append(0)
        
        # pseudoaxes solutions
        self.pseudoaxes_solns_h = 0
        self.pseudoaxes_solns_k = 0
        self.pseudoaxes_solns_l = 0
        self.pseudoaxes_solns_psi = 0
        self.pseudoaxes_solns_q = 0
        self.pseudoaxes_solns_incidence = 0
        self.pseudoaxes_solns_azimuth1 = 0
        self.pseudoaxes_solns_emergence = 0
        self.pseudoaxes_solns_azimuth2 = 0
       
    def start(self):
        self.detector = Hkl.Detector.factory_new(Hkl.DetectorType(0))
        factory       = Hkl.factories()[self.geom_name]
        self.geometry = factory.create_new_geometry()
        # set real axes
        # then create sample? why is it this way?
        self.geometry.wavelength_set(self.wavelength, Hkl.UnitEnum.USER)
        self.engines = factory.create_new_engine_list()
        
        self.sample = Hkl.Sample.new("toto")
        lattice     = Hkl.Lattice.new(*self.latt)
        self.sample.lattice_set(lattice)
        
       
    def forward(self):
        '''
        forward hkl calculation, real -> reciprocal
        inputs
            wavelength :float:
            lattice :: basis vectors of crystal lattice in radians
            geom :str: instrument specific geometry. Options: E4CV (4-circle) ...
            value_w :list: takes in a list of 6 elements corresponding to ...

        outputs (for E4CV)
            values_hkl :list: (h,k,l)
        '''
        #TODO, make it so that I don't need to run all of start() for every forward()
        print("Forward function start")
        self.reset_pseudoaxes_solns()
        self.detector = Hkl.Detector.factory_new(Hkl.DetectorType(0))
        factory       = Hkl.factories()[self.geom_name]
        self.geometry = factory.create_new_geometry()
 
        values_w = [self.axes_omega, self.axes_chi, self.axes_phi, self.axes_tth] 
        try:
            self.geometry.axis_values_set(values_w, Hkl.UnitEnum.USER)
        except:
            print("invalid axes values")
            #TODO catch different types of errors
            return

        self.geometry.wavelength_set(self.wavelength, Hkl.UnitEnum.USER)
        
        self.sample = Hkl.Sample.new("toto")
        lattice     = Hkl.Lattice.new(*self.latt)
        self.sample.lattice_set(lattice)
         
        self.engines = factory.create_new_engine_list()
        self.engines.init(self.geometry, self.detector, self.sample)
        self.engines.get()
        hkl = self.engines.engine_get_by_name("hkl")
       
        values_hkl = hkl.pseudo_axis_values_get(Hkl.UnitEnum.USER)
        self.pseudoaxes_solns_h, self.pseudoaxes_solns_k, self.pseudoaxes_solns_l = values_hkl

    def backward(self):
        '''
        backward hkl calculation, reciprocal -> real
        inputs
            wavelength :float:
            latt :: basis vectors of crystal lattice in radians
            geometry :str: instrument specific geometry. Options: E4CV (4-circle) ...
            value_hkl :list: 

        outputs
            dependent on geometry, for 4-circle: [omega, chi, phi, 2theta]
        '''
        print("Backward function start")
        self.reset_axes_solns()
        values_hkl = [self.pseudoaxes_h, self.pseudoaxes_k, self.pseudoaxes_l]
        self.engines.init(self.geometry, self.detector, self.sample)
        self.engines.get()
        hkl = self.engines.engine_get_by_name("hkl") 
        solutions = hkl.pseudo_axis_values_set(values_hkl, Hkl.UnitEnum.USER)
        values_w_all = []
        for i, item in enumerate(solutions.items()):
            #print("id: ", i) # for kappa 6-circle
            read = item.geometry_get().axis_values_get(Hkl.UnitEnum.USER)
            #print("motor axes solution values: ", read)
            values_w_all.append(read)
        for i in range(self.num_axes_solns):
            self.axes_solns_omega[i], self.axes_solns_chi[i], \
            self.axes_solns_phi[i], self.axes_solns_tth[i] = values_w_all[i]           

    def get_axes(self):
        axes = (self.axes_omega, self.axes_chi, self.axes_phi, self.axes_tth)
        return axes

    def get_pseudoaxes(self):
        pseudoaxes = (self.pseudoaxes_h, \
                      self.pseudoaxes_k, \
                      self.pseudoaxes_l)
        return pseudoaxes

    def get_axes_solns(self, cleanprint=False):
        '''
        Which rotations are returned depends on geometry, 4-circle currently
        inputs
            cleanprint :bool: 
                if true, dataframe 
                if false, dict
        outputs
            axes :dict: 
                key is rotation name, value is list of solutions up to num_axes_solns
        '''
        axes = {}
        axes['omega'] = self.axes_solns_omega
        axes['chi']   = self.axes_solns_chi
        axes['phi']   = self.axes_solns_phi
        axes['tth']   = self.axes_solns_tth
        if(cleanprint==False):
            return axes
        elif(cleanprint==True):
            axes_df = pd.DataFrame(axes)
            return axes_df

    def get_pseudoaxes_solns(self):
        pseudoaxes_solns = (self.pseudoaxes_solns_h, \
                      self.pseudoaxes_solns_k, \
                      self.pseudoaxes_solns_l)
        return pseudoaxes_solns

    def reset_pseudoaxes_solns(self):
        self.pseudoaxes_solns_h = 0
        self.pseudoaxes_solns_k = 0
        self.pseudoaxes_solns_l = 0
        self.pseudoaxes_solns_psi = 0
        self.pseudoaxes_solns_q = 0
        self.pseudoaxes_solns_incidence = 0
        self.pseudoaxes_solns_azimuth1 = 0
        self.pseudoaxes_solns_emergence = 0
        self.pseudoaxes_solns_azimuth2 = 0

    def reset_axes_solns(self):
        self.axes_solns_omega = []
        self.axes_solns_chi = []
        self.axes_solns_phi = []
        self.axes_solns_tth = []
        for _ in range(self.num_axes_solns):
            self.axes_solns_omega.append(0)
            self.axes_solns_chi.append(0)
            self.axes_solns_phi.append(0)
            self.axes_solns_tth.append(0)
 
       
    def compute_UB_matrix(self):
        print("Computing UB matrix")
        self.add_reflection1()
        self.add_reflection2()
        self.sample.compute_UB_busing_levy(self.refl1, self.refl2)
        UB = self.sample.UB_get()
        #TODO probably a better way to do this
        for i in range(3):
            for j in range(3):
                self.UB_matrix[i,j] = UB.get(i,j)

    def add_reflection1(self):
        self.axes_omega = self.refl1_input[3]
        self.axes_chi   = self.refl1_input[4]
        self.axes_phi   = self.refl1_input[5]
        self.axes_tth   = self.refl1_input[6]
        self.forward() # replace with an update of sample with motor positions
        self.refl1 = Hkl.SampleReflection(self.geometry, self.detector, self.refl1_input[0], \
                                         self.refl1_input[1], self.refl1_input[2])

    def add_reflection2(self):
        self.axes_omega = self.refl2_input[3]
        self.axes_chi   = self.refl2_input[4]
        self.axes_phi   = self.refl2_input[5]
        self.axes_tth   = self.refl2_input[6]
        self.forward()
        self.refl2 = Hkl.SampleReflection(self.geometry, self.detector, self.refl2_input[0], \
                                         self.refl2_input[1], self.refl2_input[2])

    def test(self):
        # starting sample, instrument parameters
        self.wavelength = 1.54 #Angstrom
        self.geom_name = 'E4CV' # 4-circle
        self.latt = [1.54, 1.54, 1.54,
                math.radians(90.0),
                math.radians(90.0),
                math.radians(90.)] # cubic
        self.start() # start hkl
       
        # forward test
        self.axes_omega = 30.
        self.axes_chi   = 0.
        self.axes_phi   = 0.
        self.axes_tth   = 60.
        print("Running test - Initial conditions:")
        print("##################################")
        print("wavelength: ", self.wavelength)
        print("geometry: ", self.geom_name)
        print("lattice: ", self.latt)
        self.forward() 
        f_results = self.get_pseudoaxes_solns()
        print("input motor values: ", self.get_axes()) 
        print("forward function results:\n", f_results)  
        # backward test
        self.pseudoaxes_h = 0.
        self.pseudoaxes_k = 0.
        self.pseudoaxes_l = 1.
        values_hkl = [self.pseudoaxes_h, self.pseudoaxes_k, self.pseudoaxes_l]
        self.backward()
        b_results = self.get_axes_solns(cleanprint=True)
        print("input hkl values: ", values_hkl)
        print("backward function results:\n", b_results)
        #UB matrix test
        # Hkl.SampleReflection(self.geometry, self.detector, h, k, l)
        # When this function takes in self.geometry, it pulls the axes positions from there
        # So, need to run a forward() with reflection1 motor positions to capture reflections
        # reflection #1
        self.refl1_input[3] = -145.451 # omega
        self.refl1_input[4] = 0 # chi
        self.refl1_input[5] = 0 # phi
        self.refl1_input[6] = 69.0966 # tth
        self.refl1_input[0] = 4 # h
        self.refl1_input[1] = 0 # k
        self.refl1_input[2] = 0 # l
        # reflection #2
        self.refl2_input[3] = -145.451 # omega
        self.refl2_input[4] = 90 # chi
        self.refl2_input[5] = 0 # phi
        self.refl2_input[6] = 69.0966 # tth
        self.refl2_input[0] = 0 # h
        self.refl2_input[1] = 4 # k
        self.refl2_input[2] = 0 # l
        # Finally, compute the UB matrix #TODO calculated matrix not correct
        self.compute_UB_matrix()
        print("Testing UB matrix calculation")
        print(f'reflection #1: {self.refl1}')
        print(f'reflection #2: {self.refl2}')
        print(f'Resulting UB matrix: {self.UB_matrix}')