import math
import numpy as np
import datetime
import gi
from gi.repository import GLib
gi.require_version('Hkl', '5.0')
from gi.repository import Hkl
import intensities
from util import energy2wavelength_neutron

class hklCalculator():
    def __init__(self, num_axes_solns=30, num_reflections = 10, geom=1, geom_name = 'E4CV'):
        # initials
        self.wavelength = 0.
        self.geom_name = geom_name
        self.geom = geom
        self.geometry = np.nan # hkl object placeholder
        self.detector = np.nan # hkl object placeholder
        self.factory = np.nan # hkl object placeholder
        self.sample = np.nan # hkl object placeholder
        self.engines = np.nan # hkl object placeholder
        self.engine_hkl = np.nan # hkl object placeholder
        self.engine_psi = np.nan # hkl object placeholder
        self.engine_q = np.nan # hkl object placeholder
        self.engine_incidence = np.nan # hkl object placeholder
        self.engine_emergence = np.nan # hkl object placeholder
        self.engine_eulerians = np.nan # hkl object placeholder
        self.engine_psi = np.nan # hkl object placeholder
        self.engine_q2 = np.nan # hkl object placeholder
        self.engine_qper_qpar = np.nan # hkl object placeholder
        self.engine_tth2 = np.nan # hkl object placeholder
        self.mode_4c = 0 # bissector, constant omega...
        self.mode_6c = 0 # bissector_vertical, constant_omega_vertical
        self.latt = [0., 0., 0., 0., 0., 0.] 
        # ^ [a1, a2, a3, alpha, beta, gamma], angstroms and radians
        self.lattice = np.nan 
        self.lattice_vol = 0.
        
        ct = datetime.datetime.now().isoformat()
        self.errors = [ord(c) for c in str(ct)]
        self.intensities = ''
        self.cif_path=''        

        self.energy = 0.
        self.wavelength_result = 0.
        #self.particle_type = 1 # 0 photon, 1 neutron, ... #TODO
        self.neutron = 0.
        self.velocity = 0. # sqrt(2E/m) 
        self.e = 1.6021766300e-19 # [C]
        self.h = 6.6260701500e-34 # [m^2*kg/s]
        self.c = 299792458 # [m/s^2]
        self.m_neutron = 1.6749274710e-27 #[kg]
        self.m_proton = 1.6726219200e-27 #[kg]
        self.m_electron = 9.1093837e-27 #[kg]

        # sample orientation
        # initial 2 reflections
        self.num_reflections = num_reflections
        self.refl1_input_e4c = [0., 0., 0., 0., 0., 0., 0.]
        self.refl2_input_e4c = [0., 0., 0., 0., 0., 0., 0.]
        self.refl1_input_k4c = [0., 0., 0., 0., 0., 0., 0.]
        self.refl2_input_k4c = [0., 0., 0., 0., 0., 0., 0.]
        self.refl1_input_e6c = [0., 0., 0., 0., 0., 0., 0., 0., 0.]
        self.refl2_input_e6c = [0., 0., 0., 0., 0., 0., 0., 0., 0.]
        self.refl1_input_k6c = [0., 0., 0., 0., 0., 0., 0., 0., 0.]
        self.refl2_input_k6c = [0., 0., 0., 0., 0., 0., 0., 0., 0.]
        self.refl1 = np.nan
        self.refl2 = np.nan
        

        # refine with reflections
        self.refl_refine_input_e4c = [0., 0., 0., 0., 0., 0., 0.]
        self.refl_refine_input_k4c = [0., 0., 0., 0., 0., 0., 0.]
        self.refl_refine_input_e6c = [0., 0., 0., 0., 0., 0., 0., 0., 0.]
        self.refl_refine_input_k6c = [0., 0., 0., 0., 0., 0., 0., 0., 0.]
        self.refl_refine_input_list_e4c = []
        self.refl_refine_input_list_k4c = []
        self.refl_refine_input_list_e6c = []
        self.refl_refine_input_list_k6c = []

        for i in range(self.num_reflections):
            self.refl_refine_input_list_e4c.append([0., 0., 0., 0., 0., 0., 0.])
            self.refl_refine_input_list_k4c.append([0., 0., 0., 0., 0., 0., 0.])
            self.refl_refine_input_list_e6c.append([0., 0., 0., 0., 0., 0., 0., 0., 0.])
            self.refl_refine_input_list_k6c.append([0., 0., 0., 0., 0., 0., 0., 0., 0.])

        self.refl_list_e4c = [] # store reflection objects, for deleting
        self.refl_list_k4c = [] # store reflection objects, for deleting
        self.refl_list_e6c = [] # store reflection objects, for deleting
        self.refl_list_k6c = [] # store reflection objects, for deleting
        self.curr_refl_index_e4c = 0
        self.curr_refl_index_k4c = 0
        self.curr_refl_index_e6c = 0
        self.curr_refl_index_k6c = 0
 
        self.selected_refl_index = 0 # reflection to delete, index of refl_list
        # UB
        self.UB_matrix = np.zeros((3,3), dtype=float)
        self.UB_matrix_input = np.zeros((3,3), dtype=float)

        #self.sample_rot_matrix = np.zeros((8,8), dtype=float)
        self.u_matrix = np.zeros((3,3), dtype=float)
        
        # U vector
        self.ux = 0.
        self.uy = 0.
        self.uz = 0.

        ### axes
        self.num_axes_solns = num_axes_solns

        # Eulerian 4-circle (omega, chi, phi, tth)
        self.axes_e4c = [0.,0.,0.,0.]
        self.axes_e4c_min = [-180.,-180.,-180.,-180.]
        self.axes_e4c_max = [180.,180.,180.,180.]

        # Kappa 4-circle (komega, kappa, kphi, tth)
        self.axes_k4c = [0.,0.,0.,0.]
        self.axes_k4c_min = [-180.,-180.,-180.,-180.]
        self.axes_k4c_max = [180.,180.,180.,180.]
        
        # Eulerian 6-circle (mu, omega, chi, phi, gamma, delta)
        self.axes_e6c = [0.,0.,0.,0.,0.,0.]
        self.axes_e6c_min = [-180.,-180.,-180.,-180.,-180.,-180.]
        self.axes_e6c_max = [180.,180.,180.,180.,180.,180.]

        # Kappa 6-circle (mu, komega, kappa, kphi, gamma, delta)
        self.axes_k6c = [0.,0.,0.,0.,0.,0.]
        self.axes_k6c_min = [-180.,-180.,-180.,-180.,-180.,-180.] 
        self.axes_k6c_max = [180.,180.,180.,180.,180.,180.] 

        ### axes for UB calculation - only used internally - avoids setting on calculation
        self.axes_UB_e4c = [0., 0., 0., 0.]
        self.axes_UB_k4c = [0., 0., 0., 0.]
        self.axes_UB_e6c = [0., 0., 0., 0., 0., 0.]
        self.axes_UB_k6c = [0., 0., 0., 0., 0., 0.]

        ### pseduoaxes 
        self.pseudoaxes_h = 0.
        self.pseudoaxes_k = 0.
        self.pseudoaxes_l = 0.
        self.pseudoaxes_dspacing = 0. #TODO get from cif2hkl waveform PV
        self.pseudoaxes_fc2 = 0. #TODO set from cif2hkl waveform PV
        self.pseudoaxes_psi = 0.
        self.pseudoaxes_q = 0.
        self.pseudoaxes_alpha = 0.
        self.pseudoaxes_qper = 0.
        self.pseudoaxes_qpar = 0.
        self.pseudoaxes_alpha = 0.
        self.pseudoaxes_alpha2 = 0.
        self.pseudoaxes_incidence = 0.
        self.pseudoaxes_azimuth1 = 0.
        self.pseudoaxes_emergence = 0.
        self.pseudoaxes_azimuth2 = 0.
        self.pseudoaxes_omega = 0.
        self.pseudoaxes_chi = 0.
        self.pseudoaxes_phi = 0.
        self.pseudoaxes_tth = 0.

        ### axes solutions 
        # Eulerian 4-circle
        self.axes_solns_omega_e4c = []
        self.axes_solns_chi_e4c = []
        self.axes_solns_phi_e4c = []
        self.axes_solns_tth_e4c = []

        # Kappa 4-circle
        self.axes_solns_komega_k4c = []
        self.axes_solns_kappa_k4c = []
        self.axes_solns_kphi_k4c = []
        self.axes_solns_tth_k4c = []

        # Eulerian 6-circle
        self.axes_solns_mu_e6c = []
        self.axes_solns_omega_e6c = []
        self.axes_solns_chi_e6c = []
        self.axes_solns_phi_e6c = []
        self.axes_solns_gamma_e6c = []
        self.axes_solns_delta_e6c = []

        # Kappa 6-circle
        self.axes_solns_mu_k6c = []
        self.axes_solns_komega_k6c = []
        self.axes_solns_kappa_k6c = []
        self.axes_solns_kphi_k6c = []
        self.axes_solns_gamma_k6c = []
        self.axes_solns_delta_k6c = []

        for _ in range(self.num_axes_solns):
            # Eulerian 4-circle
            self.axes_solns_omega_e4c.append(0)
            self.axes_solns_chi_e4c.append(0)
            self.axes_solns_phi_e4c.append(0)
            self.axes_solns_tth_e4c.append(0)
            # Kappa 4-circle
            self.axes_solns_komega_k4c.append(0)
            self.axes_solns_kappa_k4c.append(0)
            self.axes_solns_kphi_k4c.append(0)
            self.axes_solns_tth_k4c.append(0)
            # Eulerian 6-circle
            self.axes_solns_mu_e6c.append(0)
            self.axes_solns_omega_e6c.append(0)
            self.axes_solns_chi_e6c.append(0)
            self.axes_solns_phi_e6c.append(0)
            self.axes_solns_gamma_e6c.append(0)
            self.axes_solns_delta_e6c.append(0)
            # Kappa 6-circle
            self.axes_solns_mu_k6c.append(0)
            self.axes_solns_komega_k6c.append(0)
            self.axes_solns_kappa_k6c.append(0)
            self.axes_solns_kphi_k6c.append(0)
            self.axes_solns_gamma_k6c.append(0)
            self.axes_solns_delta_k6c.append(0)
        
        # pseudoaxes solutions
        self.pseudoaxes_solns_h = 0.
        self.pseudoaxes_solns_k = 0.
        self.pseudoaxes_solns_l = 0.
        self.pseudoaxes_solns_dspacing = 0. #TODO set from cif2hkl waveform PV
        self.pseudoaxes_solns_fc2 = 0. #TODO set from cif2hkl waveform PV
        self.pseudoaxes_solns_psi = 0.
        self.pseudoaxes_solns_q = 0.
        self.pseudoaxes_solns_incidence = 0.
        self.pseudoaxes_solns_azimuth1 = 0.
        self.pseudoaxes_solns_emergence = 0.
        self.pseudoaxes_solns_azimuth2 = 0.
        self.pseudoaxes_solns_omega = 0.
        self.pseudoaxes_solns_chi = 0.
        self.pseudoaxes_solns_phi = 0.
        self.pseudoaxes_solns_tth = 0.
        self.pseudoaxes_solns_alpha = 0.
        self.pseudoaxes_solns_alpha2 = 0.
        self.pseudoaxes_solns_qper = 0.
        self.pseudoaxes_solns_qpar = 0.
 
    def start(self):
        self.reset_pseudoaxes_solns()
        self.detector = Hkl.Detector.factory_new(Hkl.DetectorType(0))
        self.factory  = Hkl.factories()[self.geom_name]
        self.geometry = self.factory.create_new_geometry()
        self.geometry.wavelength_set(self.wavelength, Hkl.UnitEnum.USER)
        
        self.sample = Hkl.Sample.new("toto")
        a,b,c,alpha,beta,gamma = [i for i in self.latt]
        alpha=math.radians(alpha)
        beta=math.radians(beta)
        gamma=math.radians(gamma)
        self.lattice = Hkl.Lattice.new(a,b,c,alpha,beta,gamma)
        self.sample.lattice_set(self.lattice)             

        self.engines = self.factory.create_new_engine_list()
        self.engines.init(self.geometry, self.detector, self.sample)
        self.engines.get()

        if (self.geom == 0) or (self.geom == 1):
            self.engine_hkl = self.engines.engine_get_by_name("hkl")
            self.engine_psi = self.engines.engine_get_by_name("psi")
            self.engine_q = self.engines.engine_get_by_name("q")
            self.engine_incidence = self.engines.engine_get_by_name("incidence")
            self.engine_emergence = self.engines.engine_get_by_name("emergence")
        elif (self.geom == 2): 
            self.engine_hkl = self.engines.engine_get_by_name("hkl")
            self.engine_eulerians = self.engines.engine_get_by_name("eulerians")
            self.engine_psi = self.engines.engine_get_by_name("psi")
            self.engine_q = self.engines.engine_get_by_name("q")
            self.engine_incidence = self.engines.engine_get_by_name("incidence")
            self.engine_emergence = self.engines.engine_get_by_name("emergence")
        elif (self.geom == 3): 
            self.engine_hkl = self.engines.engine_get_by_name("hkl")
            self.engine_psi = self.engines.engine_get_by_name("psi")
            self.engine_q2 = self.engines.engine_get_by_name("q2")
            self.engine_qper_qpar = self.engines.engine_get_by_name("qper_qpar")
            self.engine_tth2 = self.engines.engine_get_by_name("tth2")
            self.engine_incidence = self.engines.engine_get_by_name("incidence")
            self.engine_emergence = self.engines.engine_get_by_name("emergence")
        elif (self.geom == 4): 
            self.engine_hkl = self.engines.engine_get_by_name("hkl")
            self.engine_eulerians = self.engines.engine_get_by_name("eulerians")
            self.engine_psi = self.engines.engine_get_by_name("psi")
            self.engine_q2 = self.engines.engine_get_by_name("q2")
            self.engine_qper_qpar = self.engines.engine_get_by_name("qper_qpar")
            self.engine_incidence = self.engines.engine_get_by_name("incidence")
            self.engine_tth2 = self.engines.engine_get_by_name("tth2")
            self.engine_emergence = self.engines.engine_get_by_name("emergence")
        self.get_UB_matrix()    
        self.get_latt_vol()
        status = self.get_info()
        ct = datetime.datetime.now().isoformat()
        status_string = f'initialized geometry {status}'
        self.errors = [ord(c) for c in str(ct + '\n' + status_string)]

    def energy_to_wavelength_neutron(self):
        self.wavelength_result = energy2wavelength_neutron(self.energy)

    def switch_geom(self):
        if self.geom == 0:
            print("switching to E4CH")
            self.geom_name = "E4CH"
            self.start()

        if self.geom == 1:
            print("switching to E4CV")
            self.geom_name = "E4CV"
            self.start()

        if self.geom == 2:
            print("switching to K4CV")
            self.geom_name = "K4CV"
            self.start()

        if self.geom == 3:
            print("switching to E6C")
            self.geom_name = "E6C"
            self.start()

        if self.geom == 4:
            print("switching to K6C")
            self.geom_name = "K6C"
            self.start()

    def forward(self):
        print("Forward function start")
        self.reset_pseudoaxes_solns()
        if (self.geom_name == 'E4CH') or (self.geom_name == 'E4CV'):
            values_w = [float(self.axes_e4c[0]), \
                        float(self.axes_e4c[1]), \
                        float(self.axes_e4c[2]), \
                        float(self.axes_e4c[3])] 
        elif self.geom_name == "K4CV":
            values_w = [float(self.axes_k4c[0]), \
                        float(self.axes_k4c[1]), \
                        float(self.axes_k4c[2]), \
                        float(self.axes_k4c[3])] 
        elif self.geom_name == "E6C":
            values_w = [float(self.axes_e6c[0]), \
                        float(self.axes_e6c[1]), \
                        float(self.axes_e6c[2]), \
                        float(self.axes_e6c[3]), \
                        float(self.axes_e6c[4]), \
                        float(self.axes_e6c[5])] 
        elif self.geom_name == "K6C":
            values_w = [float(self.axes_k6c[0]), \
                        float(self.axes_k6c[1]), \
                        float(self.axes_k6c[2]), \
                        float(self.axes_k6c[3]), \
                        float(self.axes_k6c[4]), \
                        float(self.axes_k6c[5])] 
        try:
            self.set_axis_limits()
            self.geometry.axis_values_set(values_w, Hkl.UnitEnum.USER)
            ct = datetime.datetime.now().isoformat()
            self.errors=[ord(c) for c in str(ct)] # success
        except Exception as e:
            ct = datetime.datetime.now().isoformat()
            self.errors=[ord(c) for c in str(ct + '\n' + str(e))] # failure
            print(f'forward() error: {e}') 
            return

        self.engines.get()
        ### common to all geoms
        # hkl
        values_hkl = self.engine_hkl.pseudo_axis_values_get(Hkl.UnitEnum.USER)
        self.pseudoaxes_solns_h, self.pseudoaxes_solns_k, self.pseudoaxes_solns_l = values_hkl
        # psi
        values_psi = self.engine_psi.pseudo_axis_values_get(Hkl.UnitEnum.USER)
        self.pseudoaxes_solns_psi = values_psi[0] # [0] required otherwise assigns as list
        # incidence
        values_incidence = self.engine_incidence.pseudo_axis_values_get(Hkl.UnitEnum.USER)
        self.pseudoaxes_solns_incidence, self.pseudoaxes_solns_azimuth1 = values_incidence
        # emergence
        values_emergence = self.engine_emergence.pseudo_axis_values_get(Hkl.UnitEnum.USER)
        self.pseudoaxes_solns_emergence, self.pseudoaxes_solns_azimuth2 = values_emergence
        if (self.geom == 0) or (self.geom == 1):
            # q
            values_q = self.engine_q.pseudo_axis_values_get(Hkl.UnitEnum.USER)
            self.pseudoaxes_solns_q = values_q[0]
        elif self.geom == 2:
            # q
            values_q = self.engine_q.pseudo_axis_values_get(Hkl.UnitEnum.USER)
            self.pseudoaxes_solns_q = values_q[0]
            # eulerians
            values_eulerians = self.engine_eulerians.pseudo_axis_values_get(Hkl.UnitEnum.USER)
            self.pseudoaxes_solns_omega, self.pseudoaxes_solns_chi, self.pseudoaxes_solns_phi = \
                values_eulerians
        elif self.geom == 3:
            # q2
            values_q2 = self.engine_q2.pseudo_axis_values_get(Hkl.UnitEnum.USER)
            self.pseudoaxes_solns_q, self.pseudoaxes_solns_alpha = values_q2
            #qper, qpar
            values_qper_qpar = self.engine_qper_qpar.pseudo_axis_values_get(Hkl.UnitEnum.USER)
            self.pseudoaxes_solns_qper, self.pseudoaxes_solns_qpar = values_qper_qpar
            # tth2
            values_tth2 = self.engine_tth2.pseudo_axis_values_get(Hkl.UnitEnum.USER)
            self.pseudoaxes_solns_tth2, self.pseudoaxes_solns_alpha2 = values_tth2
        elif self.geom == 4:
            # q2
            values_q2 = self.engine_q2.pseudo_axis_values_get(Hkl.UnitEnum.USER)
            self.pseudoaxes_solns_q, self.pseudoaxes_solns_alpha = values_q2
            #qper, qpar
            values_qper_qpar = self.engine_qper_qpar.pseudo_axis_values_get(Hkl.UnitEnum.USER)
            self.pseudoaxes_solns_qper, self.pseudoaxes_solns_qpar = values_qper_qpar
            # tth2
            values_tth2 = self.engine_tth2.pseudo_axis_values_get(Hkl.UnitEnum.USER)
            self.pseudoaxes_solns_tth2, self.pseudoaxes_solns_alpha2 = values_tth2
            #eulerians
            values_eulerians = self.engine_eulerians.pseudo_axis_values_get(Hkl.UnitEnum.USER)
            self.pseudoaxes_solns_omega, self.pseudoaxes_solns_chi, self.pseudoaxes_solns_phi = \
                values_eulerians
        self.get_UB_matrix()

    def forward_UB(self):
        #TODO rename? adds a reflection, kind of confusing.
        #TODO runs twice when  pressing "calculate UB
        print("Forward UB function start")
        if self.geom == 0 or 1:
            values_w = [float(self.axes_UB_e4c[0]), \
                        float(self.axes_UB_e4c[1]), \
                        float(self.axes_UB_e4c[2]), \
                        float(self.axes_UB_e4c[3])] 
        elif self.geom == 2:
            values_w = [float(self.axes_UB_k4c[0]), \
                        float(self.axes_UB_k4c[1]), \
                        float(self.axes_UB_k4c[2]), \
                        float(self.axes_UB_k4c[3])] 
        elif self.geom == 3:
            values_w = [float(self.axes_UB_e6c[0]), \
                        float(self.axes_UB_e6c[1]), \
                        float(self.axes_UB_e6c[2]), \
                        float(self.axes_UB_e6c[3]), \
                        float(self.axes_UB_e6c[4]), \
                        float(self.axes_UB_e6c[5])] 
        elif self.geom == 4:
            values_w = [float(self.axes_UB_k6c[0]), \
                        float(self.axes_UB_k6c[1]), \
                        float(self.axes_UB_k6c[2]), \
                        float(self.axes_UB_k6c[3]), \
                        float(self.axes_UB_k6c[4]), \
                        float(self.axes_UB_k6c[5])] 
        try:
            #TODO check if this is how UB calculation is done in hkl package
            # currently need to run this function when adding reflections for some reason
            self.geometry.axis_values_set(values_w, Hkl.UnitEnum.USER)
        except Exception as e:
            ct = datetime.datetime.now().isoformat()
            self.errors = [ord(c) for c in str(ct + '\n' + str(e))]
            print(f'forward_UB() error: {e}')
            return

    def backward(self):
        print("Backward function start")
        self.reset_axes_solns()
        # set mode
        if (self.geom == 0) or (self.geom==1) or (self.geom==2):
            if self.mode_4c == 0:
                self.engine_hkl.current_mode_set('bissector')
            elif self.mode_4c == 1:
                self.engine_hkl.current_mode_set('constant_omega')
            elif self.mode_4c == 2:
                self.engine_hkl.current_mode_set('constant_chi')
            elif self.mode_4c == 3:
                self.engine_hkl.current_mode_set('constant_phi')
        elif self.geom == (3 or 4):
            if self.mode_6c == 0: 
                self.engine_hkl.current_mode_set('bissector_vertical')
            if self.mode_6c == 1: 
                self.engine_hkl.current_mode_set('constant_omega_vertical')
            if self.mode_6c == 2: 
                self.engine_hkl.current_mode_set('constant_chi_vertical')
            if self.mode_6c == 3: 
                self.engine_hkl.current_mode_set('constant_phi_vertical')

        values_hkl = [float(self.pseudoaxes_h), \
                      float(self.pseudoaxes_k), \
                      float(self.pseudoaxes_l)]
        try:
            self.set_axis_limits()
            solutions = self.engine_hkl.pseudo_axis_values_set(values_hkl, Hkl.UnitEnum.USER)
            values_w_all = []
            len_solns = len(solutions.items())
            for i, item in enumerate(solutions.items()):
                read = item.geometry_get().axis_values_get(Hkl.UnitEnum.USER)
                values_w_all.append(read)
            if len_solns > self.num_axes_solns: # truncate if above max available soln slots
                len_solns = self.num_axes_solns -1
            if (self.geom == 0) or (self.geom==1):
                for i in range(len_solns): 
                    self.axes_solns_omega_e4c[i], \
                    self.axes_solns_chi_e4c[i], \
                    self.axes_solns_phi_e4c[i], \
                    self.axes_solns_tth_e4c[i] = values_w_all[i]         
            elif self.geom == 2:
                for i in range(len_solns): 
                    self.axes_solns_komega_k4c[i], \
                    self.axes_solns_kappa_k4c[i], \
                    self.axes_solns_kphi_k4c[i], \
                    self.axes_solns_tth_k4c[i] = values_w_all[i]         
            elif self.geom == 3:
                for i in range(len_solns): 
                    self.axes_solns_mu_e6c[i], \
                    self.axes_solns_omega_e6c[i], \
                    self.axes_solns_chi_e6c[i], \
                    self.axes_solns_phi_e6c[i], \
                    self.axes_solns_gamma_e6c[i], \
                    self.axes_solns_delta_e6c[i] = values_w_all[i]       
            elif self.geom==4:  
                 for i in range(len_solns): 
                    self.axes_solns_mu_k6c[i], \
                    self.axes_solns_komega_k6c[i], \
                    self.axes_solns_kappa_k6c[i], \
                    self.axes_solns_kphi_k6c[i], \
                    self.axes_solns_gamma_k6c[i], \
                    self.axes_solns_delta_k6c[i] = values_w_all[i]
            ct = datetime.datetime.now().isoformat()
            self.errors = [ord(c) for c in str(ct)] # success
        except Exception as e:
            ct = datetime.datetime.now().isoformat()
            self.errors = [ord(c) for c in str(ct + '\n' + str(e))] # failure
            print(f'backward() error: {e}')

    def set_axis_limits(self):
        axes = self.geometry.axis_names_get()
        if (self.geom==0) or (self.geom==1):
            for i, axis in enumerate(axes): 
                tmp = self.geometry.axis_get(axis)
                tmp.min_max_set(self.axes_e4c_min[i], \
                                self.axes_e4c_max[i], \
                                Hkl.UnitEnum.USER)
                self.geometry.axis_set(axis, tmp)
        elif self.geom==2:
            for i, axis in enumerate(axes): 
                tmp = self.geometry.axis_get(axis)
                tmp.min_max_set(self.axes_k4c_min[i], \
                                self.axes_k4c_max[i], \
                                Hkl.UnitEnum.USER)
                self.geometry.axis_set(axis, tmp)
        elif self.geom==3:
            for i, axis in enumerate(axes): 
                tmp = self.geometry.axis_get(axis)
                tmp.min_max_set(self.axes_e6c_min[i], \
                                self.axes_e6c_max[i], \
                                Hkl.UnitEnum.USER)
                self.geometry.axis_set(axis, tmp)
        elif self.geom==4:
            for i, axis in enumerate(axes): 
                tmp = self.geometry.axis_get(axis)
                tmp.min_max_set(self.axes_k6c_min[i], \
                                self.axes_k6c_max[i], \
                                Hkl.UnitEnum.USER)
                self.geometry.axis_set(axis, tmp)

    def reset_pseudoaxes_solns(self):
        self.pseudoaxes_solns_h = 0
        self.pseudoaxes_solns_k = 0
        self.pseudoaxes_solns_l = 0
        self.pseudoaxes_solns_dpsacing = 0
        self.pseudoaxes_solns_fc2 = 0
        self.pseudoaxes_solns_psi = 0
        self.pseudoaxes_solns_q = 0
        self.pseudoaxes_solns_incidence = 0
        self.pseudoaxes_solns_azimuth1 = 0
        self.pseudoaxes_solns_emergence = 0
        self.pseudoaxes_solns_azimuth2 = 0

    def reset_axes_solns(self):
        if self.geom == 0 or 1:
            for lst in [self.axes_solns_omega_e4c, \
                        self.axes_solns_chi_e4c, \
                        self.axes_solns_phi_e4c, \
                        self.axes_solns_tth_e4c]:
                lst[:] = [0 for _ in lst]
        elif self.geom == 2:
            for lst in [self.axes_solns_komega_k4c, \
                        self.axes_solns_kappa_k4c, \
                        self.axes_solns_kphi_k4c, \
                        self.axes_solns_tth_k4c ]:
                lst[:] = [0 for _ in lst]
        elif self.geom == 3:
            for lst in [self.axes_solns_mu_e6c, \
                        self.axes_solns_omega_e6c, \
                        self.axes_solns_chi_e6c, \
                        self.axes_solns_phi_e6c, \
                        self.axes_solns_gamma_e6c, \
                        self.axes_solns_delta_e6c]:
                lst[:] = [0 for _ in lst]
        elif self.geom == 4:
            for lst in [self.axes_solns_mu_k6c, \
                        self.axes_solns_komega_k6c, \
                        self.axes_solns_kappa_k6c, \
                        self.axes_solns_kphi_k6c, \
                        self.axes_solns_gamma_k6c, \
                        self.axes_solns_delta_k6c]:
                lst[:] = [0 for _ in lst]

    def add_reflection1(self):
        '''
        adds reflection #1 to sample for busing-levy calculation
        '''
        print("add BL reflection 1")
        if self.geom == 0 or 1:
            for i in range(4):
                self.axes_UB_e4c[i] = self.refl1_input_e4c[(i+3)]
        elif self.geom == 2:
            for i in range(4):
                self.axes_UB_k4c[i] = self.refl1_input_k4c[(i+3)]
        elif self.geom == 3:
            for i in range(6):
                self.axes_UB_e6c[i] = self.refl1_input_e6c[(i+3)]
        elif self.geom == 4:
            for i in range(6):
                self.axes_UB_k6c[i] = self.refl1_input_k6c[(i+3)]
        self.forward_UB()
        if self.geom == 0 or 1:
            self.refl1 = self.sample.add_reflection(self.geometry, \
                                                    self.detector, \
                                                    self.refl1_input_e4c[0], \
                                                    self.refl1_input_e4c[1], \
                                                    self.refl1_input_e4c[2])
            self.refl_refine_input_list_e4c[self.curr_refl_index_e4c] = self.refl1_input_e4c.copy()
            self.refl_list_e4c.append(self.refl1)
            self.curr_refl_index_e4c += 1
        elif self.geom == 2:
            self.refl1 = self.sample.add_reflection(self.geometry, \
                                                    self.detector, \
                                                    self.refl1_input_k4c[0], \
                                                    self.refl1_input_k4c[1], \
                                                    self.refl1_input_k4c[2])
            self.refl_refine_input_list_k4c[self.curr_refl_index_k4c] = self.refl1_input_k4c.copy()
            self.refl_list_k4c.append(self.refl1)
            self.curr_refl_index_k4c += 1
        elif self.geom == 3:
            self.refl1 = self.sample.add_reflection(self.geometry, \
                                                    self.detector, \
                                                    self.refl1_input_e6c[0], \
                                                    self.refl1_input_e6c[1], \
                                                    self.refl1_input_e6c[2])
            self.refl_refine_input_list_e6c[self.curr_refl_index_e6c] = self.refl1_input_e6c.copy()
            self.refl_list_e6c.append(self.refl1)
            self.curr_refl_index_e6c += 1
        elif self.geom == 4:
            self.refl1 = self.sample.add_reflection(self.geometry, \
                                                    self.detector, \
                                                    self.refl1_input_k6c[0], \
                                                    self.refl1_input_k6c[1], \
                                                    self.refl1_input_k6c[2])
            self.refl_refine_input_list_k6c[self.curr_refl_index_k6c] = self.refl1_input_k6c.copy()
            self.refl_list_k6c.append(self.refl1)
            self.curr_refl_index_k6c += 1

    def add_reflection2(self):
        '''
        adds reflection #2 to sample for busing levy calculation
        '''
        print("add BL reflection 2")
        if self.geom == 0 or 1:
            for i in range(4):
                self.axes_UB_e4c[i] = self.refl2_input_e4c[(i+3)]
        elif self.geom == 2:
            for i in range(4):
                self.axes_UB_k4c[i] = self.refl2_input_k4c[(i+3)]
        elif self.geom == 3:
            for i in range(6):
                self.axes_UB_e6c[i] = self.refl2_input_e6c[(i+3)]
        elif self.geom == 4:
            for i in range(6):
                self.axes_UB_k6c[i] = self.refl2_input_k6c[(i+3)]
        self.forward_UB()
        # Hkl.SampleReflection(self.geometry, self.detector, h, k, l)
        if self.geom == 0 or 1:
            self.refl2 = self.sample.add_reflection(self.geometry, \
                                                    self.detector, \
                                                    self.refl2_input_e4c[0], \
                                                    self.refl2_input_e4c[1], \
                                                    self.refl2_input_e4c[2])
            self.refl_refine_input_list_e4c[self.curr_refl_index_e4c] = self.refl2_input_e4c.copy()
            self.refl_list_e4c.append(self.refl2)
            self.curr_refl_index_e4c += 1
        elif self.geom == 2:
            self.refl2 = self.sample.add_reflection(self.geometry, \
                                                    self.detector, \
                                                    self.refl2_input_k4c[0], \
                                                    self.refl2_input_k4c[1], \
                                                    self.refl2_input_k4c[2])
            self.refl_refine_input_list_k4c[self.curr_refl_index_k4c] = self.refl2_input_k4c.copy()
            self.refl_list_k4c.append(self.refl2)
            self.curr_refl_index_k4c += 1
        elif self.geom == 3:
            self.refl2 = self.sample.add_reflection(self.geometry, \
                                                    self.detector, \
                                                    self.refl2_input_e6c[0], \
                                                    self.refl2_input_e6c[1], \
                                                    self.refl2_input_e6c[2])
            self.refl_refine_input_list_e6c[self.curr_refl_index_e6c] = self.refl2_input_e6c.copy()
            self.refl_list_e6c.append(self.refl2)
            self.curr_refl_index_e6c += 1
        elif self.geom == 4:
            self.refl2 = self.sample.add_reflection(self.geometry, \
                                                    self.detector, \
                                                    self.refl2_input_k6c[0], \
                                                    self.refl2_input_k6c[1], \
                                                    self.refl2_input_k6c[2])
            self.refl_refine_input_list_k6c[self.curr_refl_index_k6c] = self.refl2_input_k6c.copy()
            self.refl_list_k6c.append(self.refl2)
            self.curr_refl_index_k6c += 1
 
    def compute_set_UB_matrix(self):
        #TODO replace by feeding into user input UB
        '''
        same thing as compute_UB_matrix, but without start()
        '''
        self.add_reflection1()
        self.add_reflection2()
        try:
            self.sample.compute_UB_busing_levy(self.refl1, self.refl2)
            UB = self.sample.UB_get()
            for i in range(3):
                for j in range(3):
                    self.UB_matrix[i,j] = UB.get(i,j)
            ct = datetime.datetime.now().isoformat()
            self.errors = [ord(c) for c in str(ct)] # success
        except Exception as e:
            ct = datetime.datetime.now().isoformat()
            self.errors = [ord(c) for c in str(ct + '\n' + str(e))] # failure
            print(f'compute_set_UB_matrix() error: {e}')

    def get_UB_matrix(self):
        UB = self.sample.UB_get()
        for i in range(3):
            for j in range(3):
                self.UB_matrix[i,j] = UB.get(i,j)

    def set_input_UB(self):
        UB_temp = self.sample.UB_get()
        for i in range(3):
            for j in range(3):
                print(UB_temp.get(i,j))
        self.UB_matrix = self.UB_matrix_input    
        Hkl.Matrix.init(UB_temp, *self.UB_matrix.ravel())         
        for i in range(3):
            for j in range(3):
                print(UB_temp.get(i,j))
        try:
            self.sample.UB_set(UB_temp)
            ct = datetime.datetime.now().isoformat()
            self.errors = [ord(c) for c in str(ct)] # success
        except Exception as e:
            ct = datetime.datetime.now().isoformat()
            self.errors = [ord(c) for c in str(ct + '\n' + str(e))] # failure
            print(f'set_input_UB() error: {e}')
        self.get_UB_matrix()

    def affine_set(self):
        '''
        takes in >2 reflections to refine lattice parameters and UB matrix
        '''
        try:
            self.sample.affine()
            #TODO just route PVs/python variables, don't re-run
            UB = self.sample.UB_get()
            for i in range(3):
                for j in range(3):
                    self.UB_matrix[i,j] = UB.get(i,j)
            self.latt[0], self.latt[1], self.latt[2], self.latt[3], self.latt[4], \
                self.latt[5] = self.lattice.get(Hkl.UnitEnum.USER)
            #self.latt = [self.lattice.get(Hkl.UnitEnum.USER)] #TODO test this
            # or maybe starred expression
        except Exception as e:
            ct = datetime.datetime.now().isoformat()
            self.errors = [ord(c) for c in str(ct + '\n' + str(e))]
            print(f'affine_set() error: {e}')
               
    def add_refl_refine(self):
        if self.geom == 0 or 1:
            for i in range(4):
                self.axes_UB_e4c[i] = self.refl_refine_input_e4c[(i+3)]
        elif self.geom == 2:
            for i in range(4):
                self.axes_UB_k4c[i] = self.refl_refine_input_k4c[(i+3)]
        elif self.geom == 3:
            for i in range(6):
                self.axes_UB_e6c[i] = self.refl_refine_input_e6c[(i+3)]
        elif self.geom == 4:
            for i in range(6):
                self.axes_UB_k6c[i] = self.refl_refine_input_k6c[(i+3)]
        self.forward_UB()
        if (self.geom == 0) or (self.geom==1):
            #TODO issue with self.refl_refine
            self.refl_refine = self.sample.add_reflection(self.geometry, \
                self.detector, self.refl_refine_input_e4c[0], \
                self.refl_refine_input_e4c[1], self.refl_refine_input_e4c[2])
            self.sample.add_reflection(self.geometry, \
                self.detector, self.refl_refine_input_e4c[0], \
                self.refl_refine_input_e4c[1], self.refl_refine_input_e4c[2])
            self.refl_refine_input_list_e4c[self.curr_refl_index_e4c] = self.refl_refine_input_e4c.copy()
            self.refl_list_e4c.append(self.refl_refine)
            self.curr_refl_index_e4c += 1
        elif self.geom == 2:
            self.refl_refine = self.sample.add_reflection(self.geometry, \
                self.detector, self.refl_refine_input_k4c[0], \
                self.refl_refine_input_k4c[1], self.refl_refine_input_k4c[2]) 
            self.refl_refine_input_list_k4c[self.curr_refl_index_k4c] = self.refl_refine_input_k4c.copy()
            self.refl_list_k4c.append(self.refl_refine)
            self.curr_refl_index_k4c += 1
        elif self.geom == 3:
            self.refl_refine = self.sample.add_reflection(self.geometry, \
                self.detector, self.refl_refine_input_e6c[0], \
                self.refl_refine_input_e6c[1], self.refl_refine_input_e6c[2])
            self.refl_refine_input_list_e6c[self.curr_refl_index_e6c] = self.refl_refine_input_e6c.copy()
            self.refl_list_e6c.append(self.refl_refine)
            self.curr_refl_index_e6c += 1
        elif self.geom == 4:
            self.refl_refine = self.sample.add_reflection(self.geometry, \
                self.detector, self.refl_refine_input_k6c[0], \
                self.refl_refine_input_k6c[1], self.refl_refine_input_k6c[2])
            self.refl_refine_input_list_k6c[self.curr_refl_index_k6c] = self.refl_refine_input_k6c.copy()
            self.refl_list_k6c.append(self.refl_refine)
            self.curr_refl_index_k6c += 1

    def del_refl_refine(self):
        i = self.selected_refl_index
        the_refl_list = self.sample.reflections_get()
        print(f'\nindex: {i}\n')
        #print(the_refl_list)
        #print(the_refl_list[i]) #todo
        if (self.geom==0) or (self.geom==1):
            print(self.refl_list_e4c[i])
            #self.sample.del_reflection(self.refl_list_e4c[i])
            self.sample.del_reflection(the_refl_list[i])
            self.refl_refine_input_list_e4c[i] = [0.,0.,0.,0.,0.,0.,0.]
            for j in range(i, 9):
                self.refl_refine_input_list_e4c[j] = self.refl_refine_input_list_e4c[j+1]
            self.refl_refine_input_list_e4c[9] = [0.,0.,0.,0.,0.,0.,0.]
            self.curr_refl_index_e4c -= 1

        elif self.geom==2:
            self.sample.del_reflection(self.refl_list_k4c[i])
            self.refl_refine_input_list_k4c[i] = [0.,0.,0.,0.,0.,0.,0.]
            for j in range(i, 9):
                self.refl_refine_input_list_k4c[j] = self.refl_refine_input_list_k4c[j+1]
            self.refl_refine_input_list_k4c[9] = [0.,0.,0.,0.,0.,0.,0.]
            self.curr_refl_index_k4c -= 1

        elif self.geom==3:
            self.sample.del_reflection(self.refl_list_e6c[i])
            self.refl_refine_input_list_e6c[i] = [0.,0.,0.,0.,0.,0.,0.,0.,0.]
            for j in range(i, 9):
                self.refl_refine_input_list_e6c[j] = self.refl_refine_input_list_e6c[j+1]
            self.refl_refine_input_list_e6c[9] = [0.,0.,0.,0.,0.,0.,0.,0.,0.]
            self.curr_refl_index_e6c -= 1

        elif self.geom==4:
            self.sample.del_reflection(self.refl_list_k6c[i])
            self.refl_refine_input_list_k6c[i] = [0.,0.,0.,0.,0.,0.,0.,0.,0.]
            for j in range(i, 9):
                self.refl_refine_input_list_k6c[j] = self.refl_refine_input_list_k6c[j+1]
            self.refl_refine_input_list_k6c[9] = [0.,0.,0.,0.,0.,0.,0.,0.,0.]
            self.curr_refl_index_k6c -= 1

    def reset_reflections(self):
        self.refl_refine_input_e4c = [0., 0., 0., 0., 0., 0., 0.]
        self.refl_refine_input_k4c = [0., 0., 0., 0., 0., 0., 0.]
        self.refl_refine_input_e6c = [0., 0., 0., 0., 0., 0., 0., 0., 0.]
        self.refl_refine_input_k6c = [0., 0., 0., 0., 0., 0., 0., 0., 0.]
        self.refl_refine_input_list_e4c = []
        self.refl_refine_input_list_k4c = []
        self.refl_refine_input_list_e6c = []
        self.refl_refine_input_list_k6c = []
        for i in range(self.num_reflections):
            self.refl_refine_input_list_e4c.append([0., 0., 0., 0., 0., 0., 0.])
            self.refl_refine_input_list_k4c.append([0., 0., 0., 0., 0., 0., 0.])
            self.refl_refine_input_list_e6c.append([0., 0., 0., 0., 0., 0., 0., 0., 0.])
            self.refl_refine_input_list_k6c.append([0., 0., 0., 0., 0., 0., 0., 0., 0.])

    def read_cif(self, givenpath):
        self.cif_path = givenpath

    def run_cif(self):
        #TODO better error handling within intensities function
        try:
            temp_intensities, templatt = intensities.intensity_calc(self.wavelength, self.cif_path)
        except Exception as e:
            self.errors = f'run_cif error (input file): {e}'
            return
        self.intensities = [ord(c) for c in str(temp_intensities)]
        try:
            self.latt[0] = templatt['a']
            self.latt[1] = templatt['b']
            self.latt[2] = templatt['c']
            self.latt[3] = templatt['alpha']
            self.latt[4] = templatt['beta']
            self.latt[5] = templatt['gamma']
            self.start()
        except Exception as e:
            self.errors = f'run_cif error (lattice parameters): {e}'
            return

    def get_sample_rotation(self):
        rot = self.geometry.sample_rotation_get(self.sample).to_matrix()
        dim = len(self.sample_rot_matrix)
        for i in range(dim):
            for j in range(dim):
                self.sample_rot_matrix[i,j] = rot.get(i,j)     
        return self.sample_rot_matrix  

    def get_u_matrix(self):
        rot = self.sample.U_get()
        dim = len(self.u_matrix)
        for i in range(dim):
            for j in range(dim):
                self.u_matrix[i,j] = rot.get(i,j)     
        return self.u_matrix

    def get_u_xyz(self):
        self.ux = self.sample.ux_get().value_get(0)
        self.uy = self.sample.uy_get().value_get(0)
        self.uz = self.sample.uz_get().value_get(0)
        return self.ux, self.uy, self.uz

    def set_wavelength(self, wlen):
        self.wavelength = wlen

    def get_latt_vol(self):
        self.lattice_vol = self.lattice.volume_get().value_get(0)
        print(self.lattice_vol)

    def get_info(self):
        lines = []
        #lines.append(f'self.geom: {self.geom}')
        diff_geom = self.factory.name_get()
        lines.append(str(diff_geom))
        x = self.engine_hkl.axis_names_get(Hkl.EngineAxisNamesGet.READ)
        lines.append(f'diffractometer axes:\n{x}')
        samp = self.sample.lattice_get()
        a, b, c, alpha, beta, gamma = samp.get(Hkl.UnitEnum.USER)
        lines.append(f'a: {a}, b: {b}, c: {c}, alpha: {alpha}, beta: {beta}, gamma: {gamma}')
        lines.append(f'UB: {self.UB_matrix}')
        lines.append(f'latt vol: {self.lattice_vol}')
        return '\n'.join(lines)
