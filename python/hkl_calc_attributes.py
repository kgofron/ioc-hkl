#TODO

# Then in the class: 

# my_class.py
#from hkl_calc_attributes import hkl_calc_defaults

#class hkl_calc:
#    def __init__(self):
#        for k, v in hkl_calc_defaults.items():
#            setattr(self, k, v)



hkl_calc_defaults={
        # initials
        "wavelength": 0.,
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
        self.businglevyflag = 0
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

        self.refl_list_e4c = [] #TODO maybe delete - depends if the list returned and stored in the_refl_list is in the correct order everytime it is received - store reflection objects, for deleting
        self.refl_list_k4c = [] #TODO maybe delete store reflection objects, for deleting
        self.refl_list_e6c = [] #TODO maybe delete store reflection objects, for deleting
        self.refl_list_k6c = [] #TODO maybe delete store reflection objects, for deleting
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
}
