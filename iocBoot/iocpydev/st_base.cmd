#!../../bin/linux-x86_64/hklApp

< envPaths

# PYTHONPATH points to folders where Python modules are.
epicsEnvSet("PYTHONPATH","$(TOP)/python")

# Designate beamline:IOC_name as conventional prefix
epicsEnvSet("PREFIX", "abtest:ioc-hkl:")

cd ${TOP}

## Register all support components
dbLoadDatabase "${TOP}/dbd/hklApp.dbd"
hklApp_registerRecordDeviceDriver pdbbase

## Load record instances
dbLoadRecords("db/hkl_main.db")

cd ${TOP}/iocBoot/${IOC}

pydev("import hklApp")
pydev("hkl_calc = hklApp.hklCalcs()")

iocInit

epicsThreadSleep(1)

dbpf("$(PREFIX)geom","1")
dbpf("$(PREFIX)wlen","5.431")
dbpf("$(PREFIX)latt_a","5.431")
dbpf("$(PREFIX)latt_b","5.431")
dbpf("$(PREFIX)latt_c","5.431")
dbpf("$(PREFIX)latt_alpha","90.")
dbpf("$(PREFIX)latt_beta","90.")
dbpf("$(PREFIX)latt_gamma","90.")
dbpf("$(PREFIX)omega_e4c","30.")
dbpf("$(PREFIX)chi_e4c","20.")
dbpf("$(PREFIX)phi_e4c","10.")
dbpf("$(PREFIX)tth_e4c","10.")
dbpf("$(PREFIX)refl1_h_e4c","0.")
dbpf("$(PREFIX)refl1_k_e4c","0.")
dbpf("$(PREFIX)refl1_l_e4c","4.")
dbpf("$(PREFIX)refl1_omega_e4c","-145.")
dbpf("$(PREFIX)refl1_chi_e4c","0.")
dbpf("$(PREFIX)refl1_phi_e4c","0.")
dbpf("$(PREFIX)refl1_tth_e4c","60.")
dbpf("$(PREFIX)refl2_h_e4c","0.")
dbpf("$(PREFIX)refl2_k_e4c","4.")
dbpf("$(PREFIX)refl2_l_e4c","0.")
dbpf("$(PREFIX)refl2_omega_e4c","-145.")
dbpf("$(PREFIX)refl2_chi_e4c","90.")
dbpf("$(PREFIX)refl2_phi_e4c","0.")
dbpf("$(PREFIX)refl2_tth_e4c","69.")
dbpf("$(PREFIX)omega_e4c_min","-180.")
dbpf("$(PREFIX)chi_e4c_min","-180.")
dbpf("$(PREFIX)phi_e4c_min","-180.")
dbpf("$(PREFIX)tth_e4c_min","-180.")
dbpf("$(PREFIX)omega_e4c_max","180.")
dbpf("$(PREFIX)chi_e4c_max","180.")
dbpf("$(PREFIX)phi_e4c_max","180.")
dbpf("$(PREFIX)tth_e4c_max","180.")
dbpf("$(PREFIX)komega_k4c_min","-180.")
dbpf("$(PREFIX)kappa_k4c_min","-180.")
dbpf("$(PREFIX)kphi_k4c_min","-180.")
dbpf("$(PREFIX)tth_k4c_min","-180.")
dbpf("$(PREFIX)komega_k4c_max","180.")
dbpf("$(PREFIX)kappa_k4c_max","180.")
dbpf("$(PREFIX)kphi_k4c_max","180.")
dbpf("$(PREFIX)tth_k4c_max","180.")
dbpf("$(PREFIX)mu_e6c_min","-180.")
dbpf("$(PREFIX)omega_e6c_min","-180.")
dbpf("$(PREFIX)chi_e6c_min","-180.")
dbpf("$(PREFIX)phi_e6c_min","-180.")
dbpf("$(PREFIX)gamma_e6c_min","-180.")
dbpf("$(PREFIX)delta_e6c_min","-180.")
dbpf("$(PREFIX)mu_e6c_max","180.")
dbpf("$(PREFIX)omega_e6c_max","180.")
dbpf("$(PREFIX)chi_e6c_max","180.")
dbpf("$(PREFIX)phi_e6c_max","180.")
dbpf("$(PREFIX)gamma_e6c_max","180.")
dbpf("$(PREFIX)delta_e6c_max","180.")
dbpf("$(PREFIX)mu_k6c_min","-180.")
dbpf("$(PREFIX)komega_k6c_min","-180.")
dbpf("$(PREFIX)kappa_k6c_min","-180.")
dbpf("$(PREFIX)kphi_k6c_min","-180.")
dbpf("$(PREFIX)gamma_k6c_min","-180.")
dbpf("$(PREFIX)delta_k6c_min","-180.")
dbpf("$(PREFIX)mu_k6c_max","180.")
dbpf("$(PREFIX)komega_k6c_max","180.")
dbpf("$(PREFIX)kappa_k6c_max","180.")
dbpf("$(PREFIX)kphi_k6c_max","180.")
dbpf("$(PREFIX)gamma_k6c_max","180.")
dbpf("$(PREFIX)delta_k6c_max","180.")
dbpf("$(PREFIX)refl_h_si_e4c","4.")
dbpf("$(PREFIX)refl_omega_si_e4c","-145.")
dbpf("$(PREFIX)refl_phi_si_e4c","90.")
dbpf("$(PREFIX)refl_tth_si_e4c","65.")
dbpf("$(PREFIX)pseudoaxes_psi_h2","1.")
dbpf("$(PREFIX)pseudoaxes_psi_k2","1.")
dbpf("$(PREFIX)pseudoaxes_psi_l2","1.")
dbpf("$(PREFIX)pseudoaxes_incidence_x","0.")
dbpf("$(PREFIX)pseudoaxes_incidence_y","1.")
dbpf("$(PREFIX)pseudoaxes_incidence_z","0.")
dbpf("$(PREFIX)pseudoaxes_emergence_x","0.")
dbpf("$(PREFIX)pseudoaxes_emergence_y","1.")
dbpf("$(PREFIX)pseudoaxes_emergence_z","0.")
dbpf("$(PREFIX)pseudoaxes_eulerians_solns","1.")
dbpf("$(PREFIX)pseudoaxes_qperqpar_x","0.")
dbpf("$(PREFIX)pseudoaxes_qperqpar_y","1.")
dbpf("$(PREFIX)pseudoaxes_qperqpar_z","0.")



dbl > pvlist.dbl
