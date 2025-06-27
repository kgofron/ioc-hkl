import numpy as np


def businglevy(wlen, lattice, refl1, refl2):
    omega = refl1['omega']
    tth = refl1['tth']
    
    def scattering_unit_vec():
    return UB

def simplex(wlen, lattice, refl_list):
    return lattice, UB

def motor2hkl(wlen, lattice, omega, tth):
    return h, k, l

def hkl2motor(wlen, lattice, h, k, l):
    return omega, tth



def B_reciplatt_calc(a,b,c,alpha,beta,gamma):
    V_abg = np.sqrt(1 - np.cos^2(alpha) - np.cos^(beta) - np.cos^2(gamma) + \
                    2*np.cos(alpha)*np.cos(beta)*np.cos(gamma)) # V_alphabetagamma

    a_star = np.sin(alpha)/(V_abg*a)
    b_star = np.sin(beta)/(V_abg*b)
    c_star = np.sin(gamma)/V_abg*c)
    alpha_star = np.arccos((np.cos(beta)*np.cos(gamma)-np.cos(alpha))/(np.sin(beta)*np.sin(gamma))
    beta_star  = np.arccos((np.cos(gamma)*np.cos(alpha)-np.cos(beta))/(np.sin(gamma)*np.sin(alpha))
    gamma_star = np.arccos((np.cos(alpha)*np.cos(beta)-np.cos(gamma))/(np.sin(alpha)*np.sin(beta))

    B = np.ndarray([[a_star, 0., 0.],
                    [b_star*np.cos(gamma_star), b_star*np.sin(gamma_star), 0.],
                    [c_star*np.cos(beta_star), -c_star*sin(beta_star)*np.cos(alpha), 1/c])

    recip_latt =  a_star, b_star, c_star, alpha_star, beta_star, gamma_star

    return B, recip_latt



# Example known reflections
reflections = [
    {"hkl": [1, 0, 0], "theta": 20.0, "twotheta": 40.0},  # in degrees
    {"hkl": [0, 1, 0], "theta": 30.0, "twotheta": 60.0},
]
wavelength = 1.5406  # Angstrom, for Cu Kα

# Calculate B matrix (reciprocal lattice vectors)
import numpy as np
def b_matrix_from_lattice(a, b, c, alpha, beta, gamma):
    alpha, beta, gamma = np.radians([alpha, beta, gamma])
    volume = (
        a * b * c
        * np.sqrt(
            1
            - np.cos(alpha)**2
            - np.cos(beta)**2
            - np.cos(gamma)**2
            + 2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma)
        )
    )
    B = np.zeros((3, 3))
    B[0, 0] = 2 * np.pi / a
    B[0, 1] = 0
    B[0, 2] = 0
    B[1, 0] = -2 * np.pi * np.cos(gamma) / (a * np.sin(gamma))
    B[1, 1] = 2 * np.pi / (b * np.sin(gamma))
    B[1, 2] = 0
    B[2, 0] = 2 * np.pi * (
        np.cos(alpha) * np.cos(gamma) - np.cos(beta)
    ) / (a * volume * np.sin(gamma))
    B[2, 1] = -2 * np.pi * np.cos(alpha) / (b * volume * np.sin(gamma))
    B[2, 2] = 2 * np.pi / (c * volume)
    return B


#calculate scattered vectors \vec{Q} from each reflection
# |\vec{Q}| = \frac{4\pi}{lambda}sin(\theta)
# \vec{Q} = |\vec{Q}| \dot \heat{k_f} - \hat{k_i}
def q_vector(theta_deg, twotheta_deg, wavelength):
    theta = np.radians(theta_deg)
    twotheta = np.radians(twotheta_deg)

    ki = np.array([np.cos(theta), 0, -np.sin(theta)])
    kf = np.array([np.cos(twotheta - theta), 0, np.sin(twotheta - theta)])

    return (2 * np.pi / wavelength) * (kf - ki)


# example usage
G = []  # hkl
Q = []  # measured Q

for refl in reflections:
    G.append(refl["hkl"])
    Q.append(q_vector(refl["theta"], refl["twotheta"], wavelength))

G = np.array(G).T  # 3x2
Q = np.array(Q).T  # 3x2
UB = Q @ np.linalg.inv(G)  # Solve Q = UB * hkl

def hkl_to_q(hkl, UB):
    return UB @ np.array(hkl)


