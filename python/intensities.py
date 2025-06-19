import subprocess
import re

def parse_cif_lattice_params(cif_path):
    lattice = {}
    pattern = re.compile(r"([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)(?:\(\d+\))?")
    with open(cif_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('_cell_length_a'):
                match = pattern.search(line)
                if match:
                    lattice['a'] = float(match.group(1))
            elif line.startswith('_cell_length_b'):
                match = pattern.search(line)
                if match:
                    lattice['b'] = float(match.group(1))
            elif line.startswith('_cell_length_c'):
                match = pattern.search(line)
                if match:
                    lattice['c'] = float(match.group(1))
            elif line.startswith('_cell_angle_alpha'):
                match = pattern.search(line)
                if match:
                    lattice['alpha'] = float(match.group(1))
            elif line.startswith('_cell_angle_beta'):
                match = pattern.search(line)
                if match:
                    lattice['beta'] = float(match.group(1))
            elif line.startswith('_cell_angle_gamma'):
                match = pattern.search(line)
                if match:
                    lattice['gamma'] = float(match.group(1))
    return lattice


def intensity_calc(wavelength, cif_path):
    hkl_path = cif_path.replace(".cif", ".hkl")
    print(cif_path)
    print(hkl_path)
    lattice = parse_cif_lattice_params(cif_path)
    print(lattice)
    ##### scattering intensities calculation #####
    cif2hkl_bin = '/usr/bin/cif2hkl'
    cmd = [cif2hkl_bin, '--mode', 'NUC', '--out', hkl_path, '--lambda', str(wavelength), '--xtal', cif_path]
    print(cmd)
    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = proc.communicate()
        #TODO propogate error/success text
        print(out)
        print(err)
    except Exception as e:
        print("Exception running cif2hkl:", e)
        output = e
        return output, lattice
    try:
        with open(hkl_path, "r") as f:
            lines = f.readlines()
        intensity_lines = []
        found_data_start = False
        for line in lines:
            line = line.strip()
            if not line:
                continue
            if line.startswith("# H") and "|Fc|^2" in line:
                found_data_start = True
                continue
            if not found_data_start or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 6:
                h, k, l, mult, d, intensity = parts[:6]
                intensity_lines.append("%3s %3s %3s %8s %12s" % (h, k, l, d, intensity)) 
            if len(intensity_lines) >= 50:
                break
    except Exception as e:
        intensity_lines = ["Error: " + str(e)]
    output = "\n".join(intensity_lines)
    print(output)
    print(lattice)
    return output, lattice
