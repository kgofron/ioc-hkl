#TODO move from CSS PV updater to python (pyepics)
from org.csstudio.display.builder.runtime.script import ScriptUtil
from java.lang import Runtime
from javax.swing import JFileChooser
import java.io.BufferedReader as BufferedReader
import java.io.InputStreamReader as InputStreamReader
import subprocess
from org.csstudio.display.builder.runtime.script import PVUtil
import re

prefix = "abtest:ioc-hkl:"

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

# Embedded python script
#from org.csstudio.display.builder.runtime.script import PVUtil, ScriptUtil
#print 'Hello'
# widget.setPropertyValue('text', PVUtil.getString(pvs[0]))

#print(dir(ScriptUtil))
#print(dir(PVUtil))

#TODO For now, I'm pasting this script into the "import sample .cif" action script in CSS

chooser = JFileChooser()
chooser.setDialogTitle("Select a CIF file")
result = chooser.showOpenDialog(None)

if result == JFileChooser.APPROVE_OPTION:
    file = chooser.getSelectedFile()
    cif_path = file.getAbsolutePath()
    hkl_path = cif_path.replace(".cif", ".hkl")

    print(cif_path)
    print(hkl_path)

    ##### parse .cif and set lattice parameters to IOC #####
    lattice = parse_cif_lattice_params(cif_path)
    PVUtil.writePV(f"{prefix}latt_a", lattice['a'], 1000)
    PVUtil.writePV(f"{prefix}latt_b", lattice['b'], 1000)
    PVUtil.writePV(f"{prefix}latt_c", lattice['c'], 1000)
    PVUtil.writePV(f"{prefix}latt_alpha", lattice['alpha'], 1000)
    PVUtil.writePV(f"{prefix}latt_beta", lattice['beta'], 1000)
    PVUtil.writePV(f"{prefix}latt_gamma", lattice['gamma'], 1000)

    ##### scattering intensities calculation #####
    cif2hkl_bin = '/usr/bin/cif2hkl'
    wlenpv = PVUtil.createPV(f"{prefix}wlen_RBV", 1000)
    wavelength = PVUtil.getDouble(wlenpv)
    #TODO if wavelength value is entered in box but "Set Wavelength & Sample Lattice" isn't pressed and assigned to the underlying hkl smaple, the resulting intensities from cif2hkl will not be correct
    cmd = [cif2hkl_bin, '--mode', 'NUC', '--out', hkl_path, '--lambda', str(wavelength), '--xtal', cif_path]
    #cmd = [cif2hkl_bin, '--xtal', cif_path]
    print("cmd:", cmd)
    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = proc.communicate()
        print("cif2hkl stdout:\n", out)
        print("cif2hkl stderr:\n", err)
    except Exception as e:
        print("Exception running cif2hkl:", e)
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
                #intensity_lines.append("%3s %3s %3s   %12s" % (h, k, l, intensity))
            if len(intensity_lines) >= 50:
                break
    except Exception as e:
        intensity_lines = ["Error: " + str(e)]
    output = "\n".join(intensity_lines)
    #print("Parsed reflections:\n", output)
    PVUtil.writePV(f"{prefix}intensities", output, 1000)
