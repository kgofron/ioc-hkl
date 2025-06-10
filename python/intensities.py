from org.csstudio.display.builder.runtime.script import ScriptUtil
from java.lang import Runtime
from javax.swing import JFileChooser
import java.io.BufferedReader as BufferedReader
import java.io.InputStreamReader as InputStreamReader
#import pyepics
import subprocess
from org.csstudio.display.builder.runtime.script import PVUtil



#print(dir(ScriptUtil))

#pv = ScriptUtil.getPV("loc://abtest:ioc-hkl:intensities")
#pv.write("Waiting for intensities...")


#TODO For now, I'm pasting this script into the "import sample .cif" action script in CSS

#TODO get wavelength from IOC with pyepics, feed into cif2hkl, for now
wavelength=0.4

cif2hkl_bin = '/usr/bin/cif2hkl'


chooser = JFileChooser()
chooser.setDialogTitle("Select a CIF file")
result = chooser.showOpenDialog(None)

if result == JFileChooser.APPROVE_OPTION:
    file = chooser.getSelectedFile()
    cif_path = file.getAbsolutePath()
    hkl_path = cif_path.replace(".cif", ".hkl")

    print(cif_path)
    print(hkl_path)
    #cmd = [cif2hkl_bin, '--mode', 'NUC', '--lambda', str(wavelength), '--xtal', cif_path, '--out', hkl_path]
    cmd = [cif2hkl_bin, '--mode', 'NUC', '--out', hkl_path, '--lambda', str(wavelength), '--xtal', cif_path]
    #cmd = [cif2hkl_bin, '--xtal', cif_path]
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

            # Wait until we hit the reflection header
            if line.startswith("# H") and "|Fc|^2" in line:
                found_data_start = True
                continue

            # Skip all lines before the reflection table
            if not found_data_start or line.startswith("#"):
                continue

            parts = line.split()
            if len(parts) >= 6:
                h, k, l, mult, d, intensity = parts[:6]
                intensity_lines.append("%3s %3s %3s   %12s" % (h, k, l, intensity))

            if len(intensity_lines) >= 50:
                break

    except Exception as e:
        intensity_lines = ["Error: " + str(e)]
        
    output = "\n".join(intensity_lines)
    print("Parsed reflections:\n", output)

    # Push to local PV for display
    #ScriptUtil.getPVByName(widget, "loc://abtest:ioc-hkl:intensities").write(output)
    PVUtil.writePV("loc://abtest:ioc-hkl:intensities", output, 1000)

