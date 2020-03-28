#!/usr/bin/env python3

import subprocess
import sys
import os.path

directory = os.path.dirname(__file__)

tests = []

try:
    import pyflatsurf
    tests += ["pyflatsurf_conversion.py"]
    tests += [os.path.join("gl2r_orbit_closure", name) for name in ["discriminat_loci_H_1_1.py", "gothic_locus.py",
              "rank2_quadrilaterals.py",  "veech_n_gons_and_double_n_gons.py", "veech_surfaces_H2.py"]]
except ImportError:
    pass

output = 0
for test in tests:
    test = os.path.join(directory, test)
    output += subprocess.call([sys.executable, test] + sys.argv[1:])

sys.exit(output)
