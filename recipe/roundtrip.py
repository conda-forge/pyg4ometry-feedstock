"""Build a CSG geometry, write it to GDML and read it back.

The boolean and the mesh info come from the compiled CGAL bindings, the reader
parses the GDML expressions with the ANTLR runtime, so a clean run covers the
parts of the package that a bare import does not.
"""

from __future__ import annotations

import math
import pathlib

import pyg4ometry.gdml as gd
import pyg4ometry.geant4 as g4

reg = g4.Registry()

wx = gd.Constant("wx", "200", reg, True)
bx = gd.Constant("bx", "100", reg, True)
tr = gd.Constant("tr", "30", reg, True)

world_s = g4.solid.Box("world_s", wx, wx, wx, reg, "mm")
box_s = g4.solid.Box("box_s", bx, bx, bx, reg, "mm")
tub_s = g4.solid.Tubs("tub_s", 0, tr, 2 * bx, 0, "2*pi", reg, "mm")
sub_s = g4.solid.Subtraction("sub_s", box_s, tub_s, [[0, 0, 0], [0, 0, 0]], reg)

world_l = g4.LogicalVolume(world_s, g4.MaterialPredefined("G4_Galactic"), "world_l", reg)
sub_l = g4.LogicalVolume(sub_s, g4.MaterialPredefined("G4_Fe"), "sub_l", reg)
g4.PhysicalVolume([0, 0, 0], [0, 0, 0], sub_l, "sub_pv", world_l, reg)
reg.setWorld(world_l.name)

info = sub_l.mesh.localmesh.info()
assert not info["null"] and info["closed"], info
assert not info["selfintersect"], info

# the mesher approximates the cylinder with a 16-sided prism, so compare
# against the analytic volume loosely -- a boolean that silently did nothing
# would still be 39% out
expected = 100**3 - math.pi * 30**2 * 100
assert math.isclose(info["volume"], expected, rel_tol=0.03), (info["volume"], expected)

out = pathlib.Path("roundtrip.gdml")
writer = gd.Writer()
writer.addDetector(reg)
writer.write(out)
assert out.stat().st_size > 0

reg2 = gd.Reader(str(out)).getRegistry()
assert reg2.getWorldVolume().name == "world_l"
assert {"world_s", "box_s", "tub_s", "sub_s"} <= set(reg2.solidDict), sorted(reg2.solidDict)

info2 = reg2.logicalVolumeDict["sub_l"].mesh.localmesh.info()
assert math.isclose(info2["volume"], info["volume"], rel_tol=1e-9), (
    info2["volume"],
    info["volume"],
)

print("round-trip OK, mesh volume =", info["volume"])
