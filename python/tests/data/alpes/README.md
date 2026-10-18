# Open ALPES geometry fixture

Attribution: **The ALPES Team**, [ALPES](https://github.com/itpplasma/ALPES),
revision `eb31dc8f4725a5bc355a558b4357f9234c312264`.
These data are licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
The complete upstream license is retained in `LICENSE-CC-BY-4.0`.
See the upstream README for the attribution and data licensing statement.

Unmodified files from that revision (about 225 KB total):

- `wout.nc`: `design/CoilOptimization/wout_0021326_1_000_000043.nc`, renamed only.
- `coil_coordinates0.txt` through `coil_coordinates11.txt`:
  `coils/coils_xyz/coilData/coil_coordinates*.txt`, unchanged XYZ filament samples.

`SHA256SUMS` records the vendored data checksums. No network access is needed
for tests. To verify the copies, run `sha256sum -c SHA256SUMS` in this directory.
Pinned original files are available under
`https://raw.githubusercontent.com/itpplasma/ALPES/eb31dc8f4725a5bc355a558b4357f9234c312264/`
followed by the upstream paths above.

Tests use the numerical coordinates as supplied, without applying the separate
engineering scaling in ALPES scripts. The wout major radius is about 0.99 in
these units. Unit currents are assigned only to construct libneo Filaments;
these tests do not compute a magnetic field or assert that these separate
upstream coil and equilibrium files are a self-consistent field solution.
The wall offset is 0.05 in the wout units, measured in each RZ plane. Winding
packs use 0.01 by 0.02 in the filament coordinate units. Collision detection
uses GECKO's existing approximate surface-proximity algorithm, not an exact
solid-intersection oracle. A translated duplicate exercises positive detection.
