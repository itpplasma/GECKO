# STEP round-trip source

`alpes_coil0.txt` is the complete (160-point) `coil_coordinates0.txt` from
[itpplasma/ALPES](https://github.com/itpplasma/ALPES/blob/eb31dc8f4725a5bc355a558b4357f9234c312264/coils/coils_xyz/coilData/coil_coordinates0.txt),
revision `eb31dc8f4725a5bc355a558b4357f9234c312264`, shared under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Attribution:
ALPES contributors. The data is unchanged; tests multiply coordinates by 330
(the generator's `100 * 0.330 * 10`) to obtain its millimeter model.

`step_helpers.loft_coil` adapts the equations and sequential rectangular
loft/union construction in ALPES `windingCoordinateGenerator.py` and
`CadQueryMultiCoilGenerator.py` at that revision, using coil 0's Y steering
vector and the original 40 by 67 mm section. It exports only to pytest's
temporary directory, and does not execute the original export script.
Those source files are MIT licensed; the notice is reproduced in
`ALPES-LICENSE.txt`. STEP files are generated at test time, not committed.
