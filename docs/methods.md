# Methods and Algorithms

This document describes the algorithms and methods used in GECKO, with references
to the underlying literature for reproducibility and citation purposes.

## Coil Tube Mesh Generation

### Overview

GECKO generates 3D tube meshes from coil centerlines with rectangular cross-sections.
This is essential for visualizing stellarator and tokamak coil winding packs.

### Centroid Frame

The key challenge in tube mesh extrusion is defining an orthonormal frame
$(\mathbf{t}, \mathbf{n}, \mathbf{b})$ at each point along the curve that:

1. Varies smoothly along the curve
2. Closes seamlessly for closed loops (coils)
3. Has no twist accumulation

GECKO uses the **centroid frame** approach from {cite:t}`singh2020coils`, originally
developed for stellarator coil optimization. This defines the frame as:

$$
\mathbf{t} = \frac{d\mathbf{r}/ds}{|d\mathbf{r}/ds|}
$$

$$
\mathbf{n} = \frac{\boldsymbol{\delta} - (\boldsymbol{\delta} \cdot \mathbf{t})\mathbf{t}}{|\boldsymbol{\delta} - (\boldsymbol{\delta} \cdot \mathbf{t})\mathbf{t}|}
$$

$$
\mathbf{b} = \mathbf{t} \times \mathbf{n}
$$

where:
- $\mathbf{t}$ is the unit tangent vector
- $\mathbf{n}$ is the normal vector pointing toward the coil centroid
- $\mathbf{b}$ is the binormal completing the right-handed frame
- $\boldsymbol{\delta} = \mathbf{R} - \mathbf{r}$ points from the curve point toward the centroid $\mathbf{R}$

### Why Centroid Frame?

Alternative frame definitions have limitations:

| Method | Issue for Closed Loops |
|--------|----------------------|
| **Frenet frame** | Undefined at inflection points; twist accumulation |
| **Global up vector** | Discontinuous jumps when tangent aligns with up |
| **Rotation-minimizing frame** | Numerical drift; residual twist at closure |
| **Centroid frame** | Natural closure; no twist accumulation |

The centroid frame works because all points reference the same centroid,
ensuring the frame at the last point naturally aligns with the first point.

### Implementation

```python
from gecko.coils import coil_to_tube_mesh

# Create tube mesh with 40cm x 80cm rectangular cross-section
mesh = coil_to_tube_mesh(coords, width=0.4, height=0.8, closed=True)
```

See {py:func}`gecko.coils.coil_to_tube_mesh` for the full API.

## Wall Mesh Generation

### VMEC Surface Extraction

GECKO extracts flux surfaces from VMEC equilibria using the Fourier representation:

$$
R(\theta, \zeta) = \sum_{m,n} R_{mn} \cos(m\theta - n N_p \zeta)
$$

$$
Z(\theta, \zeta) = \sum_{m,n} Z_{mn} \sin(m\theta - n N_p \zeta)
$$

where $N_p$ is the number of field periods.

The surface is triangulated using a structured grid in $(\theta, \zeta)$ space.

### Wall Offset

For creating first-wall geometries offset from the LCFS, GECKO uses polygon buffering
via Shapely {cite:p}`shapely`. In each poloidal plane ($\zeta = \text{const}$):

1. Extract the $(R, Z)$ cross-section
2. Apply Shapely's `buffer()` operation with the desired offset
3. Triangulate the resulting surface

This approach handles:
- Convex and concave surfaces
- Self-intersecting offset curves (automatically resolved)
- Non-uniform curvature

## Visualization

GECKO uses PyVista {cite:p}`pyvista` for interactive 3D visualization, providing:

- GPU-accelerated rendering via VTK
- Interactive rotation, zoom, and pan
- Off-screen rendering for headless servers
- Export to screenshots and various 3D formats

## References

```{bibliography}
:filter: docname in docnames
```
