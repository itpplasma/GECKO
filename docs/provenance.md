# STL provenance

`gecko wall-from-vmec`, `gecko stl-convert` and
`gecko.stl.export_mesh_to_stl(mesh, path)` write a sidecar file
`<filename>.stl.json` next to each STL file. It records the GECKO package
version, the output of `git describe --tags --always --dirty` and the SHA-256
checksum of the STL file. Keep the sidecar with meshes used in publications.

Installed wheels keep the revision they were built from; source archives built
without Git record a null revision. Writing a mesh directly with a third-party
call such as `mesh.export(...)` does not create a sidecar.
