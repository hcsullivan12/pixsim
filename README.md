# pixsim - Field response calculation for pixel geometries using BEM
# Overview
This tool was inspired by Brett Viren's tool [larf](https://github.com/brettviren/larf) for wire LArTPCs.

What this package can do:
- Generate 2D/3D geometries using high level information. Has dedicated classes for different TPC geometries. Explicit construction carried out by [pygmsh](https://link). 
- Load surface meshes and corresponding domain mapping, apply user-defined Dirichlet boundary conditions, and solve for Nuemann boundary conditions using BEM++.
- Evalute electric potential on a volume of points.
- Plot potential and electric field.
- Drift charge through electric field.
- Compute weighting fields.
- Compute instantaneous current seen seen for different drift paths.
- Export results into *.vtk* format for visualization with [Paraview](https://link).

This package is a command-line tool that has seperate commands for each stage of the simulation. Data from each result is persisted in a database that gets updated at the end of each stage and can referenced by name or ID in subsequent stages. 

# Geometry
Currently, [pixsim](https://link) uses the [pygmsh](https://link) api to generate geometries and the corresponding *.geo* files. The readout electrodes are defined at a high level (e.g. position, shape) and are used in constructing the TPC geometry. In principle, any geometry with arbitrary complexity can be constructed, and the [pygmsh](https://link) api provides a relatively simple way of building these geometries. While [larf](https://github.com/brettviren/larf) has dedicated routines for constructing the meshes, this layer has been removed in pixsim. Rather, meshes are constructed from the *.geo* files using [GMSH](https://link). The high level information for the readout electrodes (pixels) is saved and the *.msh* file can be passed at runtime for subsequent algorithms. 

There are a few naming conventions throughout the simulation. The geometry construction has a few hard-coded names for the domains that are assumed when solving for the boundary conditions. 