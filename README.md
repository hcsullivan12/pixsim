# pixsim - Field response calculation for pixel geometries using BEM
# Overview
This tool was inspired by Brett Viren's tool [larf](https://github.com/brettviren/larf) for wire LArTPCs. 

What this package can do:
- Generate 2D/3D geometries using high level information. Has dedicated classes for different TPC geometries. Explicit construction carried out by [pygmsh](https://pypi.org/project/pygmsh/). 
- Load surface meshes and corresponding domain mapping, apply user-defined Dirichlet boundary conditions, and solve for Nuemann boundary conditions using BEM++.
- Evalute electric potential on a volume of points.
- Plot potential and electric field.
- Drift charge through electric field.
- Compute weighting fields.
- Compute instantaneous current on electrodes for different drift paths.
- Export results into *.vtk* format for visualization with [ParaView](https://www.paraview.org/).

This package is a command-line tool that has seperate commands for each stage of the simulation. Data from each result is persisted in a database that gets updated at the end of each stage and can be referenced by name or ID in subsequent stages. 

See [install](README/INSTALL.md) for installation instructions.
See [driftsim](README/DRIFT_SIM.md) for details on integrating with a drift simulation.

![alt text](https://raw.githubusercontent.com/hcsullivan12/picture_store/master/steps2.png "Steps")


# Geometry
Currently, *pixsim* uses the [pygmsh](https://pypi.org/project/pygmsh/) api to generate geometries and the corresponding *.geo* files. The readout electrodes are defined at a high level (e.g. position, shape) and are used in constructing the TPC geometry. In principle, any geometry with arbitrary complexity can be constructed, and the [pygmsh](https://pypi.org/project/pygmsh/) api provides a relatively simple way of building these geometries. While [larf](https://github.com/brettviren/larf) has dedicated routines for constructing the meshes, this layer has been removed in pixsim. Rather, meshes are constructed from the *.geo* files using [GMSH](http://gmsh.info/). The high level information for the readout electrodes (pixels) is saved and the *.msh* file can be passed at runtime for subsequent algorithms. 

There are a few naming conventions throughout the simulation. The geometry construction has a few hard-coded names for the domains that are assumed when solving for the boundary conditions. 

# Example
This section will walk through a simulation using the *boxtpc3/sphere_pads/with_grid* configuration and will use the alias `pix` defined as
```
$ alias pix='pixsim -c /path/to/boxtpc.cfg -s test.db -m tpcgeometry.msh'
```
The configuration file *boxtpc.cfg* defines the parameters to use in the various stages of the simulation. See *drift_sim/with_grid/boxtpc.cfg* as an example. The results will be saved in a database *test.db*. Although not every algorithm needs it, the mesh file *tpcgeometry.msh* is passed in using the `-m` flag. We will generate this file shortly. Each command requires the path to the configuration file, database, and mesh file.

To view a complete list of available commands, run 
```
$ pix --help
```

## Geometry
There should already be a geometry and mesh file. If you want to create a new geometry/mesh, you can add or change the implentation in *pixsim/geometry.py*, then run
```
$ pix gen geo -c <geometry_section_name> -o <name_of_geo_file>
```
Many of the commands offer the `-c` flag as an option to specify the section name in the configuration file and should have default names if not used. Although this command was meant to generate a *.geo* file, it will also generate a *.msh* file. The *.msh* can be used directly or the *.geo* file can be imported into [GMSH](http://gmsh.info/) where a new mesh can be constructed. 

The geometry construction should have also uniquely defined the domains. To save this in our store, run
```
$ pix gen dmap
```

## Boundary
Next stage is to solve for the Neumann boundary conditions using BEM++. The anode Dirichlet boundary condition can be specified in the configuration section (e.g. boundary) along with the electric field. If using predefined *tpcgeometry.msh*, this geometry is a simple box TPC with an anode, cathode, walls, pixels, and a focusing grid. The anode, cathode, and wall Dirichlet boundary conditions are determined by specifying anode voltage and electric field in the configuration section. To solve, run
```
$ pix boundary
```
To should have saved 2 arrays, dirichlet and neumann. To view contents of the store, *pixsim* offers the `dump` command:
```
$ pixsim dump all 
Results...
id: 1   name: domain_map  typename: geometry    data: 1   parent: None
id: 2   name: solution    typename: boundary    data: 2   parent: None
Arrays...
id: 1   name: domain_map  typename: tuples      shape: (36, 3)   
id: 2   name: dirichlet   typename: scalar      shape: (20080,)  
id: 3   name: neumann     typename: scalar      shape: (20080,)
```
Run with the `--help` flag to see other ways of dumping the contents of the store. 

## Raster
Next stage will evaluate solution on a volume of points determined by the *linspace* variable in the corresponding configuration section.
```
$ pix raster
```
This will compute the potential and gradient and save other relevent information. An option is provided to pass a specific boundary result ID or name using the `-b` flag. If none is passed, the algorithm will look for the most recent boundary result. As an example,
```
$ pix raster -b 2
```
or 
```
$ pix raster -b my_boundary_result_name
```

## Velocity 
To compute the drift velocity on the raster, run
```
$ pix velocity
```

## Stepping 
The stepping algorithm will use the results from the velocity stage and step a particle through the velocity field until termination. Stepping will terminate when the step sizes become too small or when the path intersects a electrode (e.g. pixel). The stepping algorithm needs vertices to step from. There are two ways to initialize the stepping:
1. Generate vertices on a rectangle patch given by variables in configuration section. 
2. Read vertices from a text file.

The stepping algorithm expects the following format from the text file:
```
vtx1 x y z 
vtx2 x y z 
vtx3 x y z 
etc...
```
The text file can be specified in the configuration file or passed at runtime as
```
$ pix -p step:stepfile=/path/to/step/file.txt step 
```
In general, one can override variables in the configuration file with the syntax
```
$ pix -p <section_name>:<variable_name>=<new_value> <other_commands>
```
Running stepping stage with a stepfile containing one vertex named *genvtx1*
```
$ pix -p step:stepfile=gen_step_vtx/pixel11_genvtx.txt step
```
will save 2 arrays:
1. List of vertices (named vtxs)
2. List of steps from each vertex (named genvtx1, [ x y z vx vy vz t])

## Weighting field
To compute the weighting field used in computing the instantaneous current, we will solve again for the Neumann boundary conditions with different Dirichlet boundary conditions. The weighting field of a particular electrode is defined as the field that is produced when holding that electrode at unit potential and all others at 0 V.

The only difference between this calculation and the regular boundary calculation is that we must specify *weighting* as the method and the domain ID of the particular electrode in the configuration section. For a section name *weight*, we can run
```
$ pix boundary -c weight -n my_first_weight
```
where the section name has been passed and name has been explicity given to the result using the `-n` flag. 

Now we must rerun the raster algorithm to calculate the field on our volume. 
```
$ pix raster -b my_first_weight -n my_first_weight_raster
```

## Instantaneous current
With the weighting field now in hand, we can calculate the current waveforms using the Shockley-Ramo relation. The ingredients are:
- Steps (position and drift velocity)
- Weighting field

For a step result ID of 5 and weight result name of *my_first_weight_raster*, we can run the current calculation:
```
$ pix current -s 5 -r my_first_weight_raster
```
This should save a waveform for each path from the step result. The waveform is a list of pairs [time, amplitude].

## Utilities
*pixsim* also provides functions for plotting, exporting, and renaming and deleting results. 

### Plotting
One can plot the potential from a raster result by running the command
```
$ pix plot -r <raster_result_id> -q potential
```
CAUTION: depending on how fine the parameters are in the configuration section and in the code itself, this could take a while.

### Exporting 
To export results to *.vtk* format for viewing with [ParaView](https://www.paraview.org/), run
```
$ pix export -s <result_type> -o <filename_prefix>
```

### Editing the store
Rename
```
$ pix rename --help
```
Remove
```
$ pix rm --help
```

# Generation 
To help automate the stages in the simulation, *pixsim* provides multiple subcommands under the *gen* command for running identical simulation for closely related input. As an example, one may like to simulate the current reponse for many electrodes, taking an average at the end. This requires simulations for identical relative drift vertices for multiple electrodes using their corresponding weighting field. 

## Weighting field
Run the weighting field stage for multiple domains through the command
```
$ pix gen weight
```
which will save a weighting field result for (configurable) each domain in the domain map. 

## Raster
Run the raster stage for these domains
```
$ pix gen raster -b <first_weight_res_id>:<last_weight_res_id>
```

## Vertices
Generate similar vertices using a template file containing the electrode IDs and the relative position to the electrode. This will create a new directory containing stepfiles for identical paths around each of the specific electrodes. 
```
$ pix gen vtx -s step_template.txt -d genvtxs
```

## Stepping
Run the stepping stage for each of the stepfiles created above.
```
$ pix gen step -d genvtxs
```

## Exporting
Export a range of results
```
$ pix gen export -s <result_type> -o <name_prefix> -r <first_result_id>:<last_result_id>
```
where result_type is *boundary*, *raster*, or *step*.
![alt text](https://raw.githubusercontent.com/hcsullivan12/picture_store/master/steps1.png "Steps")
