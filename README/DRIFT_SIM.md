# Integrating pixsim with a drift simulation
Currently, pixsim does not offer a drift simulation. Instead, a tool like [pixgen](https://github.com/hcsullivan/pixgen) can be used to derive the charge distribution near a readout plane taking into account attenuation and diffusion effects. [pixgen](https://github.com/hcsullivan/pixgen) is a modified version of [larsim](https://cdcvs.fnal.gov/redmine/projects/larsim) that uses [Geant4](https://geant4.web.cern.ch/) to simulate various particles inside a liquid argon time projection chamber. One of the standard routines is an ionization charge drift simulation. The items that are needed for pixsim are the (x,y,z) coordinates of the drifted charge and the pixel coordinates. For the remainder of this section, it is assumed that a [ROOT](https://root.cern.ch/) ntuple contains this information.

Note: It may be useful to walk through the instructions here [README](https://github.com/hcsullivan/pixgen/README.md) first as some details have been left out.

# Example 
This section will walk through a simulation integrating a drift simulation. Let's use the alias `pix` defined as
```
$ alias pix='pixsim -c boxtpc.cfg -s test.db -m tpcgeometry.msh'
```

## Initial steps
Generate the geometry
```
$ pix gen geo -c <geometry_section_name> -o <name_of_geo_file>
```
Generate the domain map
```
$ pix gen dmap
```
Run BEM++
```
$ pix boundary
```
Evaluate
```
$ pix raster
```
Compute velocity
```
$ pix velocity
```

## Deriving the responses
The path chosen in pixsim is to derive the response using information from a single pixel. To do this, we select the center pixel as the pixel of interest which has a domain ID of 50. 

Calculate weighting field for pixel domain ID 50
```
$ pix boundary -c weighting
```
Generate the vertices to be used in the stepping algorithm. Due to symmetry, we only need to generate vertices in the half quadrant. 
```
$ pix gen vtx -s step_lookup.txt
```
Run the stepping algorithm for these vertices
```
$ pix gen step 
```
Calculate induced current
```
$ pix gen current 
```
Convolve field response with electronics response
```
$ pix gen response -w <current_result_id> -s <step_result_id>
```
Run the simulation 
```
$ pix sim
```

