# multi-contour-extrusion
blender 2.79 add-on that is capable to create extruded geometries following path and going exactly through given number of contours. input can be both mesh and curve or even mixed

### Usage
create a set of contours with bezier curves or but direct editing of mesh. do not try to set new points in any specific order add-on will reorder them using edge adjacency (example: add a simple plane, subdivide individual edges, distort the cotour as you like feel free to remove face, but keep the edges, make sure your does not produce non-manifold geometry)

create a path also using bezier spline or mesh. cyclic pathes are not supported yet. it should work, but you have to stirch last segment manually.

run add-on from 3D view->Object->Multi-Contour Extrusion

play with settings. if your contours have different number of vertices you will get error message and no output. just tick "experimental" checkbox in left tab (Multi-Contour extrusion settings)

### Usage
###### Basic usage:

