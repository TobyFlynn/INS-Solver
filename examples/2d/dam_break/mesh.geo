// Target element size
top_mesh_size = 0.04;
bottom_mesh_size = 0.006;

// Bottom left coords
xL = 0.0;
yL = 0.0;

// Top right coords
xR = 2.0;
yR = 2.0;

// Points defining the square
Point(1) = {xL, yL, 0.0, bottom_mesh_size};
Point(2) = {xR, yL, 0.0, bottom_mesh_size};
Point(3) = {xR, yR, 0.0, top_mesh_size};
Point(4) = {xL, yR, 0.0, top_mesh_size};

// Edges of the square
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Surface of square
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Physical Surface("fluid") = {1};
