// Target element size
lc = 0.2;

// Bottom left coords
xL = 0;
yL = 0;
zL = 0;

// Top right coords
xR = 1;
yR = 1;
zR = 1;

// Points defining the cube
// Front face
Point(1) = {xL, yL, zL, lc};
Point(2) = {xR, yL, zL, lc};
Point(3) = {xR, yR, zL, lc};
Point(4) = {xL, yR, zL, lc};
// Back face
Point(5) = {xL, yL, zR, lc};
Point(6) = {xR, yL, zR, lc};
Point(7) = {xR, yR, zR, lc};
Point(8) = {xL, yR, zR, lc};

// Edges of the cube
// Front face
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
// Back face
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
// Left connecting edges
Line(9) = {1, 5};
Line(10) = {4, 8};
// Right connecting edges
Line(11) = {2, 6};
Line(12) = {3, 7};

// Surfaces of the cube
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};
Curve Loop(3) = {2, 12, -6, -11};
Plane Surface(3) = {3};
Curve Loop(4) = {4, 9, -8, -10};
Plane Surface(4) = {4};
Curve Loop(5) = {1, 11, -5, -9};
Plane Surface(5) = {5};
Curve Loop(6) = {3, 10, -7, -12};
Plane Surface(6) = {6};

// Volume of the cube
Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};
Physical Volume("fluid") = {1};
