// Target element size
lc = 0.22;

// Bottom left coords
xL = 0.0;
yL = 0.0;
zL = 0.0;

// Top right coords
xR = 2 * Pi;
yR = 2 * Pi;
zR = 2 * Pi;

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

// Periodic stuff
Periodic Surface{1} = {2} Translate {0,0,-zR};
Periodic Surface{4} = {3} Translate {-xR,0,0};
Periodic Surface{5} = {6} Translate {0,-yR,0};

// Volume of the cube
Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};
Physical Surface("periodic_0_l") = {1};
Physical Surface("periodic_0_r") = {2};
Physical Surface("periodic_1_l") = {4};
Physical Surface("periodic_1_r") = {3};
Physical Surface("periodic_2_l") = {5};
Physical Surface("periodic_2_r") = {6};
Physical Volume("fluid") = {1};
