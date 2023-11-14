// Target element size
lc = 1.0 / 3.0;

// Bottom left coords
xL = -10.0;
yL = -10.0;

// Top right coords
xR = 10.0;
yR = 10.0;

// Points defining the square
Point(1) = {xL, yL, 0.0, lc};
Point(2) = {xR, yL, 0.0, lc};
Point(3) = {xR, yR, 0.0, lc};
Point(4) = {xL, yR, 0.0, lc};

// Edges of the square
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Periodic stuff
Periodic Curve{1} = {3} Translate {0,-2*yR,0};
Periodic Curve{2} = {4} Translate {-2*xR,0,0};

// Surface of square
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Physical Surface("fluid") = {1};
