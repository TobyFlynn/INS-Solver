// Target element size of outer boundary
grid_top = 0.025;
grid_bottom = 0.01;

// Bottom left coords
xL = 0;
yL = 0;

// Top right coords
xR = 2.0;
yR = 1.0;


// Outline points
p1 = newp; Point(p1) = {xL, yL, 0, grid_bottom};
p2 = newp; Point(p2) = {xL, yR, 0, grid_top};
p3 = newp; Point(p3) = {xR, yR, 0, grid_top};
p4 = newp; Point(p4) = {xR, yL, 0, grid_bottom};

// Outline lines
l1 = newl; Line(l1) = {p1, p2};
l2 = newl; Line(l2) = {p2, p3};
l3 = newl; Line(l3) = {p3, p4};
l4 = newl; Line(l4) = {p4, p1};

Curve Loop(1) = {l1, l2, l3, l4};

// Surface of mesh
Plane Surface(1) = {1};
