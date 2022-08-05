// Target element size of outer boundary
lc2 = (0.41 / 10) / 0.41;

// Bottom left coords
xL = 0;
yL = 0;

// Top right coords
xR = 2.2 / 0.41;
yR = 0.41 / 0.41;

p5 = newp; Point(p5) = {0.2 / 0.41, 0.2 / 0.41, 0, lc2};
p6 = newp; Point(p6) = {0.25 / 0.41, 0.2 / 0.41, 0, lc2};
p7 = newp; Point(p7) = {0.2 / 0.41, 0.15 / 0.41, 0, lc2};
p8 = newp; Point(p8) = {0.15 / 0.41, 0.2 / 0.41, 0, lc2};
p9 = newp; Point(p9) = {0.2 / 0.41, 0.25 / 0.41, 0, lc2};

Circle(1) = {p6, p5, p7};
Circle(2) = {p7, p5, p8};
Circle(3) = {p8, p5, p9};
Circle(4) = {p9, p5, p6};

Line Loop(1) = {2, 3, 4, 1};

// Outline points
p1 = newp; Point(p1) = {xL, yL, 0, lc2};
p2 = newp; Point(p2) = {xL, yR, 0, lc2};
p3 = newp; Point(p3) = {xR, yR, 0, lc2};
p4 = newp; Point(p4) = {xR, yL, 0, lc2};

// Outline lines
l1 = newl; Line(l1) = {p1, p2};
l2 = newl; Line(l2) = {p2, p3};
l3 = newl; Line(l3) = {p3, p4};
l4 = newl; Line(l4) = {p4, p1};

Curve Loop(2) = {l1, l2, l3, l4};

// Surface of mesh
Plane Surface(1) = {2, 1};
