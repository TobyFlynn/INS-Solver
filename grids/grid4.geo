// Target element size of outer boundary
coarse = 0.00652;
fine = 0.004;

// Bottom left coords
xL = 0;
yL = 0;

// Top right coords
xR = 5.5;
yR = 1.0;


p5 = newp; Point(p5) = {1.0, 0.5, 0, coarse};
p6 = newp; Point(p6) = {1.2, 0.5, 0, fine};
p7 = newp; Point(p7) = {1.0, 0.3, 0, fine};
p8 = newp; Point(p8) = {0.8, 0.5, 0, fine};
p9 = newp; Point(p9) = {1.0, 0.7, 0, fine};

Circle(1) = {p6, p5, p7};
Circle(2) = {p7, p5, p8};
Circle(3) = {p8, p5, p9};
Circle(4) = {p9, p5, p6};

p16 = newp; Point(p16) = {1.25, 0.5, 0, fine};
p17 = newp; Point(p17) = {1.0, 0.25, 0, fine};
p18 = newp; Point(p18) = {0.75, 0.5, 0, fine};
p19 = newp; Point(p19) = {1.0, 0.75, 0, fine};

Circle(11) = {p16, p5, p17};
Circle(12) = {p17, p5, p18};
Circle(13) = {p18, p5, p19};
Circle(14) = {p19, p5, p16};

p26 = newp; Point(p26) = {1.15, 0.5, 0, fine};
p27 = newp; Point(p27) = {1.0, 0.35, 0, fine};
p28 = newp; Point(p28) = {0.85, 0.5, 0, fine};
p29 = newp; Point(p29) = {1.0, 0.65, 0, fine};

Circle(21) = {p26, p5, p27};
Circle(22) = {p27, p5, p28};
Circle(23) = {p28, p5, p29};
Circle(24) = {p29, p5, p26};

// Outline points
p1 = newp; Point(p1) = {xL, yL, 0, coarse};
p2 = newp; Point(p2) = {xL, yR, 0, coarse};
p3 = newp; Point(p3) = {xR, yR, 0, coarse};
p4 = newp; Point(p4) = {xR, yL, 0, coarse};

// Outline lines
l1 = newl; Line(l1) = {p1, p2};
l2 = newl; Line(l2) = {p2, p3};
l3 = newl; Line(l3) = {p3, p4};
l4 = newl; Line(l4) = {p4, p1};

Curve Loop(1) = {l1, l2, l3, l4};

// Surface of mesh
Plane Surface(1) = {1};

Point{p5} In Surface{1};

Curve{1,2,3,4,11,12,13,14,21,22,23,24} In Surface{1};
