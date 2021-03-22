// Target element size of naca0012
lc = 0.01;

// Number of points describing the upper surface
num = 100;

deltaX = (Pi / 2.0) / num;
xRef = 0.0;

// Surface points
// Upper
For i In {1:num-1}
  x = Sin(xRef)^2;
  y = 0.594689181 * (0.298222773*Sqrt(x) - 0.127125232*x - 0.357907906*(x^2) + 0.291984971*(x^3) - 0.105174696*(x^4));
  Point(i) = {x, y, 0, lc};
  xRef = xRef + deltaX;
EndFor
x = Sin(xRef)^2;
Point(num) = {x, 0, 0, lc};

// Lower
xRef = xRef - deltaX;
For i In {num+1:2*num - 2}
  x = Sin(xRef)^2;
  y = 0.594689181 * (0.298222773*Sqrt(x) - 0.127125232*x - 0.357907906*(x^2) + 0.291984971*(x^3) - 0.105174696*(x^4));
  Point(i) = {x, -y, 0, lc};
  xRef = xRef - deltaX;
EndFor

// Lines
For i In {1:2*num-3}
  Line(i) = {i, i+1};
EndFor

Line(2*num-2) = {2*num-2, 1};

// Curve loop of airfoil
Curve Loop(1) = {1:2*num-2};

// Surface of airfoil
//Plane Surface(1) = {1};

// Target element size of outer boundary
lc2 = 1;

// Bottom left coords
xL = -5;
yL = -2.5;

// Top right coords
xR = 5;
yR = 2.5;

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
