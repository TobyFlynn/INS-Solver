mesh_size = 0.4;
inlet_size = 0.01;

pipe_radius = 1.0;
pipe_length = 5.0;
inlet_radius = 0.05;

Point(1) = {0.0, 0.0, 0.0, mesh_size};
Point(2) = {0.0, pipe_radius, 0.0, mesh_size};
Point(3) = {0.0, 0.0, pipe_radius, mesh_size};
Point(4) = {0.0, -pipe_radius, 0.0, mesh_size};
Point(5) = {0.0, 0.0, -pipe_radius, mesh_size};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
Curve Loop(1) = {1, 2, 3, 4};

Point(11) = {0.0, inlet_radius, 0.0, inlet_size};
Point(12) = {0.0, 0.0, inlet_radius, inlet_size};
Point(13) = {0.0, -inlet_radius, 0.0, inlet_size};
Point(14) = {0.0, 0.0, -inlet_radius, inlet_size};

Circle(13) = {11, 1, 12};
Circle(14) = {12, 1, 13};
Circle(15) = {13, 1, 14};
Circle(16) = {14, 1, 11};
Curve Loop(7) = {13, 14, 15, 16};

Plane Surface(1) = {1};
Curve{13,14,15,16} In Surface{1};

Point(6) = {pipe_length, 0.0, 0.0, mesh_size};
Point(7) = {pipe_length, pipe_radius, 0.0, mesh_size};
Point(8) = {pipe_length, 0.0, pipe_radius, mesh_size};
Point(9) = {pipe_length, -pipe_radius, 0.0, mesh_size};
Point(10) = {pipe_length, 0.0, -pipe_radius, mesh_size};

Circle(5) = {7, 6, 8};
Circle(6) = {8, 6, 9};
Circle(7) = {9, 6, 10};
Circle(8) = {10, 6, 7};

Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};

Line(9)  = {2, 7};
Line(10) = {3, 8};
Line(11) = {4, 9};
Line(12) = {5, 10};

Curve Loop(3) = {1, 10, -5, -9};
Curve Loop(4) = {2, 11, -6, -10};
Curve Loop(5) = {3, 12, -7, -11};
Curve Loop(6) = {4, 9, -8, -12};
Surface(3) = {3};
Surface(4) = {4};
Surface(5) = {5};
Surface(6) = {6};

Surface Loop(1) = {4, 1, 3, 2, 5, 6};
Volume(1) = {1};
