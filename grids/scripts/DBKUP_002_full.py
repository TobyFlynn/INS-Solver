import csv

x = []
y = []
x_top = []
y_top = []
x_bottom = []
y_bottom = []

with open('DBKUP_002.csv') as csvfile:
  csvreader = csv.reader(csvfile, delimiter=',')
  i = 0
  for row in csvreader:
    x.append(row[0])
    y.append(row[1])
    i = i + 1

with open('DBKUP_002_full.geo', 'w') as geofile:
  # Actual airfoil
  line = "airfoil_tol = 1;\n"
  geofile.write(line)
  line = "factor = 0.47 / 4.9e-4;\n"
  geofile.write(line)
  for i in range(0, len(x)):
    line = "Point(%i) = { %8.8f * factor, %8.8f * factor, 0, airfoil_tol};\n" % (i,float(x[i]),float(y[i]))
    geofile.write(line)
  line = "Spline(0) = { 0 : %i, 0 };\nCurve Loop(1) = {0};\n" % (len(x) - 1)
  geofile.write(line)

  # Larger outline for making mesh fine near airfoil plus outline of droplet path
  line = "factor_outline = (0.47 / 4.9e-4) * 1.1;\n"
  geofile.write(line)
  line = "x_offset_outline = (factor_outline - factor) / 2.0;\n"
  geofile.write(line)
  for i in range(0, 143 - len(x)):
    p_id = i + len(x)
    line = "Point(%i) = { %8.8f * factor_outline - x_offset_outline, %8.8f * factor_outline, 0, airfoil_tol};\n" % (p_id,float(x[i]),float(y[i]))
    geofile.write(line)
  
  p_id = 143
  line = "droplet_release = -0.17 * factor;\n"
  geofile.write(line)
  line = "droplet_path_height = 10;\n"
  geofile.write(line)
  line = "Point(%i) = { droplet_release, -droplet_path_height, 0, airfoil_tol };\n" % (p_id)
  geofile.write(line)
  line = "Point(%i) = { droplet_release, droplet_path_height, 0, airfoil_tol };\n" % (p_id + 1)
  geofile.write(line)
  
  for i in range(145 - len(x) + 1, len(x)):
    p_id = i + len(x) - 1
    line = "Point(%i) = { %8.8f * factor_outline - x_offset_outline, %8.8f * factor_outline, 0, airfoil_tol};\n" % (p_id,float(x[i]),float(y[i]))
    geofile.write(line)
  line = "Spline(1) = { %i : %i, %i };\n" % (len(x), 2 * len(x) - 2, len(x))
  geofile.write(line)

  boundaryPts = \
"""boundary_tol = 100;
xL = -2 * factor;
yL = -1.5 * factor;
xR = 6 * factor;
yR = 1.5 * factor;
// Outline points
p1 = newp; Point(p1) = {xL, yL, 0, boundary_tol};
p2 = newp; Point(p2) = {xL, yR, 0, boundary_tol};
p3 = newp; Point(p3) = {xR, yR, 0, boundary_tol};
p4 = newp; Point(p4) = {xR, yL, 0, boundary_tol};

// Outline lines
l1 = newl; Line(l1) = {p1, p2};
l2 = newl; Line(l2) = {p2, p3};
l3 = newl; Line(l3) = {p3, p4};
l4 = newl; Line(l4) = {p4, p1};

Curve Loop(2) = {l1, l2, l3, l4};
Plane Surface(1) = {2, 1};
Curve{1} In Surface{1};\n"""
  geofile.write(boundaryPts)


