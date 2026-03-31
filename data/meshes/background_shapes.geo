SetFactory("OpenCASCADE");
lc = 0.5;       // reduce to increase number of elements
lc2 = 0.001;    // reduce to increase number of elements near point where domain is cut

// Domain dimensions
Dx = 3.0;       // domain length
Dy = 2.0;       // domain depth
Point(1) = {0.0, -Dy, 0, 0.01*lc};
Point(2) = {Dx, -Dy, 0, 0.5*lc};
Point(3) = {Dx, 0.0, 0, 0.01*lc};
Point(4) = {0.0, 0.0, 0, 0.005*lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(5) = {3,4,1,2};
Plane Surface(6) = {5};

Point(5) = {0.1, -0.1, 0, lc2};
Point{5} In Surface{6};

Physical Curve("walls", 7) = {2};
//+
Physical Curve("symmetry", 8) = {4};
//+
Physical Curve("surface", 9) = {3};
//+
Physical Curve("seabed", 10) = {1};