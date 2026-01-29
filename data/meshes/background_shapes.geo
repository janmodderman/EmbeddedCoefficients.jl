SetFactory("OpenCASCADE");
lc = 1.0;

//Point(1) = {0.0, -2.0, 0, 0.04*lc};
//Point(2) = {4.0, -2.0, 0, 0.07*lc};
//Point(3) = {4.0, 0.0, 0, 0.1*lc};
//Point(4) = {0.0, 0.0, 0, 0.001*lc};


Point(1) = {0.0, -8.0, 0, 0.01*lc};
Point(2) = {12.0, -8.0, 0, 0.5*lc};
Point(3) = {12.0, 0.0, 0, 0.01*lc};
Point(4) = {0.0, 0.0, 0, 0.003*lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(5) = {3,4,1,2};
Plane Surface(6) = {5};

//lc2 = 0.005;
lc2 = 0.003;

//Point(6) = {0.1, -0.0001, 0, lc2};
//Point{6} In Surface{6};
Point(5) = {0.1, -0.1, 0, lc2};
Point{5} In Surface{6};
//Point(7) = {0.001, -0.1, 0, lc2};
//Point{7} In Surface{6};
//Point(8) = {0.1, -0.01, 0, lc2};
//Point{8} In Surface{6};
//Point(9) = {0.01, -0.1, 0, lc2};
//Point{9} In Surface{6};

//+
//Transfinite Curve {3} = 400 Using Progression 0.97; // surface
//+
//Transfinite Curve {2} = 50 Using Progression 1; // outlet
//+
//Transfinite Curve {1} = 50 Using Progression 1; // seabed
//+
//Transfinite Curve {4} = 150 Using Progression 1.05; // inlet
//+
Physical Curve("outlet", 7) = {2};
//+
Physical Curve("inlet", 8) = {4};
//+
Physical Curve("surface", 9) = {3};
//+
Physical Curve("seabed", 10) = {1};
