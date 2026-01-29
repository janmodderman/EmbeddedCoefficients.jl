// Gmsh project created on Wed Jan 15 10:19:34 2025
// Gmsh project created on Tue Jan 14 09:27:40 2025
SetFactory("OpenCASCADE");
//+
Cylinder(1) = {0, 0, -200, 0, 0, 200, 300, 2*Pi};
//+
//Point(3) = {0, 0, 0, 0.5};
//+
//Point(4) = {0, 0, -125, 0.5};
//+
//Point(5) = {0, 0, -200, 1.0};
//Point{3} In Surface{2};
//Point{5} In Surface{3};
//Point{4} In Volume{1};

//Transfinite Curve {1} = 200 Using Progression 1.0;
//+
//Transfinite Curve {3} = 20 Using Progression 1.0;
//+
//Transfinite Curve {2} = 20 Using Progression 1.0;
//+
Physical Volume("water", 4) = {1};
//+
Physical Surface("surface", 5) = {2};
//+
Physical Surface("seabed", 6) = {3};
//+
Physical Surface("walls", 7) = {1};//+
Field[1] = Cylinder;
//+
Field[1].VIn = 0.4;
//+
Field[1].ZAxis = -22;
//+
Field[1].ZCenter = 0;
Field[1].XAxis = 0;
//+
Field[1].XCenter = 0;
Field[1].YAxis = 0;
//+
Field[1].YCenter = 0;
//+
Field[1].Radius = 50;
//+
Field[1].VOut = 100.0;
//+
Background Field = 1;
//+//+
Field[2] = Cylinder;
//+
Field[2].Radius = 1500;
//+
Field[2].VIn = 2.5;
Field[2].ZAxis = 1;
//+
Field[2].ZCenter = 0;
Field[2].XAxis = 0;
//+
Field[2].XCenter = 0;
Field[2].YAxis = 0;
//+
Field[2].YCenter = 0;
//+
Background Field = 2;
//+
Field[3] = Min;
//+
Field[3].FieldsList = {1, 2};
//+
Background Field = 3;
