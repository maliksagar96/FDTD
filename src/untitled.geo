//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, 1, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {1, 0, 0, 1.0};
//+
Line(1) = {2, 1};
//+
Line(2) = {3, 2};
//+
Line(3) = {4, 3};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {2, 1, -4, 3};
//+
Plane Surface(1) = {1};
//+
Physical Curve("Inlet", 5) = {1};
//+
Physical Curve("Outlet", 6) = {3};
//+
Physical Curve("Top", 7) = {2};
//+
Physical Curve("Bottom", 8) = {4};
//+
Transfinite Surface {1} = {3, 4, 1, 2};
//+
Transfinite Curve {1, 3, 2, 4} = 20 Using Progression 1;
//+
Recombine Surface {1};
