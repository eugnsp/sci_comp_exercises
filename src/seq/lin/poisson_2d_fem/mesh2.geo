Point(1) = {0, 0, 0, .75};
Point(2) = {3, 0, 0, .75};
Point(3) = {3, 3, 0, .75};
Point(4) = {0, 3, 0, .75};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = 1;
