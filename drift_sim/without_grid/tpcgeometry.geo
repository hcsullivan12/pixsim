// This code was created by pygmsh v5.0.2.
SetFactory("OpenCASCADE");
vol0 = newv;
Box(vol0) = {0, -1.75, -1.75, 3.0, 3.5, 3.5};
vol1 = newv;
Sphere(vol1) = {0.25, 1.2, -1.2, 0.04};
vol2 = newv;
Sphere(vol2) = {0.25, 1.2, -0.8, 0.04};
vol3 = newv;
Sphere(vol3) = {0.25, 1.2, -0.4, 0.04};
vol4 = newv;
Sphere(vol4) = {0.25, 1.2, 0.0, 0.04};
vol5 = newv;
Sphere(vol5) = {0.25, 1.2, 0.4, 0.04};
vol6 = newv;
Sphere(vol6) = {0.25, 1.2, 0.8, 0.04};
vol7 = newv;
Sphere(vol7) = {0.25, 1.2, 1.2, 0.04};
vol8 = newv;
Sphere(vol8) = {0.25, 0.8, -1.2, 0.04};
vol9 = newv;
Sphere(vol9) = {0.25, 0.8, -0.8, 0.04};
vol10 = newv;
Sphere(vol10) = {0.25, 0.8, -0.4, 0.04};
vol11 = newv;
Sphere(vol11) = {0.25, 0.8, 0.0, 0.04};
vol12 = newv;
Sphere(vol12) = {0.25, 0.8, 0.4, 0.04};
vol13 = newv;
Sphere(vol13) = {0.25, 0.8, 0.8, 0.04};
vol14 = newv;
Sphere(vol14) = {0.25, 0.8, 1.2, 0.04};
vol15 = newv;
Sphere(vol15) = {0.25, 0.4, -1.2, 0.04};
vol16 = newv;
Sphere(vol16) = {0.25, 0.4, -0.8, 0.04};
vol17 = newv;
Sphere(vol17) = {0.25, 0.4, -0.4, 0.04};
vol18 = newv;
Sphere(vol18) = {0.25, 0.4, 0.0, 0.04};
vol19 = newv;
Sphere(vol19) = {0.25, 0.4, 0.4, 0.04};
vol20 = newv;
Sphere(vol20) = {0.25, 0.4, 0.8, 0.04};
vol21 = newv;
Sphere(vol21) = {0.25, 0.4, 1.2, 0.04};
vol22 = newv;
Sphere(vol22) = {0.25, 0.0, -1.2, 0.04};
vol23 = newv;
Sphere(vol23) = {0.25, 0.0, -0.8, 0.04};
vol24 = newv;
Sphere(vol24) = {0.25, 0.0, -0.4, 0.04};
vol25 = newv;
Sphere(vol25) = {0.25, 0.0, 0.0, 0.04};
vol26 = newv;
Sphere(vol26) = {0.25, 0.0, 0.4, 0.04};
vol27 = newv;
Sphere(vol27) = {0.25, 0.0, 0.8, 0.04};
vol28 = newv;
Sphere(vol28) = {0.25, 0.0, 1.2, 0.04};
vol29 = newv;
Sphere(vol29) = {0.25, -0.4, -1.2, 0.04};
vol30 = newv;
Sphere(vol30) = {0.25, -0.4, -0.8, 0.04};
vol31 = newv;
Sphere(vol31) = {0.25, -0.4, -0.4, 0.04};
vol32 = newv;
Sphere(vol32) = {0.25, -0.4, 0.0, 0.04};
vol33 = newv;
Sphere(vol33) = {0.25, -0.4, 0.4, 0.04};
vol34 = newv;
Sphere(vol34) = {0.25, -0.4, 0.8, 0.04};
vol35 = newv;
Sphere(vol35) = {0.25, -0.4, 1.2, 0.04};
vol36 = newv;
Sphere(vol36) = {0.25, -0.8, -1.2, 0.04};
vol37 = newv;
Sphere(vol37) = {0.25, -0.8, -0.8, 0.04};
vol38 = newv;
Sphere(vol38) = {0.25, -0.8, -0.4, 0.04};
vol39 = newv;
Sphere(vol39) = {0.25, -0.8, 0.0, 0.04};
vol40 = newv;
Sphere(vol40) = {0.25, -0.8, 0.4, 0.04};
vol41 = newv;
Sphere(vol41) = {0.25, -0.8, 0.8, 0.04};
vol42 = newv;
Sphere(vol42) = {0.25, -0.8, 1.2, 0.04};
vol43 = newv;
Sphere(vol43) = {0.25, -1.2, -1.2, 0.04};
vol44 = newv;
Sphere(vol44) = {0.25, -1.2, -0.8, 0.04};
vol45 = newv;
Sphere(vol45) = {0.25, -1.2, -0.4, 0.04};
vol46 = newv;
Sphere(vol46) = {0.25, -1.2, 0.0, 0.04};
vol47 = newv;
Sphere(vol47) = {0.25, -1.2, 0.4, 0.04};
vol48 = newv;
Sphere(vol48) = {0.25, -1.2, 0.8, 0.04};
vol49 = newv;
Sphere(vol49) = {0.25, -1.2, 1.2, 0.04};
//+
Physical Surface("cathode") = {2};
//+
Physical Surface("anode") = {1};
//+
Physical Surface("pixel8") = {14};
//+
Physical Surface("pixel9") = {15};
//+
Physical Surface("pixel6") = {12};
//+
Physical Surface("pixel7") = {13};
//+
Physical Surface("pixel4") = {10};
//+
Physical Surface("pixel5") = {11};
//+
Physical Surface("pixel2") = {8};
//+
Physical Surface("pixel3") = {9};
//+
Physical Surface("pixel1") = {7};
//+
Physical Surface("pixel36") = {42};
//+
Physical Surface("pixel37") = {43};
//+
Physical Surface("pixel34") = {40};
//+
Physical Surface("pixel35") = {41};
//+
Physical Surface("pixel32") = {38};
//+
Physical Surface("pixel33") = {39};
//+
Physical Surface("pixel30") = {36};
//+
Physical Surface("pixel31") = {37};
//+
Physical Surface("pixel38") = {44};
//+
Physical Surface("pixel39") = {45};
//+
Physical Surface("pixel18") = {24};
//+
Physical Surface("pixel19") = {25};
//+
Physical Surface("pixel14") = {20};
//+
Physical Surface("pixel15") = {21};
//+
Physical Surface("pixel16") = {22};
//+
Physical Surface("pixel17") = {23};
//+
Physical Surface("pixel10") = {16};
//+
Physical Surface("pixel11") = {17};
//+
Physical Surface("pixel12") = {18};
//+
Physical Surface("pixel13") = {19};
//+
Physical Surface("walls") = {3, 4, 5, 6};
//+
Physical Surface("pixel43") = {49};
//+
Physical Surface("pixel42") = {48};
//+
Physical Surface("pixel41") = {47};
//+
Physical Surface("pixel40") = {46};
//+
Physical Surface("pixel47") = {53};
//+
Physical Surface("pixel46") = {52};
//+
Physical Surface("pixel45") = {51};
//+
Physical Surface("pixel44") = {50};
//+
Physical Surface("pixel49") = {55};
//+
Physical Surface("pixel48") = {54};
//+
Physical Surface("pixel29") = {35};
//+
Physical Surface("pixel28") = {34};
//+
Physical Surface("pixel21") = {27};
//+
Physical Surface("pixel20") = {26};
//+
Physical Surface("pixel23") = {29};
//+
Physical Surface("pixel22") = {28};
//+
Physical Surface("pixel25") = {31};
//+
Physical Surface("pixel24") = {30};
//+
Physical Surface("pixel27") = {33};
//+
Physical Surface("pixel26") = {32};//+
Transfinite Line {12, 6, 10, 2, 1, 4, 3, 11, 8, 7, 5, 9} = 15 Using Progression 1;
//+
Transfinite Surface {6} = {3, 7, 5, 1};
//+
Transfinite Surface {4} = {4, 8, 7, 3};
//+
Transfinite Surface {5} = {4, 2, 6, 8};
//+
Transfinite Surface {3} = {1, 5, 6, 2};
//+
Transfinite Surface {1} = {4, 3, 1, 2};
//+
Transfinite Surface {2} = {8, 6, 5, 7};
