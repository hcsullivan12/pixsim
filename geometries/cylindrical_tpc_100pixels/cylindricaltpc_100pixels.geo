// This code was created by pygmsh v5.0.2.
SetFactory("OpenCASCADE");
vol0 = newv;
Cylinder(vol0) = {-0.05, 0, 0, 0.1, 0, 0, 5.0};
vol1 = newv;
Cylinder(vol1) = {4.95, 0, 0, 0.1, 0, 0, 5.0};
vol2 = newv;
Cylinder(vol2) = {0.15, 0, 0, 4.85, 0, 0, 4.9};
vol3 = newv;
Cylinder(vol3) = {0.15, 0, 0, 4.85, 0, 0, 5.0};
bo1[] = BooleanDifference{ Volume{vol3}; Delete; } { Volume{vol2}; Delete;};
vol4 = newv;
Cylinder(vol4) = {0.25, 1.8, -1.8, 0.04, 0, 0, 0.04};
vol5 = newv;
Cylinder(vol5) = {0.25, 1.8, -1.4, 0.04, 0, 0, 0.04};
vol6 = newv;
Cylinder(vol6) = {0.25, 1.8, -1.0, 0.04, 0, 0, 0.04};
vol7 = newv;
Cylinder(vol7) = {0.25, 1.8, -0.6, 0.04, 0, 0, 0.04};
vol8 = newv;
Cylinder(vol8) = {0.25, 1.8, -0.2, 0.04, 0, 0, 0.04};
vol9 = newv;
Cylinder(vol9) = {0.25, 1.8, 0.2, 0.04, 0, 0, 0.04};
vol10 = newv;
Cylinder(vol10) = {0.25, 1.8, 0.6, 0.04, 0, 0, 0.04};
vol11 = newv;
Cylinder(vol11) = {0.25, 1.8, 1.0, 0.04, 0, 0, 0.04};
vol12 = newv;
Cylinder(vol12) = {0.25, 1.8, 1.4, 0.04, 0, 0, 0.04};
vol13 = newv;
Cylinder(vol13) = {0.25, 1.8, 1.8, 0.04, 0, 0, 0.04};
vol14 = newv;
Cylinder(vol14) = {0.25, 1.4, -1.8, 0.04, 0, 0, 0.04};
vol15 = newv;
Cylinder(vol15) = {0.25, 1.4, -1.4, 0.04, 0, 0, 0.04};
vol16 = newv;
Cylinder(vol16) = {0.25, 1.4, -1.0, 0.04, 0, 0, 0.04};
vol17 = newv;
Cylinder(vol17) = {0.25, 1.4, -0.6, 0.04, 0, 0, 0.04};
vol18 = newv;
Cylinder(vol18) = {0.25, 1.4, -0.2, 0.04, 0, 0, 0.04};
vol19 = newv;
Cylinder(vol19) = {0.25, 1.4, 0.2, 0.04, 0, 0, 0.04};
vol20 = newv;
Cylinder(vol20) = {0.25, 1.4, 0.6, 0.04, 0, 0, 0.04};
vol21 = newv;
Cylinder(vol21) = {0.25, 1.4, 1.0, 0.04, 0, 0, 0.04};
vol22 = newv;
Cylinder(vol22) = {0.25, 1.4, 1.4, 0.04, 0, 0, 0.04};
vol23 = newv;
Cylinder(vol23) = {0.25, 1.4, 1.8, 0.04, 0, 0, 0.04};
vol24 = newv;
Cylinder(vol24) = {0.25, 1.0, -1.8, 0.04, 0, 0, 0.04};
vol25 = newv;
Cylinder(vol25) = {0.25, 1.0, -1.4, 0.04, 0, 0, 0.04};
vol26 = newv;
Cylinder(vol26) = {0.25, 1.0, -1.0, 0.04, 0, 0, 0.04};
vol27 = newv;
Cylinder(vol27) = {0.25, 1.0, -0.6, 0.04, 0, 0, 0.04};
vol28 = newv;
Cylinder(vol28) = {0.25, 1.0, -0.2, 0.04, 0, 0, 0.04};
vol29 = newv;
Cylinder(vol29) = {0.25, 1.0, 0.2, 0.04, 0, 0, 0.04};
vol30 = newv;
Cylinder(vol30) = {0.25, 1.0, 0.6, 0.04, 0, 0, 0.04};
vol31 = newv;
Cylinder(vol31) = {0.25, 1.0, 1.0, 0.04, 0, 0, 0.04};
vol32 = newv;
Cylinder(vol32) = {0.25, 1.0, 1.4, 0.04, 0, 0, 0.04};
vol33 = newv;
Cylinder(vol33) = {0.25, 1.0, 1.8, 0.04, 0, 0, 0.04};
vol34 = newv;
Cylinder(vol34) = {0.25, 0.6, -1.8, 0.04, 0, 0, 0.04};
vol35 = newv;
Cylinder(vol35) = {0.25, 0.6, -1.4, 0.04, 0, 0, 0.04};
vol36 = newv;
Cylinder(vol36) = {0.25, 0.6, -1.0, 0.04, 0, 0, 0.04};
vol37 = newv;
Cylinder(vol37) = {0.25, 0.6, -0.6, 0.04, 0, 0, 0.04};
vol38 = newv;
Cylinder(vol38) = {0.25, 0.6, -0.2, 0.04, 0, 0, 0.04};
vol39 = newv;
Cylinder(vol39) = {0.25, 0.6, 0.2, 0.04, 0, 0, 0.04};
vol40 = newv;
Cylinder(vol40) = {0.25, 0.6, 0.6, 0.04, 0, 0, 0.04};
vol41 = newv;
Cylinder(vol41) = {0.25, 0.6, 1.0, 0.04, 0, 0, 0.04};
vol42 = newv;
Cylinder(vol42) = {0.25, 0.6, 1.4, 0.04, 0, 0, 0.04};
vol43 = newv;
Cylinder(vol43) = {0.25, 0.6, 1.8, 0.04, 0, 0, 0.04};
vol44 = newv;
Cylinder(vol44) = {0.25, 0.2, -1.8, 0.04, 0, 0, 0.04};
vol45 = newv;
Cylinder(vol45) = {0.25, 0.2, -1.4, 0.04, 0, 0, 0.04};
vol46 = newv;
Cylinder(vol46) = {0.25, 0.2, -1.0, 0.04, 0, 0, 0.04};
vol47 = newv;
Cylinder(vol47) = {0.25, 0.2, -0.6, 0.04, 0, 0, 0.04};
vol48 = newv;
Cylinder(vol48) = {0.25, 0.2, -0.2, 0.04, 0, 0, 0.04};
vol49 = newv;
Cylinder(vol49) = {0.25, 0.2, 0.2, 0.04, 0, 0, 0.04};
vol50 = newv;
Cylinder(vol50) = {0.25, 0.2, 0.6, 0.04, 0, 0, 0.04};
vol51 = newv;
Cylinder(vol51) = {0.25, 0.2, 1.0, 0.04, 0, 0, 0.04};
vol52 = newv;
Cylinder(vol52) = {0.25, 0.2, 1.4, 0.04, 0, 0, 0.04};
vol53 = newv;
Cylinder(vol53) = {0.25, 0.2, 1.8, 0.04, 0, 0, 0.04};
vol54 = newv;
Cylinder(vol54) = {0.25, -0.2, -1.8, 0.04, 0, 0, 0.04};
vol55 = newv;
Cylinder(vol55) = {0.25, -0.2, -1.4, 0.04, 0, 0, 0.04};
vol56 = newv;
Cylinder(vol56) = {0.25, -0.2, -1.0, 0.04, 0, 0, 0.04};
vol57 = newv;
Cylinder(vol57) = {0.25, -0.2, -0.6, 0.04, 0, 0, 0.04};
vol58 = newv;
Cylinder(vol58) = {0.25, -0.2, -0.2, 0.04, 0, 0, 0.04};
vol59 = newv;
Cylinder(vol59) = {0.25, -0.2, 0.2, 0.04, 0, 0, 0.04};
vol60 = newv;
Cylinder(vol60) = {0.25, -0.2, 0.6, 0.04, 0, 0, 0.04};
vol61 = newv;
Cylinder(vol61) = {0.25, -0.2, 1.0, 0.04, 0, 0, 0.04};
vol62 = newv;
Cylinder(vol62) = {0.25, -0.2, 1.4, 0.04, 0, 0, 0.04};
vol63 = newv;
Cylinder(vol63) = {0.25, -0.2, 1.8, 0.04, 0, 0, 0.04};
vol64 = newv;
Cylinder(vol64) = {0.25, -0.6, -1.8, 0.04, 0, 0, 0.04};
vol65 = newv;
Cylinder(vol65) = {0.25, -0.6, -1.4, 0.04, 0, 0, 0.04};
vol66 = newv;
Cylinder(vol66) = {0.25, -0.6, -1.0, 0.04, 0, 0, 0.04};
vol67 = newv;
Cylinder(vol67) = {0.25, -0.6, -0.6, 0.04, 0, 0, 0.04};
vol68 = newv;
Cylinder(vol68) = {0.25, -0.6, -0.2, 0.04, 0, 0, 0.04};
vol69 = newv;
Cylinder(vol69) = {0.25, -0.6, 0.2, 0.04, 0, 0, 0.04};
vol70 = newv;
Cylinder(vol70) = {0.25, -0.6, 0.6, 0.04, 0, 0, 0.04};
vol71 = newv;
Cylinder(vol71) = {0.25, -0.6, 1.0, 0.04, 0, 0, 0.04};
vol72 = newv;
Cylinder(vol72) = {0.25, -0.6, 1.4, 0.04, 0, 0, 0.04};
vol73 = newv;
Cylinder(vol73) = {0.25, -0.6, 1.8, 0.04, 0, 0, 0.04};
vol74 = newv;
Cylinder(vol74) = {0.25, -1.0, -1.8, 0.04, 0, 0, 0.04};
vol75 = newv;
Cylinder(vol75) = {0.25, -1.0, -1.4, 0.04, 0, 0, 0.04};
vol76 = newv;
Cylinder(vol76) = {0.25, -1.0, -1.0, 0.04, 0, 0, 0.04};
vol77 = newv;
Cylinder(vol77) = {0.25, -1.0, -0.6, 0.04, 0, 0, 0.04};
vol78 = newv;
Cylinder(vol78) = {0.25, -1.0, -0.2, 0.04, 0, 0, 0.04};
vol79 = newv;
Cylinder(vol79) = {0.25, -1.0, 0.2, 0.04, 0, 0, 0.04};
vol80 = newv;
Cylinder(vol80) = {0.25, -1.0, 0.6, 0.04, 0, 0, 0.04};
vol81 = newv;
Cylinder(vol81) = {0.25, -1.0, 1.0, 0.04, 0, 0, 0.04};
vol82 = newv;
Cylinder(vol82) = {0.25, -1.0, 1.4, 0.04, 0, 0, 0.04};
vol83 = newv;
Cylinder(vol83) = {0.25, -1.0, 1.8, 0.04, 0, 0, 0.04};
vol84 = newv;
Cylinder(vol84) = {0.25, -1.4, -1.8, 0.04, 0, 0, 0.04};
vol85 = newv;
Cylinder(vol85) = {0.25, -1.4, -1.4, 0.04, 0, 0, 0.04};
vol86 = newv;
Cylinder(vol86) = {0.25, -1.4, -1.0, 0.04, 0, 0, 0.04};
vol87 = newv;
Cylinder(vol87) = {0.25, -1.4, -0.6, 0.04, 0, 0, 0.04};
vol88 = newv;
Cylinder(vol88) = {0.25, -1.4, -0.2, 0.04, 0, 0, 0.04};
vol89 = newv;
Cylinder(vol89) = {0.25, -1.4, 0.2, 0.04, 0, 0, 0.04};
vol90 = newv;
Cylinder(vol90) = {0.25, -1.4, 0.6, 0.04, 0, 0, 0.04};
vol91 = newv;
Cylinder(vol91) = {0.25, -1.4, 1.0, 0.04, 0, 0, 0.04};
vol92 = newv;
Cylinder(vol92) = {0.25, -1.4, 1.4, 0.04, 0, 0, 0.04};
vol93 = newv;
Cylinder(vol93) = {0.25, -1.4, 1.8, 0.04, 0, 0, 0.04};
vol94 = newv;
Cylinder(vol94) = {0.25, -1.8, -1.8, 0.04, 0, 0, 0.04};
vol95 = newv;
Cylinder(vol95) = {0.25, -1.8, -1.4, 0.04, 0, 0, 0.04};
vol96 = newv;
Cylinder(vol96) = {0.25, -1.8, -1.0, 0.04, 0, 0, 0.04};
vol97 = newv;
Cylinder(vol97) = {0.25, -1.8, -0.6, 0.04, 0, 0, 0.04};
vol98 = newv;
Cylinder(vol98) = {0.25, -1.8, -0.2, 0.04, 0, 0, 0.04};
vol99 = newv;
Cylinder(vol99) = {0.25, -1.8, 0.2, 0.04, 0, 0, 0.04};
vol100 = newv;
Cylinder(vol100) = {0.25, -1.8, 0.6, 0.04, 0, 0, 0.04};
vol101 = newv;
Cylinder(vol101) = {0.25, -1.8, 1.0, 0.04, 0, 0, 0.04};
vol102 = newv;
Cylinder(vol102) = {0.25, -1.8, 1.4, 0.04, 0, 0, 0.04};
vol103 = newv;
Cylinder(vol103) = {0.25, -1.8, 1.8, 0.04, 0, 0, 0.04};
vol104 = newv;
Box(vol104) = {0.25, -2.02, -2.02, 0.04, 4.04, 4.04};
vol105 = newv;
Box(vol105) = {0.25, 1.62, -1.98, 0.04, 0.36, 0.36};
bo2[] = BooleanDifference{ Volume{vol104}; Delete; } { Volume{vol105}; Delete;};
vol106 = newv;
Box(vol106) = {0.25, 1.62, -1.58, 0.04, 0.36, 0.36};
bo3[] = BooleanDifference{ Volume{bo2[]}; Delete; } { Volume{vol106}; Delete;};
vol107 = newv;
Box(vol107) = {0.25, 1.62, -1.18, 0.04, 0.36, 0.36};
bo4[] = BooleanDifference{ Volume{bo3[]}; Delete; } { Volume{vol107}; Delete;};
vol108 = newv;
Box(vol108) = {0.25, 1.62, -0.78, 0.04, 0.36, 0.36};
bo5[] = BooleanDifference{ Volume{bo4[]}; Delete; } { Volume{vol108}; Delete;};
vol109 = newv;
Box(vol109) = {0.25, 1.62, -0.38, 0.04, 0.36, 0.36};
bo6[] = BooleanDifference{ Volume{bo5[]}; Delete; } { Volume{vol109}; Delete;};
vol110 = newv;
Box(vol110) = {0.25, 1.62, 0.02, 0.04, 0.36, 0.36};
bo7[] = BooleanDifference{ Volume{bo6[]}; Delete; } { Volume{vol110}; Delete;};
vol111 = newv;
Box(vol111) = {0.25, 1.62, 0.42, 0.04, 0.36, 0.36};
bo8[] = BooleanDifference{ Volume{bo7[]}; Delete; } { Volume{vol111}; Delete;};
vol112 = newv;
Box(vol112) = {0.25, 1.62, 0.82, 0.04, 0.36, 0.36};
bo9[] = BooleanDifference{ Volume{bo8[]}; Delete; } { Volume{vol112}; Delete;};
vol113 = newv;
Box(vol113) = {0.25, 1.62, 1.22, 0.04, 0.36, 0.36};
bo10[] = BooleanDifference{ Volume{bo9[]}; Delete; } { Volume{vol113}; Delete;};
vol114 = newv;
Box(vol114) = {0.25, 1.62, 1.62, 0.04, 0.36, 0.36};
bo11[] = BooleanDifference{ Volume{bo10[]}; Delete; } { Volume{vol114}; Delete;};
vol115 = newv;
Box(vol115) = {0.25, 1.22, -1.98, 0.04, 0.36, 0.36};
bo12[] = BooleanDifference{ Volume{bo11[]}; Delete; } { Volume{vol115}; Delete;};
vol116 = newv;
Box(vol116) = {0.25, 1.22, -1.58, 0.04, 0.36, 0.36};
bo13[] = BooleanDifference{ Volume{bo12[]}; Delete; } { Volume{vol116}; Delete;};
vol117 = newv;
Box(vol117) = {0.25, 1.22, -1.18, 0.04, 0.36, 0.36};
bo14[] = BooleanDifference{ Volume{bo13[]}; Delete; } { Volume{vol117}; Delete;};
vol118 = newv;
Box(vol118) = {0.25, 1.22, -0.78, 0.04, 0.36, 0.36};
bo15[] = BooleanDifference{ Volume{bo14[]}; Delete; } { Volume{vol118}; Delete;};
vol119 = newv;
Box(vol119) = {0.25, 1.22, -0.38, 0.04, 0.36, 0.36};
bo16[] = BooleanDifference{ Volume{bo15[]}; Delete; } { Volume{vol119}; Delete;};
vol120 = newv;
Box(vol120) = {0.25, 1.22, 0.02, 0.04, 0.36, 0.36};
bo17[] = BooleanDifference{ Volume{bo16[]}; Delete; } { Volume{vol120}; Delete;};
vol121 = newv;
Box(vol121) = {0.25, 1.22, 0.42, 0.04, 0.36, 0.36};
bo18[] = BooleanDifference{ Volume{bo17[]}; Delete; } { Volume{vol121}; Delete;};
vol122 = newv;
Box(vol122) = {0.25, 1.22, 0.82, 0.04, 0.36, 0.36};
bo19[] = BooleanDifference{ Volume{bo18[]}; Delete; } { Volume{vol122}; Delete;};
vol123 = newv;
Box(vol123) = {0.25, 1.22, 1.22, 0.04, 0.36, 0.36};
bo20[] = BooleanDifference{ Volume{bo19[]}; Delete; } { Volume{vol123}; Delete;};
vol124 = newv;
Box(vol124) = {0.25, 1.22, 1.62, 0.04, 0.36, 0.36};
bo21[] = BooleanDifference{ Volume{bo20[]}; Delete; } { Volume{vol124}; Delete;};
vol125 = newv;
Box(vol125) = {0.25, 0.82, -1.98, 0.04, 0.36, 0.36};
bo22[] = BooleanDifference{ Volume{bo21[]}; Delete; } { Volume{vol125}; Delete;};
vol126 = newv;
Box(vol126) = {0.25, 0.82, -1.58, 0.04, 0.36, 0.36};
bo23[] = BooleanDifference{ Volume{bo22[]}; Delete; } { Volume{vol126}; Delete;};
vol127 = newv;
Box(vol127) = {0.25, 0.82, -1.18, 0.04, 0.36, 0.36};
bo24[] = BooleanDifference{ Volume{bo23[]}; Delete; } { Volume{vol127}; Delete;};
vol128 = newv;
Box(vol128) = {0.25, 0.82, -0.78, 0.04, 0.36, 0.36};
bo25[] = BooleanDifference{ Volume{bo24[]}; Delete; } { Volume{vol128}; Delete;};
vol129 = newv;
Box(vol129) = {0.25, 0.82, -0.38, 0.04, 0.36, 0.36};
bo26[] = BooleanDifference{ Volume{bo25[]}; Delete; } { Volume{vol129}; Delete;};
vol130 = newv;
Box(vol130) = {0.25, 0.82, 0.02, 0.04, 0.36, 0.36};
bo27[] = BooleanDifference{ Volume{bo26[]}; Delete; } { Volume{vol130}; Delete;};
vol131 = newv;
Box(vol131) = {0.25, 0.82, 0.42, 0.04, 0.36, 0.36};
bo28[] = BooleanDifference{ Volume{bo27[]}; Delete; } { Volume{vol131}; Delete;};
vol132 = newv;
Box(vol132) = {0.25, 0.82, 0.82, 0.04, 0.36, 0.36};
bo29[] = BooleanDifference{ Volume{bo28[]}; Delete; } { Volume{vol132}; Delete;};
vol133 = newv;
Box(vol133) = {0.25, 0.82, 1.22, 0.04, 0.36, 0.36};
bo30[] = BooleanDifference{ Volume{bo29[]}; Delete; } { Volume{vol133}; Delete;};
vol134 = newv;
Box(vol134) = {0.25, 0.82, 1.62, 0.04, 0.36, 0.36};
bo31[] = BooleanDifference{ Volume{bo30[]}; Delete; } { Volume{vol134}; Delete;};
vol135 = newv;
Box(vol135) = {0.25, 0.42, -1.98, 0.04, 0.36, 0.36};
bo32[] = BooleanDifference{ Volume{bo31[]}; Delete; } { Volume{vol135}; Delete;};
vol136 = newv;
Box(vol136) = {0.25, 0.42, -1.58, 0.04, 0.36, 0.36};
bo33[] = BooleanDifference{ Volume{bo32[]}; Delete; } { Volume{vol136}; Delete;};
vol137 = newv;
Box(vol137) = {0.25, 0.42, -1.18, 0.04, 0.36, 0.36};
bo34[] = BooleanDifference{ Volume{bo33[]}; Delete; } { Volume{vol137}; Delete;};
vol138 = newv;
Box(vol138) = {0.25, 0.42, -0.78, 0.04, 0.36, 0.36};
bo35[] = BooleanDifference{ Volume{bo34[]}; Delete; } { Volume{vol138}; Delete;};
vol139 = newv;
Box(vol139) = {0.25, 0.42, -0.38, 0.04, 0.36, 0.36};
bo36[] = BooleanDifference{ Volume{bo35[]}; Delete; } { Volume{vol139}; Delete;};
vol140 = newv;
Box(vol140) = {0.25, 0.42, 0.02, 0.04, 0.36, 0.36};
bo37[] = BooleanDifference{ Volume{bo36[]}; Delete; } { Volume{vol140}; Delete;};
vol141 = newv;
Box(vol141) = {0.25, 0.42, 0.42, 0.04, 0.36, 0.36};
bo38[] = BooleanDifference{ Volume{bo37[]}; Delete; } { Volume{vol141}; Delete;};
vol142 = newv;
Box(vol142) = {0.25, 0.42, 0.82, 0.04, 0.36, 0.36};
bo39[] = BooleanDifference{ Volume{bo38[]}; Delete; } { Volume{vol142}; Delete;};
vol143 = newv;
Box(vol143) = {0.25, 0.42, 1.22, 0.04, 0.36, 0.36};
bo40[] = BooleanDifference{ Volume{bo39[]}; Delete; } { Volume{vol143}; Delete;};
vol144 = newv;
Box(vol144) = {0.25, 0.42, 1.62, 0.04, 0.36, 0.36};
bo41[] = BooleanDifference{ Volume{bo40[]}; Delete; } { Volume{vol144}; Delete;};
vol145 = newv;
Box(vol145) = {0.25, 0.02, -1.98, 0.04, 0.36, 0.36};
bo42[] = BooleanDifference{ Volume{bo41[]}; Delete; } { Volume{vol145}; Delete;};
vol146 = newv;
Box(vol146) = {0.25, 0.02, -1.58, 0.04, 0.36, 0.36};
bo43[] = BooleanDifference{ Volume{bo42[]}; Delete; } { Volume{vol146}; Delete;};
vol147 = newv;
Box(vol147) = {0.25, 0.02, -1.18, 0.04, 0.36, 0.36};
bo44[] = BooleanDifference{ Volume{bo43[]}; Delete; } { Volume{vol147}; Delete;};
vol148 = newv;
Box(vol148) = {0.25, 0.02, -0.78, 0.04, 0.36, 0.36};
bo45[] = BooleanDifference{ Volume{bo44[]}; Delete; } { Volume{vol148}; Delete;};
vol149 = newv;
Box(vol149) = {0.25, 0.02, -0.38, 0.04, 0.36, 0.36};
bo46[] = BooleanDifference{ Volume{bo45[]}; Delete; } { Volume{vol149}; Delete;};
vol150 = newv;
Box(vol150) = {0.25, 0.02, 0.02, 0.04, 0.36, 0.36};
bo47[] = BooleanDifference{ Volume{bo46[]}; Delete; } { Volume{vol150}; Delete;};
vol151 = newv;
Box(vol151) = {0.25, 0.02, 0.42, 0.04, 0.36, 0.36};
bo48[] = BooleanDifference{ Volume{bo47[]}; Delete; } { Volume{vol151}; Delete;};
vol152 = newv;
Box(vol152) = {0.25, 0.02, 0.82, 0.04, 0.36, 0.36};
bo49[] = BooleanDifference{ Volume{bo48[]}; Delete; } { Volume{vol152}; Delete;};
vol153 = newv;
Box(vol153) = {0.25, 0.02, 1.22, 0.04, 0.36, 0.36};
bo50[] = BooleanDifference{ Volume{bo49[]}; Delete; } { Volume{vol153}; Delete;};
vol154 = newv;
Box(vol154) = {0.25, 0.02, 1.62, 0.04, 0.36, 0.36};
bo51[] = BooleanDifference{ Volume{bo50[]}; Delete; } { Volume{vol154}; Delete;};
vol155 = newv;
Box(vol155) = {0.25, -0.38, -1.98, 0.04, 0.36, 0.36};
bo52[] = BooleanDifference{ Volume{bo51[]}; Delete; } { Volume{vol155}; Delete;};
vol156 = newv;
Box(vol156) = {0.25, -0.38, -1.58, 0.04, 0.36, 0.36};
bo53[] = BooleanDifference{ Volume{bo52[]}; Delete; } { Volume{vol156}; Delete;};
vol157 = newv;
Box(vol157) = {0.25, -0.38, -1.18, 0.04, 0.36, 0.36};
bo54[] = BooleanDifference{ Volume{bo53[]}; Delete; } { Volume{vol157}; Delete;};
vol158 = newv;
Box(vol158) = {0.25, -0.38, -0.78, 0.04, 0.36, 0.36};
bo55[] = BooleanDifference{ Volume{bo54[]}; Delete; } { Volume{vol158}; Delete;};
vol159 = newv;
Box(vol159) = {0.25, -0.38, -0.38, 0.04, 0.36, 0.36};
bo56[] = BooleanDifference{ Volume{bo55[]}; Delete; } { Volume{vol159}; Delete;};
vol160 = newv;
Box(vol160) = {0.25, -0.38, 0.02, 0.04, 0.36, 0.36};
bo57[] = BooleanDifference{ Volume{bo56[]}; Delete; } { Volume{vol160}; Delete;};
vol161 = newv;
Box(vol161) = {0.25, -0.38, 0.42, 0.04, 0.36, 0.36};
bo58[] = BooleanDifference{ Volume{bo57[]}; Delete; } { Volume{vol161}; Delete;};
vol162 = newv;
Box(vol162) = {0.25, -0.38, 0.82, 0.04, 0.36, 0.36};
bo59[] = BooleanDifference{ Volume{bo58[]}; Delete; } { Volume{vol162}; Delete;};
vol163 = newv;
Box(vol163) = {0.25, -0.38, 1.22, 0.04, 0.36, 0.36};
bo60[] = BooleanDifference{ Volume{bo59[]}; Delete; } { Volume{vol163}; Delete;};
vol164 = newv;
Box(vol164) = {0.25, -0.38, 1.62, 0.04, 0.36, 0.36};
bo61[] = BooleanDifference{ Volume{bo60[]}; Delete; } { Volume{vol164}; Delete;};
vol165 = newv;
Box(vol165) = {0.25, -0.78, -1.98, 0.04, 0.36, 0.36};
bo62[] = BooleanDifference{ Volume{bo61[]}; Delete; } { Volume{vol165}; Delete;};
vol166 = newv;
Box(vol166) = {0.25, -0.78, -1.58, 0.04, 0.36, 0.36};
bo63[] = BooleanDifference{ Volume{bo62[]}; Delete; } { Volume{vol166}; Delete;};
vol167 = newv;
Box(vol167) = {0.25, -0.78, -1.18, 0.04, 0.36, 0.36};
bo64[] = BooleanDifference{ Volume{bo63[]}; Delete; } { Volume{vol167}; Delete;};
vol168 = newv;
Box(vol168) = {0.25, -0.78, -0.78, 0.04, 0.36, 0.36};
bo65[] = BooleanDifference{ Volume{bo64[]}; Delete; } { Volume{vol168}; Delete;};
vol169 = newv;
Box(vol169) = {0.25, -0.78, -0.38, 0.04, 0.36, 0.36};
bo66[] = BooleanDifference{ Volume{bo65[]}; Delete; } { Volume{vol169}; Delete;};
vol170 = newv;
Box(vol170) = {0.25, -0.78, 0.02, 0.04, 0.36, 0.36};
bo67[] = BooleanDifference{ Volume{bo66[]}; Delete; } { Volume{vol170}; Delete;};
vol171 = newv;
Box(vol171) = {0.25, -0.78, 0.42, 0.04, 0.36, 0.36};
bo68[] = BooleanDifference{ Volume{bo67[]}; Delete; } { Volume{vol171}; Delete;};
vol172 = newv;
Box(vol172) = {0.25, -0.78, 0.82, 0.04, 0.36, 0.36};
bo69[] = BooleanDifference{ Volume{bo68[]}; Delete; } { Volume{vol172}; Delete;};
vol173 = newv;
Box(vol173) = {0.25, -0.78, 1.22, 0.04, 0.36, 0.36};
bo70[] = BooleanDifference{ Volume{bo69[]}; Delete; } { Volume{vol173}; Delete;};
vol174 = newv;
Box(vol174) = {0.25, -0.78, 1.62, 0.04, 0.36, 0.36};
bo71[] = BooleanDifference{ Volume{bo70[]}; Delete; } { Volume{vol174}; Delete;};
vol175 = newv;
Box(vol175) = {0.25, -1.18, -1.98, 0.04, 0.36, 0.36};
bo72[] = BooleanDifference{ Volume{bo71[]}; Delete; } { Volume{vol175}; Delete;};
vol176 = newv;
Box(vol176) = {0.25, -1.18, -1.58, 0.04, 0.36, 0.36};
bo73[] = BooleanDifference{ Volume{bo72[]}; Delete; } { Volume{vol176}; Delete;};
vol177 = newv;
Box(vol177) = {0.25, -1.18, -1.18, 0.04, 0.36, 0.36};
bo74[] = BooleanDifference{ Volume{bo73[]}; Delete; } { Volume{vol177}; Delete;};
vol178 = newv;
Box(vol178) = {0.25, -1.18, -0.78, 0.04, 0.36, 0.36};
bo75[] = BooleanDifference{ Volume{bo74[]}; Delete; } { Volume{vol178}; Delete;};
vol179 = newv;
Box(vol179) = {0.25, -1.18, -0.38, 0.04, 0.36, 0.36};
bo76[] = BooleanDifference{ Volume{bo75[]}; Delete; } { Volume{vol179}; Delete;};
vol180 = newv;
Box(vol180) = {0.25, -1.18, 0.02, 0.04, 0.36, 0.36};
bo77[] = BooleanDifference{ Volume{bo76[]}; Delete; } { Volume{vol180}; Delete;};
vol181 = newv;
Box(vol181) = {0.25, -1.18, 0.42, 0.04, 0.36, 0.36};
bo78[] = BooleanDifference{ Volume{bo77[]}; Delete; } { Volume{vol181}; Delete;};
vol182 = newv;
Box(vol182) = {0.25, -1.18, 0.82, 0.04, 0.36, 0.36};
bo79[] = BooleanDifference{ Volume{bo78[]}; Delete; } { Volume{vol182}; Delete;};
vol183 = newv;
Box(vol183) = {0.25, -1.18, 1.22, 0.04, 0.36, 0.36};
bo80[] = BooleanDifference{ Volume{bo79[]}; Delete; } { Volume{vol183}; Delete;};
vol184 = newv;
Box(vol184) = {0.25, -1.18, 1.62, 0.04, 0.36, 0.36};
bo81[] = BooleanDifference{ Volume{bo80[]}; Delete; } { Volume{vol184}; Delete;};
vol185 = newv;
Box(vol185) = {0.25, -1.58, -1.98, 0.04, 0.36, 0.36};
bo82[] = BooleanDifference{ Volume{bo81[]}; Delete; } { Volume{vol185}; Delete;};
vol186 = newv;
Box(vol186) = {0.25, -1.58, -1.58, 0.04, 0.36, 0.36};
bo83[] = BooleanDifference{ Volume{bo82[]}; Delete; } { Volume{vol186}; Delete;};
vol187 = newv;
Box(vol187) = {0.25, -1.58, -1.18, 0.04, 0.36, 0.36};
bo84[] = BooleanDifference{ Volume{bo83[]}; Delete; } { Volume{vol187}; Delete;};
vol188 = newv;
Box(vol188) = {0.25, -1.58, -0.78, 0.04, 0.36, 0.36};
bo85[] = BooleanDifference{ Volume{bo84[]}; Delete; } { Volume{vol188}; Delete;};
vol189 = newv;
Box(vol189) = {0.25, -1.58, -0.38, 0.04, 0.36, 0.36};
bo86[] = BooleanDifference{ Volume{bo85[]}; Delete; } { Volume{vol189}; Delete;};
vol190 = newv;
Box(vol190) = {0.25, -1.58, 0.02, 0.04, 0.36, 0.36};
bo87[] = BooleanDifference{ Volume{bo86[]}; Delete; } { Volume{vol190}; Delete;};
vol191 = newv;
Box(vol191) = {0.25, -1.58, 0.42, 0.04, 0.36, 0.36};
bo88[] = BooleanDifference{ Volume{bo87[]}; Delete; } { Volume{vol191}; Delete;};
vol192 = newv;
Box(vol192) = {0.25, -1.58, 0.82, 0.04, 0.36, 0.36};
bo89[] = BooleanDifference{ Volume{bo88[]}; Delete; } { Volume{vol192}; Delete;};
vol193 = newv;
Box(vol193) = {0.25, -1.58, 1.22, 0.04, 0.36, 0.36};
bo90[] = BooleanDifference{ Volume{bo89[]}; Delete; } { Volume{vol193}; Delete;};
vol194 = newv;
Box(vol194) = {0.25, -1.58, 1.62, 0.04, 0.36, 0.36};
bo91[] = BooleanDifference{ Volume{bo90[]}; Delete; } { Volume{vol194}; Delete;};
vol195 = newv;
Box(vol195) = {0.25, -1.98, -1.98, 0.04, 0.36, 0.36};
bo92[] = BooleanDifference{ Volume{bo91[]}; Delete; } { Volume{vol195}; Delete;};
vol196 = newv;
Box(vol196) = {0.25, -1.98, -1.58, 0.04, 0.36, 0.36};
bo93[] = BooleanDifference{ Volume{bo92[]}; Delete; } { Volume{vol196}; Delete;};
vol197 = newv;
Box(vol197) = {0.25, -1.98, -1.18, 0.04, 0.36, 0.36};
bo94[] = BooleanDifference{ Volume{bo93[]}; Delete; } { Volume{vol197}; Delete;};
vol198 = newv;
Box(vol198) = {0.25, -1.98, -0.78, 0.04, 0.36, 0.36};
bo95[] = BooleanDifference{ Volume{bo94[]}; Delete; } { Volume{vol198}; Delete;};
vol199 = newv;
Box(vol199) = {0.25, -1.98, -0.38, 0.04, 0.36, 0.36};
bo96[] = BooleanDifference{ Volume{bo95[]}; Delete; } { Volume{vol199}; Delete;};
vol200 = newv;
Box(vol200) = {0.25, -1.98, 0.02, 0.04, 0.36, 0.36};
bo97[] = BooleanDifference{ Volume{bo96[]}; Delete; } { Volume{vol200}; Delete;};
vol201 = newv;
Box(vol201) = {0.25, -1.98, 0.42, 0.04, 0.36, 0.36};
bo98[] = BooleanDifference{ Volume{bo97[]}; Delete; } { Volume{vol201}; Delete;};
vol202 = newv;
Box(vol202) = {0.25, -1.98, 0.82, 0.04, 0.36, 0.36};
bo99[] = BooleanDifference{ Volume{bo98[]}; Delete; } { Volume{vol202}; Delete;};
vol203 = newv;
Box(vol203) = {0.25, -1.98, 1.22, 0.04, 0.36, 0.36};
bo100[] = BooleanDifference{ Volume{bo99[]}; Delete; } { Volume{vol203}; Delete;};
vol204 = newv;
Box(vol204) = {0.25, -1.98, 1.62, 0.04, 0.36, 0.36};
bo101[] = BooleanDifference{ Volume{bo100[]}; Delete; } { Volume{vol204}; Delete;};//+
Physical Surface("anode") = {3, 2, 1};
Physical Surface("cathode") = {4, 5, 6};
Physical Surface("walls") = {7, 8, 9, 10};
Physical Surface("pixels") = {11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310};
Physical Surface("grid") = {311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329,330,331,332,333,334,335,336,337,338,339,340,341,342,343,344,345,346,347,348,349,350,351,352,353,354,355,356,357,358,359,360,361,362,363,364,365,366,367,368,369,370,371,372,373,374,375,376,377,378,379,380,381,382,383,384,385,386,387,388,389,390,391,392,393,394,395,396,397,398,399,400,401,402,403,404,405,406,407,408,409,410,411,412,413,414,415,416,417,418,419,420,421,422,423,424,425,426,427,428,429,430,431,432,433,434,435,436,437,438,439,440,441,442,443,444,445,446,447,448,449,450,451,452,453,454,455,456,457,458,459,460,461,462,463,464,465,466,467,468,469,470,471,472,473,474,475,476,477,478,479,480,481,482,483,484,485,486,487,488,489,490,491,492,493,494,495,496,497,498,499,500,501,502,503,504,505,506,507,508,509,510,511,512,513,514,515,516,517,518,519,520,521,522,523,524,525,526,527,528,529,530,531,532,533,534,535,536,537,538,539,540,541,542,543,544,545,546,547,548,549,550,551,552,553,554,555,556,557,558,559,560,561,562,563,564,565,566,567,568,569,570,571,572,573,574,575,576,577,578,579,580,581,582,583,584,585,586,587,588,589,590,591,592,593,594,595,596,597,598,599,600,601,602,603,604,605,606,607,608,609,610,611,612,613,614,615,616,617,618,619,620,621,622,623,624,625,626,627,628,629,630,631,632,633,634,635,636,637,638,639,640,641,642,643,644,645,646,647,648,649,650,651,652,653,654,655,656,657,658,659,660,661,662,663,664,665,666,667,668,669,670,671,672,673,674,675,676,677,678,679,680,681,682,683,684,685,686,687,688,689,690,691,692,693,694,695,696,697,698,699,700,701,702,703,704,705,706,707,708,709,710,711,712,713,714,715,716};


