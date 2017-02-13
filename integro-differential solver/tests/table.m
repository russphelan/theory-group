%Author: Russell J. Phelan 
%Date: 9-9-16

%I would like to thank John Donoghue, Basem El-Menoufi, Panayotis Kevrekidis, and William ?Bill? Barnes 
%for useful conversations and inspiration related to this project. This work has been supported in part 
%by the National Science Foundation under grants NSF PHY15-20292 and NSF PHY12-25915.

row_scale_factor = transpose(scale_factor);
row_scale_1deriv = transpose(scale_1deriv);
row_scale_2deriv = transpose(scale_2deriv);
row_scale_3deriv = transpose(scale_3deriv);
row_area_matrix = transpose(area_matrix);

table(row_scale_factor(1999:2011,1),row_scale_1deriv(1999:2011,1),row_scale_2deriv(1999:2011,1),...
    row_scale_3deriv(1999:2011,1),row_area_matrix(1999:2011))