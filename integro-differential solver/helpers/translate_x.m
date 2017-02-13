%Author: Russell J. Phelan 
%Date: 8/12/16

%I would like to thank John Donoghue, Basem El-Menoufi, Panayotis Kevrekidis, and William ?Bill? Barnes 
%for useful conversations and inspiration related to this project. This work has been supported in part 
%by the National Science Foundation under grants NSF PHY15-20292 and NSF PHY12-25915.

function [ transl_func_obj ] = translate_x( func_obj,transl_amt )
%takes a function, and translates it in the independent variable by amt. 

transl_func_obj = func_obj; %for speed, allocates memory

for i=1:length(func_obj)
    transl_func_obj(2,i) = func_obj(2,i) + transl_amt;
end
end

