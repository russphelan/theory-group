function [ transl_func_obj ] = translate_x( func_obj,transl_amt )
%takes a function, and translates it in the independent variable by amt. 

transl_func_obj = func_obj; %for speed, allocates memory

for i=1:length(func_obj)
    transl_func_obj(2,i) = func_obj(2,i) + transl_amt;
end
end

