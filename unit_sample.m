function [dn, n] = unit_sample(number_of_samples) 

dn = zeros(1, number_of_samples);
dn(1) = 1;
n = (0:1:number_of_samples-1);