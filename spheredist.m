%% Spherical Distribution Function
% This function returns an array of points that are uniformly distributed
% within the volume of a sphere of dimension dim and unit radius. 
function points = spheredist(N,dim,varargin)
r = (rand(N,1))^(1/dim);
x = 0:1e-3:1;
p = sin(x);









end