%This function uses linear interpolation to compute the value function
%v(input) : (N+1)x1 vector of values of sample points.
%xloc(input): coordinate of the query point
%k(input): the first index such that xloc<=x, where x is the vector of
%          sample points
%h(input): dx
%x0(input): first point of the grid
%
%val(output): interpolated value
%
% Author: MingYi Wang, Cornell University
% Last modified: 12/2023
%
function val=LinearInterp(v,k,xloc,h,x0) 

%use the stencil [x(k-1),x(k)] to compute the value. Note k>=2
%Here x(k)=x0+(k-1)*h
val=(((x0+(k-1)*h)-xloc)*v(k-1)+(xloc-(x0+(k-2)*h))*v(k))/h;


end