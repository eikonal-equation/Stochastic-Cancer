%This function uses bilinear interpolation to compute the value function
%for interior points in the domain
%v11,v12,v21,v22(input) : function values of at four corners of the box
%xl(input): left x coordinate of the box
%xr(input): right x coordinate of the box
%yt(input): top y coordinate of the box
%yb(input): bottom y coordinate of the box
%xloc(input): x coordinate of the query point
%yloc(input): y coordinate of the query point
%dx(input): dx
%dy(input): dy
%x0(input): x coordiante of the first point of the grid
%y0(input): y coordiante of the first point of the grid
%N(input) : dimension of the input matrix minus 1
%val(output): interpolated value
%
% Author: MingYi Wang, Cornell University
% Last Modified: 06/18/2022
%
function val=BilinearInterp(v11,v12,v21,v22,xl,xr,yt,yb,xloc,yloc,dx,dy)
%bilinear interpolation
val= (yt-yloc)/dy*( (xr-xloc)/dx*v11 + (xloc-xl)/dx*v12)...
    +(yloc-yb)/dy*( (xr-xloc)/dx*v21+(xloc-xl)/dx*v22);
end