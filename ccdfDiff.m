function [dv,dh] = ccdfDiff(x1,y1, x2,y2, n)
% CCDFDIFF  Calculate horizontal and vertical difference of two sets of
% functions. 
%
%   [DV,DH] = CCDFDIFF(X1,Y1,X2,Y2,N) calculate horizontal and vertical
%   difference of (x1,y1) and (x2,y2) according to n points.
%   x1 and x2 are increasing.
%   y1 and y2 are non-increasing.
%   Their entries are in [0,1].

dv = 0; dh = 0;
x1 = real(x1);
x2 = real(x2);
y1 = real(y1);
y2 = real(y2);
y1 = double(y1);
y2 = double(y2);

if n <= 2
    fprintf('n <3. Please choose another n.\n');
    return;
end 
if size(x1,1) == 1
    x1 = x1';
end
if size(y1,1) == 1
    y1 = y1';
end
if size(x2,1) == 1
    x2 = x2';
end
if size(y2,1) == 1
    y2= y2';
end
if sum(diff(x1)<=0) >0 || sum(diff(x2)<=0) >0
    fprintf('x1/x2 is not increasing. Please choose another x1/x2.\n')
    return;
end

if sum(diff(y1)>0) >0 || sum(diff(y2)>0) >0
    fprintf('y1/y2 is not decreasing.\n')
end

if (length(x1) ~= length(y1)) || (length(x2)~=length(y2))
   fprintf("(x1,y1) or (x2,y2) doesn't match.\n")
   return;
end

if isempty(x1) || isempty(x2)
    fprintf('x1/x2 should contain at leat one element for valid calculation.\n')
    return;
end

 


% vertical distance 
xsampledPoints = linspace(0,1,n)';  
method = 'pchip';
vy1 = interp1(x1,y1,xsampledPoints,method);
vy2 = interp1(x2,y2,xsampledPoints,method);
dv = vy1 - vy2;



% horizontal distance 
ysampledPoints = linspace(0,1,n)';  
% only keep the unique values of y 
% F^{-1}(y) = min_{F(x) = y} x
[y1C,ia,ic] = unique(y1);
x1C = x1(ia);
[y2C,ia,ic] = unique(y2);
x2C = x2(ia);

hx1 = interp1(y1C,x1C,ysampledPoints,method);   
hx2 = interp1(y2C,x2C,ysampledPoints,method); 
dh = hx1 - hx2;