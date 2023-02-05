function [xnew,ynew] = revisePtv3(x, y)
% REVISEPTV3 Revise points to make it as a cdf (without a filter, which is for
% functional approximations, and with ``making the curve start from (0,0) 
% and end at (1,1)'', which is for distance evaluations).  
% 
%   [XNEW,YNEW] = REVISEPTV3(X,Y) revise points based on (x,y) to make it as 
%   a cdf, i.e., x and y are increasing from 0 to 1, by revising points that
%   exceeding [0,1] to normal values 0 or 1, replacing points (y) that are smaller than the
%   previous one by the previous one, and finally making (x,y) a (monotone) curve 
%   starting from (0,0) and ending at (1,1).

y = reshape(y,1, size(y,1)*size(y,2));
x = reshape(x,1, size(x,1)*size(x,2)); 


y(1<y) = 1;
y(y<0) = 0;

% FOR GP, Y = NAN WITH X = 0 WE NEED TO MANUALLY CHANGE IT
if isnan(y(1)) && x(1) == 0
    y(1) = 0;
end

ynew = y;

for i = 1:length(ynew)
    if i >= 2
        if ynew(i) < ynew(i-1)
            ynew(i) = ynew(i-1); 
        end
    end
end
xnew = x;

%% FOR DISTANCE EVALUATION (BOTH VERTICAL AND HORIZONTAL)

% IF THE CURVE STARTS FROM (0,X) (NOT GOOD FOR HORIZONTAL DISTANCE EVAL.)

if (xnew(1) == 0 && ynew(1) ~= 0)
    xnew = [0, 1e-15*min(xnew(2:end)), xnew(2:end)];
    ynew = [0, ynew];
end
   
% IF THE CURVE STARTS FROM (X,0) (NOT GOOD FOR VERTICAL DISTANCE EVAL.)    

if (xnew(1) ~= 0 && ynew(1) == 0)
    xnew = [0, xnew];
    ynew = [0, 1e-15*min(ynew(2:end)), ynew(2:end)];
end 

% IF THE CURVE ENDS AT (X,1) (NOT GOOD FOR VERTICAL DISTANCE EVAL.)

if (xnew(end) ~= 1 && ynew(end) == 1)
    xnew = [xnew,1];
    ynew = [ynew(1:end-1), 1 - (1-max(ynew(1:end-1)))*1e-10,1];
end

% IF THE CURVE ENDS AT (1,X) (NOT GOOD FOR HORIZONTAL DISTANCE EVAL.)

if (xnew(end) == 1 && ynew(end) ~= 1)
    xnew = [xnew(1:end-1),1 - (1-max(xnew(1:end-1)))*1e-10,1];
    ynew = [ynew, 1];
end

% IF THE POINTS DO NOT START FROM (0,0), ADD IT

if (xnew(1) ~= 0 && ynew(1)~= 0)
   tmp = zeros(length(xnew)+1,1);
   tmp(1) = 0;
   tmp(2:end) = xnew;
   xnew = tmp;
   tmp = zeros(length(ynew)+1,1);
   tmp(1) = 0;
   tmp(2:end) = ynew;
   ynew = tmp;
end
% IF THE POINTS DO NOT END AT (1,1), ADD IT

if (xnew(end) ~= 1 && ynew(end)~= 1)
   tmp = zeros(length(xnew)+1,1);
   tmp(end) = 1;
   tmp(1:end-1) = xnew;
   xnew = tmp;
   tmp = zeros(length(ynew)+1,1);
   tmp(end) = 1;
   tmp(1:end-1) = ynew;
   ynew = tmp;
end

xnew = reshape(xnew,size(xnew,1)*size(xnew,2),1);
ynew = reshape(ynew,size(ynew,1)*size(ynew,2),1);

if length(xnew) ~= length(ynew)
    1 
end 
 
end