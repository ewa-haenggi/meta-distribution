% Algorithm for recovering the accurate cdf via a moment sequence with an analytic expression 
% The GP and FL method are used 
clear all
% input: analytic expression, point of interest (default is
% linspace(0,1,21)), and the number of integer moments
% analytic expression 
g = @(n) 1./hypergeom([n,-1/2],[1/2],-1/2); %1./(n.^1+1);
f = @(n) arrayfun(g, n);
f = @(n) exp(-sqrt(n));
% points of interest
pts = linspace(0,1,21);
recy = 0*pts;
% number of used moments if GP is not valid
N = 20; 

% GP
% GP initialization
stpsz = 0.1;
lowlim = 1e-15;
uplim = 1000;
methodVar = GPClass;
methodVar = methodVar.init(stpsz, lowlim, uplim,f);
idx = 1;
while idx <= length(pts)
    pt = pts(idx);
    value = methodVar.value( pt);
    if value > 1 || value < 0
        stpsz = methodVar.stpsz/3;
        lowlim = methodVar.lowlim;
        uplim = methodVar.uplim*3;
        if uplim/stpsz > 10000000 % memory requirement
            flag = 0;
            xnew = [0 1];
            ynew = [0 1];
            break
        end
        %                 reinitialize GPClass variable
        methodVar = GPClass;
        methodVar = methodVar.init(stpsz, lowlim, uplim,f);
        idx = 1;
    else
        recy(idx) = value;
        idx = idx + 1;
    end 
end

% FL is GP is not valid
if ~flag
    moment = f((1:N)); 
    order = length(moment) - 1;
    methodVar = FLClass;
    methodVar = methodVar.init(moment,order);
    recy = methodVar.value(pts);    
end
 

