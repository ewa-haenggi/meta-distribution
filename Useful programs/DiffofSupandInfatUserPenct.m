% Algorithm for plotting the difference of the supremum and 
% the infimum at the 95\% point as a function of the number of moments
clear all;

% input: moment and the point of interest 
n = 15;
moment = 1./(1+(1:n)').^(2);
percentile = 95; % user percentile
xsampledPoint = 1 - percentile/100;

% output: distance vs n
dist = 0*moment; 
for i = 1:n
truncMoment = moment(1:i);
pt = linspace(0,1,21);
methodVar = CMClass;
methodVar = methodVar.init(length(truncMoment));

infs = 0*pt;
sups = 0*pt;

for j = 1:length(pt)
    bounds = methodVar.CMBounds(truncMoment, pt(j));
    infs(j) = bounds(1);
    sups(j) = bounds(2);
end
[x1,y1] = revisePtv3(pt, infs);
[x2,y2] = revisePtv3(pt, sups);
figure(1)
plot(x1,y1,'r', x2,y2, 'b');

% take unique values
[y1C,ia,ic] = unique(y1);
x1C = x1(ia);
[y2C,ia,ic] = unique(y2);
x2C = x2(ia);

% interpolation for inverse function (quantile function)
method = 'pchip';
ix1 = interp1(y1C,x1C,xsampledPoint,method);
ix2 = interp1(y2C,x2C,xsampledPoint,method);

d = abs(ix1 - ix2);
dist(i) = d;
end

figure(2)
plot(1:n,dist)
grid on

