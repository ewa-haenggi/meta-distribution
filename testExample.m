% test example for using the classes - reconstruction methods
% BM, FL, CM, FC, FJ, ME, GP
% initialization - moments (order from 1 to N), point of interest in (0,1),
% and moment expression (as function handle)
N = 10;
moment = 1./(1+(1:N)');
momentFunc = @(n) 1./(1+n);
pt = linspace(0,1,21);

% BM 
n = length(moment);
b = 16;
methodVar = BMClass;
methodVar = methodVar.init(n,b,moment); 
vbm = methodVar.value(pt);

% FL 
order = length(moment) - 1;
methodVar = FLClass;
methodVar = methodVar.init(moment,order, floor(0.9*order));
vfl = methodVar.value(pt);

% CM
order = length(moment);
methodVar = CMClass;
methodVar = methodVar.init(order);
bounds = [0*pt; 0*pt];
for idx = 1:length(pt)
    bounds(:,idx) = methodVar.CMBounds(moment, pt(idx));
end
vcm = 1/2*sum(bounds,1);

% FC
order = length(moment) - 1;
methodVar = FCClass;
methodVar = methodVar.init(moment,order, floor(0.9*order));
vfc = methodVar.value(pt);

% FJ
m1 = moment(1);
m2 = moment(2);
alpha = (m1-m2)*(1-m1)/(m2-m1.^2)-1;
beta = (alpha+1)*m1/(1-m1)-1;
order = length(moment);
methodVar = FJClass;
methodVar = methodVar.init(moment, order);
vfj = 0*pt;
for idx = 1:length(pt)
    vfj(idx) = methodVar.value(pt(idx));
end


% ME
order = length(moment);
methodVar = MEClass;
methodVar = methodVar.init(moment, order);
vme = 0*pt;
for idx = 1:length(pt)
    vme(idx) = methodVar.value(pt(idx));
end

% GP
stepsz = 0.1;
lowlim = 1e-10;
uplim = 5e1;
methodVar = GPClass;
methodVar = methodVar.init(stepsz, lowlim, uplim, momentFunc);
vgp = 0*pt;
for idx = 1:length(pt)
    vgp(idx) = methodVar.value(pt(idx));
end

