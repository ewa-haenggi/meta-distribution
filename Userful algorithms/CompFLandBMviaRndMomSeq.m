% Algorithm for comparing the FL and BM method via a random moment sequence with specific decay rate

% This is the raw version of FL and BM, without any interpolation and mapping
clear all;
% number of moments
n = 10;
% methods
methodName = ["FL", "BM"]; 
%% decay variable
% specify decay type and parameters
decayName = ["subPower", "power", "softPower", "intermediate", "exp", "softExp"]; 
decaySeq = randi([1 length(decayName)]);
decay = decayName(decaySeq);
s = rand()*10;
k = randi([1 10]);  
a = rand(1,k)*10;
b = rand(1,k)*10;
c = rand(1,k);
c = c/sum(c);
if strcmpi(decay, "power") == 1
    decayVar = mixedPowerDecay; 
    decayVar = decayVar.init(s,a,c);
elseif strcmpi(decay, "softPower") == 1
    decayVar = mixedSoftPowerDecay;
    decayVar = decayVar.init(s,a,b,c);
elseif strcmpi(decay, "intermediate") == 1
    decayVar = mixedIntmDecay;
    decayVar = decayVar.init(s,a,c);
elseif strcmpi(decay, "exp") == 1
    t = rand();
    s = -log(t);
    if k > 1
        s = [s, -log(rand(1, k-1)*exp(-s))];
    end
    decayVar = mixedExpDecay;
    decayVar = decayVar.init(s,c );
elseif strcmpi(decay, "softExp") == 1
    t = rand();
    s = -log(t);
    decayVar = mixedSoftExpDecay;
    decayVar = decayVar.init(s,a,c);
elseif strcmpi(decay, "subpower") == 1
    decayVar = mixedSubPowerDecay;
    decayVar = decayVar.init(s,a,c);
end

moment = decayVar.gen((1:n));

pts = linspace(0,1,21);
recy = zeros(length(methodName), length(pts));
for i = 1:length(methodName)
    method = methodName(i);
    if strcmpi(method,"FL")
        order = length(moment) - 1; 
        methodVar = FLClass;
        methodVar = methodVar.init(moment,order);
    elseif strcmpi(method,"BM")
        n = length(moment); % 
        b = 16;
        methodVar = BMClass;
        methodVar = methodVar.init(n,b,moment); 
    end
    recy(i,:) = methodVar.value(pts); 
end

plot(pts,recy(1,:), "-x", 'DisplayName',methodName(1));
hold on;
plot(pts,recy(2,:),'-d', 'DisplayName',methodName(2));
legend()
                    