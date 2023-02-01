% test example for different decay types
% sub-power, power, soft power, intermediate, exp, soft exp

% number of moments
N = 10; 

% randomize parameters
s = rand()*10;
k = randi([1 10]);  
a = rand(1,k)*10;
b = rand(1,k)*10;
c = rand(1,k);
c = c/sum(c);

% initialize class object and generate moment
% subpower
decayVar = mixedSubPowerDecay;
decayVar = decayVar.init(s,a,c);
moment = decayVar.gen(1:N);

% power
decayVar = mixedPowerDecay;
decayVar = decayVar.init(s,a,c);
moment = decayVar.gen(1:N);

% softPower
decayVar = mixedSoftPowerDecay;
decayVar = decayVar.init(s,a,b,c);
moment = decayVar.gen(1:N);

% intermediate
decayVar = mixedIntmDecay;
decayVar = decayVar.init(s,a,c);
moment = decayVar.gen(1:N);

% exp 
t = rand();
s = -log(t);
if k > 1
    s = [s, -log(rand(1, k-1)*exp(-s))];
end
decayVar = mixedExpDecay;
decayVar = decayVar.init(s,c );
moment = decayVar.gen(1:N);

% softExp
t = rand();
s = -log(t);
decayVar = mixedSoftExpDecay;
decayVar = decayVar.init(s,a,c);
moment = decayVar.gen(1:N);




