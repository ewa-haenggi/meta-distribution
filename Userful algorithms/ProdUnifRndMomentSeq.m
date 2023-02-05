% Algorithm for producing uniformly random moment sequences
clear all;
tic;
% specify length
n = 10; % the number of sequences
k = 20; % the number of moments for each sequence

momentSeq = zeros(k,n);
for i = 1:k
    momentVar = canonicalMomentClass;
    momentSeq(i,:) = momentVar.gen(n);
end
toc;