% Learning the Ising model
% Given a set of gene presence/absence data, find the underlying Markov 
% random field
% Adapted from Bresler arXiv:1411.6156v2 [cs.LG]

%20 genomes, 10 genes
x=2*(rand(20,10)>0.5)-1;

%what's the right t???
stest=learnNbhd(1,x,0.05);
disp(stest)

%%
x=2*(rand(100,20)>0.5)-1;
%%
d=5;alpha=0.01;beta=0.05;h=beta;
delta=0.5*exp(-beta*d-h);
t=(alpha^2)*(delta^(4*d+1))/(16*d*beta);
%%
edges=zeros(20,20);
for i=1:20
    nbd=learnNbhd(i,x,t);
    edges(i,nbd)=1;
    disp([i,nbd]);
end
%this is not symmetric... neighborhoods of pairs of nodes don't overlap
%completely. 
imagesc(edges);