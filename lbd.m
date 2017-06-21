% [input] X: samples
% [input] i: node index i
% [input] xs: configuration of set S
% [input] S: subset of nodes 
% [output] l: lambda 
function l=lbd(X,i,S,xs)
    l=2*phat_cond(X,i,1,S,xs)*phat_cond(X,i,-1,S,xs);
end