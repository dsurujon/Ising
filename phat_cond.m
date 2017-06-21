% [input] X: samples
% [input] xs: configuration of set S
% [input] xt: configuration of set T
% [input] S: subset of nodes 
% [input] T: subset of nodes 
% [output] p: empirical probability of observing set T in configuration xt
% given set S is in configuration xs

function p = phat_cond(X,T,xt,S,xs)
    p=phat(X,[S, T],[xs, xt])/phat(X,S,xs);
end