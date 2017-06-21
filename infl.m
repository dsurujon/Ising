% [input] X: samples
% [input] u: node index u
% [input] i: node index i
% [input] xs: configuration of set S
% [input] S: subset of nodes 
% [output] v: conditional influence of node i on node u 
function v=infl(X,u,i,S,xs)
    v=phat_cond(X,u,1,[S, i],[xs, 1])-phat_cond(X,u,1,[S, i],[xs, -1]);

end