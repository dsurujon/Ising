% [input] X: samples
% [input] xs: configuration of set S
% [input] S: subset of nodes 
% [output] p: empirical probability of observing set S in configuration xs

function p = phat(X,S,xs)
    if (isempty(S)==1)
        p=1;
    else
        n=size(X,1);
        subset=X(:,S);
        B=repmat(xs,n,1);
        p=sum(all(B==subset,2))/n;
    end
end