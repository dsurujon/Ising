% [input] X: samples
% [input] u: node index u
% [input] i: node index i
% [input] S: subset of nodes 
% [output] v: average conditional influence of node i on node u 
function v=avg_infl(X,u,i,S)
    v=0;
    %if (isempty(S)==1)
    %    v=phat_cond(X,u,1,i,1)-phat_cond(X,u,1,i,-1);
    %elseif (any(i==S)==1)
    %    return;
    %else
        subset=X(:,S);
        configs=subset;%unique(subset,'rows');
        %sum_lbd=0;
        for c=configs'
            %l=lbd(X,i,S,c');
            %sum_lbd=sum_lbd+l;
            %v=v+l*infl(X,u,i,S,c');
            v = v + infl(X,u,i,S,c')*lbd(X,i,S,c');
        end
        v=v/size(configs,1);
    %end
end