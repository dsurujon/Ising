% Learn Neighborhood of a given node
% [input] u: Node u
% [input] X: mXn matrix of binary observations of m genomes, n genes
% [input] t: threshold value
% [output] S: neighborhood of u

function S = learnNbhd(u,X,t)
    disp(u);
    maxi=1;
    keepadding=1;
    
    Sa=[];
    [w,n]=size(X);
    V=1:n;

    %keep adding to S until the maximum influence is below the threshold
    while keepadding==1
        possible_neighbors=setdiff(V,union(Sa,u));
        maxinfl=0;
        for i=possible_neighbors
            infl=avg_infl(X,u,i,Sa);
            if infl>maxinfl
                maxinfl=infl;
                maxi=i;
            end
        end
        disp(Sa);
        if maxinfl>=t
            Sa=union(Sa,maxi);
        %elseif isequal(Sa,setdiff(V,u))==1
        %    keepadding=0;
        else
            keepadding=0;
        end
        
    end
    %prune S
    for i2=Sa
        infl2=avg_infl(X,u,i2,setdiff(Sa,i2));
        disp(infl2);
        if infl2<t
            Sa=setdiff(Sa,i2);
            disp(Sa);
        end
    end
    S=Sa;
end