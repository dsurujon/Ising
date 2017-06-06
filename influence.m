% Calculate empirical average influence
% [input] u: Node u
% [input] i: Node i
% [input] X: mXn matrix of binary observations of m genomes, n genes
% [input] S: neighborhood of u
% [output] vavg: average empirical influence of i on u

function vavg = influence(u,i,S,X)
    vavg=0;
    if(ismember(i,S)==1)
        return;
    elseif (isempty(S)==1)
         Xsiu=X(:,[i,u]);
         %a: identities of unique rows
         %b: indices of (irst occurrence of) rows of a
         %c: vector matching rows of a to Xsiu
        [xsiua,xsiub,xsiuc]=unique(Xsiu,'rows');
        %get the counts of each unique row in Xsiu
        countsxsiu=histc(xsiuc,unique(xsiuc));
        
        %[T/F, index] of row in xaiua
        [q11,w11]=ismember([1,1],xsiua,'rows');
        [q10,w10]=ismember([1,-1],xsiua,'rows');
        [q01,w01]=ismember([-1,1],xsiua,'rows');
        [q00,w00]=ismember([-1,-1],xsiua,'rows');
        if w11>0
            i1u1count=countsxsiu(w11);
        else
            i1u1count=0;
        end

        if w10>0
            i1u0count=countsxsiu(w10);
        else
            i1u0count=0;
        end

        if w01>0
            i0u1count=countsxsiu(w01);
        else
            i0u1count=0;
        end
        if w00>0
            i0u0count=countsxsiu(w00);
        else
            i0u0count=0;
        end
        xs1count=i1u0count+i1u1count;
        xs0count=i0u0count+i0u1count;
        vavg=(i1u1count/xs1count)-(i0u1count/xs0count);

    else
        Xs=X(:,S);
        Xsiu=X(:,[S,i,u]);
        [xsa,xsb,xsc]=unique(Xs,'rows');
        countsxs=histc(xsc,unique(xsc));
        [xsiua,xsiub,xsiuc]=unique(Xsiu,'rows');
        countsxsiu=histc(xsiuc,unique(xsiuc));
        %iterate over all possible xs
        for xs_ix=1:length(xsa)
            %get the xs row
            xsrow=xsa(xs_ix,:);
            %and the #times that row appears in Xs
            xscount=countsxs(xs_ix);
            
            %rows with i=1/-1, u=1/-1
            xsrowi1u1=[xsrow,1,1];
            xsrowi1u0=[xsrow,1,-1];
            xsrowi0u1=[xsrow,-1,1];
            xsrowi0u0=[xsrow,-1,-1];
            %find these rows in the uniques set (xsiua)
            [q11,w11]=ismember(xsrowi1u1,xsiua,'rows');
            [q10,w10]=ismember(xsrowi1u0,xsiua,'rows');
            [q01,w01]=ismember(xsrowi0u1,xsiua,'rows');
            [q00,w00]=ismember(xsrowi0u0,xsiua,'rows');
            %disp([w11,w10,w01,w00]);
            
            %and get the #times they appear in Xsiu
            if w11>0
                i1u1count=countsxsiu(w11);
            else
                i1u1count=0;
            end
            
            if w10>0
                i1u0count=countsxsiu(w10);
            else
                i1u0count=0;
            end
            
            if w01>0
                i0u1count=countsxsiu(w01);
            else
                i0u1count=0;
            end
            if w00>0
                i0u0count=countsxsiu(w00);
            else
                i0u0count=0;
            end
            xs1count=i1u0count+i1u1count;
            xs0count=i0u0count+i0u1count;
            lbd=2*xs1count*xs0count/(xscount^2); 
            vuixs=(i1u1count/xs1count)-(i0u1count/xs0count);
            vavg=vavg+lbd*abs(vuixs);
            
        end

    end
end