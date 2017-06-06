function newedges=learnGraph(X,n,t)
    newedges=zeros(n);
    for i=1:n
        n_i=learnNbhd(i,X,t);
        newedges(i,n_i)=1;
    end
    re_sym=newedges==newedges';
    re_final(re_sym)=newedges(re_sym);
    newedges=reshape(re_final,n,n);
end
