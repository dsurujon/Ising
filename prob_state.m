
%prob(x)~exp(sum(theta_{ij}*xi*xj)+sum(theta_i*xi))
function prob = prob_state(config,edges,field,n)
    probsum=0;
    for i=1:n
        for j=1:n
            probsum=probsum+edges(i,j)*config(i)*config(j);
        end
        probsum=probsum+field(i)*config(i);
    end
    prob=exp(probsum);
end