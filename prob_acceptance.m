
%prob(x)~exp(sum(theta_{ij}*xi*xj)+sum(theta_i*xi))
function prob = prob_acceptance(neighborhood, edges, i, config,fields,b)
    probsum=config(i)*fields(i);
    probsum=probsum+sum(edges(i,neighborhood).*config(neighborhood))*config(i);
    prob=exp(-2*b*probsum);
end