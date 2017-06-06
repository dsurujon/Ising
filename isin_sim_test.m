%Monte Carlo simulation for the Ising model
rng(1)
n=15;
%start with a small graph with n nodes
nodes=1:n;
%random edges (only show lower half). theta_{ij}=1 if (i,j) is an edge, 0
%otherwise
edges=tril(rand(n)>0.9,-1);
coords=rand(n,2);
field=zeros(1,n)*0.1;
%%
%plot the graph
gplot(edges,coords);
for i=1:n
   text(coords(i,1),coords(i,2),num2str(i)); 
end

%%
curr_config=2*(rand(1,n)>0.5)-1;
num_samples=1000;

all_configs=zeros(num_samples,n);

sampling_times = 1000;
%%
for i=1:num_samples*sampling_times
    %disp(i);
    curr_state=prob_state(curr_config,edges,field);
    %flip nodes one at a time
    newconfig=curr_config;
    flip_ix=randi(n);
    newconfig(flip_ix)=-newconfig(flip_ix);
    %%%only look at neighbours
    next_state=prob_state(newconfig,edges,field);
    acceptance=next_state/curr_state;
    if acceptance>=1
        curr_config=newconfig;
    elseif rand()<acceptance
        curr_config=newconfig;
    end
    if (mod(i,sampling_times)==0)
        disp(i);
        all_configs(i/sampling_times,:)=curr_config;
    end
end
%%
colormap('bone');
imagesc(all_configs');

%%
d=max(sum(edges));alpha=1;beta=1;h=max(field);
delta=0.5*exp(-beta*d-h);
t=(alpha^2)*(delta^(4*d+1))/(16*d*beta);
%%
recovered_edges=learnGraph(all_configs,n,t);

%%
colormap('bone');
subplot(2,2,1:2);imagesc(all_configs');title('Samples generated');
subplot(2,2,3);
gplot(edges,coords);
xlim([0,1]);
ylim([0,1]);
for i=1:n
   text(coords(i,1),coords(i,2),num2str(i)); 
end
title('Real graph');
subplot(2,2,4);
gplot(recovered_edges,coords);
xlim([0,1]);
ylim([0,1]);
for i=1:n
   text(coords(i,1),coords(i,2),num2str(i)); 
end
title('Recovered graph');