%Monte Carlo simulation for the Ising model
r=1;
rng(r)
%define number of nodes and limits of interaction strengths
n=10;
beta=0.1;
alpha=0.1;
h=0;
%inv temperature
b=0.1;
%start with a small graph with n nodes
%and random edges (with weights between alpha and beta)
edg=tril(rand(n)>0.7,-1).*((beta-alpha).*rand(n)+alpha);
ed_sym=sparse(edg+edg.');
g=graph(ed_sym);
neighbors={};
for i=1:n
   neighbors{i}=g.neighbors(i); 
end
%small field
field=ones(1,n).*rand(1,n)*h;

%%
curr_config=2*(rand(1,n)>0.5)-1;
%num=exp(beta*d)*log(n*d/4 - 1)/(4*alpha*d*exp(alpha));
num_samples=5000;

all_configs=zeros(num_samples,n);

sampling_times = 500;
%%
for i=1:num_samples*sampling_times
    %flip nodes one at a time
    newconfig=curr_config;
    flip_ix=randi(n);
    newconfig(flip_ix)=-newconfig(flip_ix);
    %%%only look at neighbours
    acceptance=prob_acceptance(neighbors{flip_ix},ed_sym,flip_ix,curr_config,field,b);
    if acceptance>=1
        curr_config=newconfig;
    elseif rand()<acceptance
        curr_config=newconfig;
    end
    if (mod(i,sampling_times)==0)
        %disp(i);
        all_configs(i/sampling_times,:)=curr_config;
    end
end
%%
colormap('bone');
imagesc(all_configs');

%%
d=max(g.degree)+1;
delta=0.5*exp(-2*(beta*d+h));
t=(alpha^2)*(delta^(4*d+1))/(16*d*beta);
%%
recovered_edges=learnGraph(all_configs,n,t);
g2=graph(recovered_edges);
%%
colormap('bone');
subplot(2,2,1:2);imagesc(all_configs');title('Samples generated');
subplot(2,2,3);
plot(g);
title('Real graph');
subplot(2,2,4);
plot(g2);
title('Recovered graph');


imgname=strcat(['a',num2str(alpha),'b',num2str(beta),'h',num2str(h),'freq',num2str(sampling_times),'n',num2str(num_samples),'rng',num2str(r),'temp',num2str(b)]);
print(imgname,'-dpng');