%% Stochastic SIR
clear all;
%run 10000 simulations
Nsim = 10000;
finalsizeepid = zeros(1,Nsim);

for j=1:Nsim

    %Initialise the number of individuals in each of the S, I and R classes in the model, and set the outbreak time t = 0.

    S(1) = 999;
    N = 1000;
    I(1) = 1;
    R(1) = 0;
    beta = 0.0003;%15;
    mu = 0.1;

    R0 = beta*N/mu;

    t(1) = 0;

    n=1;

    %Steps 2-4 should be repeated while the outbreak is still ongoing (i.e. I > 0).  First calculate two random numbers r1, r2 each uniformly distributed in (0,1).
     while I(n) > 0
        rand1 = rand(1,1);
        rand2 = rand(1,1);

        t(n+1) = t(n) + 1/ (beta*I(n)*S(n) + mu * I(n)) * log(1/rand1);

        if (rand2 < beta*I(n)*S(n)/(beta*I(n)*S(n)+mu*I(n)))
            S(n+1) = S(n)-1;
            I(n+1) = I(n)+1;
            R(n+1) = R(n);
        else
            I(n+1) = I(n)-1;
            R(n+1) = R(n)+1;
            S(n+1) = S(n);
        end
        n = n+1;
     end

     finalsizeepid(j) = R(n);
     clear I;
     clear S;
     clear R;
end 




[N,edges] = histcounts(finalsizeepid,50, 'Normalization', 'probability');
xbar = edges(1:numel(N)) + mean(diff(edges))/2;
figure
bar(xbar, N)
% grid
% yt = get(gca, 'YTick'); 
% ytix = linspace(min(yt), max(yt), 10);
% set(gca, 'YTick',ytix, 'YTickLabel',fix(ytix*200/max(yt)))
% 


xlabel('Number ever infected')
ylabel('Probability')
set(gca,'FontSize',36)
ax = gca;
grid on
box on

ProbMajorEpid = 1-N(1);
