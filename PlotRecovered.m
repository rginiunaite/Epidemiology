% plot recovered


%recovered = load('Recoveredmu0p1beta0p0003Kd1minusdsquared.csv');
recovered = load('NewKd1beta0p0003mu0p1.csv');


[N,edges] = histcounts(recovered,50, 'Normalization', 'probability');
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
%ylim([0,0.7])
ax = gca;
grid on
box on

ProbMajorEpid = 1-N(1);
