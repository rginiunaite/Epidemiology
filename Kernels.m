%% Plotting exponential kernels

d = [0:0.01:100];
a = 0.02;

figure

for i=1:5
    hold on
    a = 1*i;
    Kdist = exp(-d/a) * 1/(2*pi*a*a);
    plot(d,Kdist,'LineWidth',2)

end

%  ylim([0,50])
%  xlim([0,6])
legend(['a=0.02';'a=0.04';'a=0.06';'a=0.08';'a=0.10'])
xlabel('Distance')
ylabel('Exponential kernel')
ax = gca;  
set(gca,'FontSize',36)
set(gca,'linew',2)
grid on
box on