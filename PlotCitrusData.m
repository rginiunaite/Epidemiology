%% citrus data

len = 100;

Orch2 = load('CitrusData/Citrusa2mu0p1N10000RI100.csv');
Orch2 = Orch2(1:len);

Orch4 = load('CitrusData/Citrusa4mu0p1N10000RI100.csv');
Orch4 = Orch4(1:len);

Orch6 = load('CitrusData/Citrusa20mu0p1N10000RI100.csv');
Orch6 = Orch6(1:len);

Orch8 = load('CitrusData/Citrusa200mu0p1N10000RI100.csv');
Orch8 = Orch8(1:len);

Orch10 = load('CitrusData/Citrusa2000mu0p1N10000RI100.csv');
Orch10 = Orch10(1:len);

majorEpid(1) = length(nonzeros(isinf(Orch2)))/len;
majorEpid(2) = length(nonzeros(isinf(Orch4)))/len;
majorEpid(3) = length(nonzeros(isinf(Orch6)))/len;
majorEpid(4) = length(nonzeros(isinf(Orch8)))/len;
majorEpid(5) = length(nonzeros(isinf(Orch10)))/len;


figure

a = [2,4,6,8,10];
scatter(a, majorEpid,50,'filled')
hold on
plot(a,majorEpid,'--')


xlabel('a values in the exponential kernel')
ylabel('Probability')
set(gca,'FontSize',36)
ax = gca;
grid on
box on

