%% citrus data

len = 10000;

RandPlant0 = load('ChangingP/Changingp0p0.csv');
RandPlant0 = RandPlant0(1:len);

RandPlant2 = load('ChangingP/Changingp0p2.csv');
RandPlant2 = RandPlant2(1:len);

RandPlant4 = load('ChangingP/Changingp0p4.csv');
RandPlant4 = RandPlant4(1:len);

RandPlant6 = load('ChangingP/Changingp0p6.csv');
RandPlant6 = RandPlant6(1:len);

RandPlant8 = load('ChangingP/Changingp0p8.csv');
RandPlant8 = RandPlant8(1:len);

RandPlant10 = load('ChangingP/Changingp1.csv');
RandPlant10 = RandPlant10(1:len);

majorEpid(1) = length(nonzeros(RandPlant0>50))/len;
majorEpid(2) = length(nonzeros(RandPlant2>50))/len;
majorEpid(3) = length(nonzeros(RandPlant4>50))/len;
majorEpid(4) = length(nonzeros(RandPlant6>50))/len;
majorEpid(5) = length(nonzeros(RandPlant8>50))/len;
majorEpid(6) = length(nonzeros(RandPlant10>50))/len;


figure

a = [0.0,0.2,0.4,0.6,0.8,1.0];
scatter(a, majorEpid,50,'filled')
hold on
plot(a,majorEpid,'--')


xlabel('p values')
ylabel('Probability')
set(gca,'FontSize',36)
ax = gca;
grid on
box on

