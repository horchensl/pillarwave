
%% S11

figure(2)
clf
colororder('default')
plot(time(1:Nsamples),innerDiameter(:,2:end)*1e3)
xlabel('time [s]')
ylabel('inner diameter [mm]')
legend(labels(2:end))
grid

%%

pd = 74.4; % diastolic pressure
diastolicDiameter = 6.35e-3;
vesselRigidity = 6.1;

innerPressure = pd.*exp(vesselRigidity.*(innerDiameter.^2./diastolicDiameter.^2 - 1));

lowestDiameter=6.3e-3;
highestDiameter=6.71e-3;

lowestPressure = pd.*exp(vesselRigidity.*(lowestDiameter.^2./diastolicDiameter.^2 - 1));
highestPressure = pd.*exp(vesselRigidity.*(highestDiameter.^2./diastolicDiameter.^2 - 1));

diameterAxis = lowestDiameter:1e-4:highestDiameter;
pressureAxis = pd.*exp(vesselRigidity.*(diameterAxis.^2./diastolicDiameter.^2 - 1));

%%

ind = 5;
referenceStartIndex = 253;

referenceTimeAxis = times{5}-times{5}(1);

h=gobjects(2,1);


%%

referencePressure = pressures{ind}(((1:Nsamples)+referenceStartIndex)*4-2); % pressure measured for reference


%% S12


figure(3)
clf
ax = axes();
set(gcf,'Position',[100 100 1200 400])
h(1)=plot(referenceTimeAxis((1:Nsamples)*4-3),referencePressure,'k','LineWidth',4,'DisplayName','reference');
hold on
h(2)=plot(time(1:Nsamples),innerPressure(:,ind),'r.-','LineWidth',2,'MarkerSize',15,'DisplayName','derived pressure');
ylim([70 140]);
lim = ylim;

hold off
ylabel('pressure [mmHg]')
yyaxis(ax, 'right')
ylim(lim)

ylabel('diameter [mm]')
ax.YTick = pressureAxis;
ax.YTickLabel = diameterAxis*1e3;

xlabel('time [s]')

legend(h([2 1]),'location','east')
grid on
set(gca,'FontSize',20)
set(gca,'LineWidth',2)
colororder({'k','k'})
