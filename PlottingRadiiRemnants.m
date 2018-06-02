data=importdata('~/Work/Reviews/RMP/Solar.txt'); %courtesy of Coen

colnames=data.colheaders;
MZAMS=data.data(:,find(strcmp(colnames,'MassZAMS')));
R0sol=data.data(:,find(strcmp(colnames,'minMS')));
RMSsol=data.data(:,find(strcmp(colnames,'MS')));
RHGsol=data.data(:,find(strcmp(colnames,'HG')));
RFGBsol=data.data(:,find(strcmp(colnames,'FGB')));
RCHeBsol=data.data(:,find(strcmp(colnames,'CHeB')));
REAGBsol=data.data(:,find(strcmp(colnames,'EAGB')));
RTPAGBsol=data.data(:,find(strcmp(colnames,'TPAGB')));
RHeMSsol=data.data(:,find(strcmp(colnames,'HeMS')));
RHeHGsol=data.data(:,find(strcmp(colnames,'HeHG')));
RHeWDsol=data.data(:,find(strcmp(colnames,'HeWD')));
RCOWDsol=data.data(:,find(strcmp(colnames,'COWD')));
RONeWDsol=data.data(:,find(strcmp(colnames,'ONeWD')));
TypeCORapidsol=data.data(:,find(strcmp(colnames,'RemnantTypeRapid')));
MCORapidsol=data.data(:,find(strcmp(colnames,'MassRapid')));
TypeCODelayedsol=data.data(:,find(strcmp(colnames,'RemnantTypeDelayed')));
MCODelayedsol=data.data(:,find(strcmp(colnames,'MassDelayed')));

data=importdata('~/Work/Reviews/RMP/TenthSolar.txt'); %courtesy of Coen

colnames=data.colheaders;
MZAMS=data.data(:,find(strcmp(colnames,'MassZAMS')));
R0ten=data.data(:,find(strcmp(colnames,'minMS')));
RMSten=data.data(:,find(strcmp(colnames,'MS')));
RHGten=data.data(:,find(strcmp(colnames,'HG')));
RFGBten=data.data(:,find(strcmp(colnames,'FGB')));
RCHeBten=data.data(:,find(strcmp(colnames,'CHeB')));
REAGBten=data.data(:,find(strcmp(colnames,'EAGB')));
RTPAGBten=data.data(:,find(strcmp(colnames,'TPAGB')));
RHeMSten=data.data(:,find(strcmp(colnames,'HeMS')));
RHeHGten=data.data(:,find(strcmp(colnames,'HeHG')));
RHeWDten=data.data(:,find(strcmp(colnames,'HeWD')));
RCOWDten=data.data(:,find(strcmp(colnames,'COWD')));
RONeWDten=data.data(:,find(strcmp(colnames,'ONeWD')));
TypeCORapidten=data.data(:,find(strcmp(colnames,'RemnantTypeRapid')));
MCORapidten=data.data(:,find(strcmp(colnames,'MassRapid')));
TypeCODelayedten=data.data(:,find(strcmp(colnames,'RemnantTypeDelayed')));
MCODelayedten=data.data(:,find(strcmp(colnames,'MassDelayed')));



figure(1)
set(gca,'FontSize',14), 
semilogy(MZAMS,R0sol, ...
    MZAMS,RMSsol,...
    MZAMS(RHGsol~=0),RHGsol(RHGsol~=0),...
    MZAMS(RCHeBsol~=0),RCHeBsol(RCHeBsol~=0),...
    MZAMS(MZAMS<=23.5),REAGBsol(MZAMS<=23.5),...
    MZAMS(RTPAGBsol~=0),RTPAGBsol(RTPAGBsol~=0),...
    'LineWidth',2), 
axis([0 100 0.8*min(R0sol),1.5*max(RCHeBsol)]),
set(gca,'FontSize',20), 
xlabel('$M / M_\odot$', 'Interpreter','latex'), 
ylabel ('$R / R_\odot$', 'Interpreter','latex'),
legend('ZAMS','MS','HG','CHeB','EAGB','TPAGB','Location','SouthEast');



figure(2)
set(gca,'FontSize',14), 
semilogy(MZAMS,R0ten, ...
    MZAMS,RMSten,...
    MZAMS(RHGten~=0),RHGten(RHGten~=0),...
    MZAMS(RCHeBten~=0),RCHeBten(RCHeBten~=0),...
    MZAMS(REAGBten~=0),REAGBten(REAGBten~=0),...
    MZAMS(RTPAGBten~=0),RTPAGBten(RTPAGBten~=0),...
    'LineWidth',2), 
set(gca,'FontSize',20), 
axis([0 200 0.9*min(R0ten),1.1*max(RCHeBten)]),
xlabel('$M / M_\odot$', 'Interpreter','latex'), 
ylabel ('$R / R_\odot$', 'Interpreter','latex'),
legend('ZAMS','MS','HG','CHeB','EAGB','TPAGB','Location','SouthEast');

figure(3)
set(gca,'FontSize',14),
plot(MZAMS,MCORapidsol,...
	MZAMS,MCORapidten, 'LineWidth',2),
set(gca,'FontSize',20), 
xlabel('$M_\textrm{ZAMS} / M_\odot$', 'Interpreter','latex'), 
ylabel ('$M_\textrm{CO} / M_\odot$', 'Interpreter','latex'),
legend('$Z=Z_\odot$','$Z=0.1\,Z_\odot$','Interpreter','latex','Location','NorthWest');
legend('solar Z','0.1 solar Z','Location','NorthWest');

figure(4)
set(gca,'FontSize',14),
loglog(MZAMS,MCODelayedsol,...
	MZAMS,MCODelayedten, 'LineWidth',2),
axis([0.8*min(MZAMS(MCODelayedten~=0)) 1.2*max(MZAMS) 1.0 1.2*max(MCODelayedten)])
set(gca,'FontSize',20), 
xlabel('$M_\textrm{ZAMS} / M_\odot$', 'Interpreter','latex'), 
ylabel ('$M_\textrm{CO} / M_\odot$', 'Interpreter','latex'),
%legend('$Z=Z_\odot$','$Z=0.1\,Z_\odot$','Interpreter','latex','Location','NorthWest');
legend('solar Z','0.1 solar Z','Location','NorthWest');

%min(MCODelayedsol(TypeCODelayedsol==14))
%max(MCODelayedsol(TypeCODelayedsol==13))

figure(5)
set(gca,'FontSize',14),
plot(MZAMS,MCORapidsol,...
    MZAMS,MCORapidten,...
    MZAMS, MCODelayedsol,...
    MZAMS, MCODelayedten, 'LineWidth',2),
set(gca,'FontSize',20), 
xlabel('$M_\textrm{ZAMS} / M_\odot$', 'Interpreter','latex'), 
ylabel ('$M_\textrm{CO} / M_\odot$', 'Interpreter','latex'),
l=legend('$Z=Z_\odot$, Rapid SN','$Z=0.1\,Z_\odot$, Rapid SN','$Z=Z_\odot$, Delayed SN','$Z=0.1\,Z_\odot$, Delayed SN',...
    'Location','NorthWest');
set(l,'Interpreter','latex');


RLOF=0.49/(0.6+log(1+1));
amaxinit=(r0HubbleE0/Rsun)*(MCODelayedsol/30).^(3/4);  %from Peters file
RLOFinit=RLOF*amaxinit;
figure(6)
set(gca,'FontSize',14), 
loglog(MZAMS,R0sol, ...
    MZAMS,RMSsol,...
    MZAMS(RHGsol~=0),RHGsol(RHGsol~=0),...
    MZAMS(RCHeBsol~=0),RCHeBsol(RCHeBsol~=0),...
    MZAMS(MZAMS<=23.5),REAGBsol(MZAMS<=23.5),...
    MZAMS, RLOFinit,'k-.',...
    'LineWidth',2), 
axis([min(MZAMS(MCORapidsol~=0)) 100 0.8*min(RLOFinit(MCORapidsol~=0)),1.5*max(RCHeBsol)]),
set(gca,'FontSize',20), 
xlabel('$M / M_\odot$', 'Interpreter','latex'), 
ylabel ('$R / R_\odot$', 'Interpreter','latex'),
legend('ZAMS','MS','HG','CHeB','AGB','Roche Lobe','Location','SouthEastOutside');
