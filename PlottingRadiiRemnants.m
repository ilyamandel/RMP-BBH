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
	MZAMS,MCODelayedten, 'LineWidth',3),
axis([0.8*min(MZAMS(MCODelayedten~=0)) 1.2*max(MZAMS) 1.0 1.2*max(MCODelayedten)])
set(gca,'FontSize',22), grid on;
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
set(gca,'FontSize',18), 
loglog(MZAMS,R0sol, ...
    MZAMS,RMSsol,...
    MZAMS(RHGsol~=0),RHGsol(RHGsol~=0),...
    MZAMS(RCHeBsol~=0),RCHeBsol(RCHeBsol~=0),...
    MZAMS(MZAMS<=23.5),REAGBsol(MZAMS<=23.5),...
    MZAMS, RLOFinit,'k-.',...
    'LineWidth',4), 
axis([min(MZAMS(MCORapidsol~=0)) 100 0.8*min(RLOFinit(MCORapidsol~=0)),1.5*max(RCHeBsol)]),
set(gca,'FontSize',26), grid on; 
xlabel('{\bf $M / M_\odot$}', 'Interpreter','latex','FontSize',30), 
ylabel ('$R / R_\odot$', 'Interpreter','latex', 'FontSize',30),
legend('ZAMS','MS','HG','CHeB','AGB','Roche Lobe','Location','East');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%recreate grid file -- BSE
%circular binaries with a separation so wide that they are effectively single
m1=1.5*10.^[0:0.01:2]; N=length(m1);
file=fopen('/Users/ilyam/Work/COMPAS/COMPAS/src/BSEgridMandelFarmer.txt','w');
for(i=1:N),
    gridline=['--initial-mass-1 ', num2str(m1(i)), ' --initial-mass-2 0.1 --metallicity 0.0142 --eccentricity 0.0 --semi-major-axis 1e+20\n'];
    fprintf(file,gridline);
end;
for(i=1:N),
    gridline=['--initial-mass-1 ', num2str(m1(i)), ' --initial-mass-2 0.1 --metallicity 0.001 --eccentricity 0.0 --semi-major-axis 1e+20\n'];
    fprintf(file,gridline);
end;
fclose(file);

%recreate grid file -- SSE
file=fopen('/Users/ilyam/Work/COMPAS/COMPAS/src/SSEgridMandelFarmer.txt','w');
for(i=1:N),
    gridline=['--initial-mass ', num2str(m1(i)), ' --metallicity 0.0142\n'];
    fprintf(file,gridline);
end;
for(i=1:N),
    gridline=['--initial-mass ', num2str(m1(i)), ' --metallicity 0.0001\n'];
    fprintf(file,gridline);
end;   
fclose(file);

%Redo figures 3,4 of Mandel & Farmer using new BSE file
M=zeros(N,1); R0=zeros(N,1); RMS=zeros(N,1); RHG=zeros(N,1); RHeB=zeros(N,1); RAGB=zeros(N,1); MCO=zeros(N,1);
%solar, Z=0.0142
prefix='/Users/ilyam/Work/COMPAS/COMPAS/src/BSEgridMandelFarmer/Detailed_Output/BSE_Detailed_Output_';
for(i=1:N),
    file=[prefix,num2str(i-1),'.h5'];
    mass=h5read(file,'/Mass(1)');
    radius=h5read(file,'/Radius(1)');
    time=h5read(file,'/Time');
    ST=h5read(file,'/Stellar_Type(1)');
    M(i)=mass(1);
    R0(i)=radius(1);
    if(~isempty(max(radius(ST==1 | ST==0)))),  RMS(i)=max(radius(ST==1 | ST==0)); end;
    if(~isempty(max(radius(ST==2)))), RHG(i)=max(radius(ST==2)); end;
    if(~isempty(max(radius(ST==3 | ST==4)))), RHeB(i)=max(radius(ST==3 | ST==4)); end;
    if(~isempty(max(radius(ST==5 | ST==6)))), RAGB(i)=max(radius(ST==5 | ST==6)); end;
    if(ST(length(mass))==13 | ST(length(mass))==14), MCO(i)=mass(length(mass)); end;
end;
%Z=0.001
Mlow=zeros(N,1); MCOlow=zeros(N,1);
for(i=1:N),
    file=[prefix,num2str(i+N-1),'.h5'];
    mass=h5read(file,'/Mass(1)');
    ST=h5read(file,'/Stellar_Type(1)');
    Mlow(i)=mass(1);
    if(ST(length(mass))==13 | ST(length(mass))==14), MCOlow(i)=mass(length(mass)); end;
end;


%Redo figures 3,4 of Mandel & Farmer using SSE file
M=zeros(N,1); R0=zeros(N,1); RMS=zeros(N,1); RHG=zeros(N,1); RHeB=zeros(N,1); RAGB=zeros(N,1); MCO=zeros(N,1);
%solar, Z=0.0142
prefix='/Users/ilyam/Work/COMPAS/COMPAS/src/SSEgridMandelFarmerSmallStepZ0001/Detailed_Output/SSE_Detailed_Output_';
for(i=1:N),
    file=[prefix,num2str(i-1),'.h5'];
    mass=h5read(file,'/Mass');
    radius=h5read(file,'/Radius');
    time=h5read(file,'/Time');
    ST=h5read(file,'/Stellar_Type');
    dT=h5read(file,'/dT');
    Hemass=h5read(file,'/Mass_He_Core');
    M(i)=mass(1);
    MHe(i)=Hemass(length(mass)-1);
    R0(i)=radius(1);
    dTmin(i)=min(dT);
    if(~isempty(max(radius(ST==1 | ST==0)))),  RMS(i)=max(radius(ST==1 | ST==0)); end;
    if(~isempty(max(radius(ST==2)))), RHG(i)=max(radius(ST==2)); end;
    if(~isempty(max(radius(ST==3 | ST==4)))), RHeB(i)=max(radius(ST==3 | ST==4)); end;
    if(~isempty(max(radius(ST==5 | ST==6)))), RAGB(i)=max(radius(ST==5 | ST==6)); end;
    if(ST(length(mass))==13 | ST(length(mass))==14), MCO(i)=mass(length(mass)); end;
end;
%Z=0.001
Mlow=zeros(N,1); MCOlow=zeros(N,1);
for(i=1:N),
    file=[prefix,num2str(i+N-1),'.h5'];
    mass=h5read(file,'/Mass');
    Hemass=h5read(file,'/Mass_He_Core');
    ST=h5read(file,'/Stellar_Type');
    Mlow(i)=mass(1);
    MHelow(i)=Hemass(length(mass)-1);
    if(ST(length(mass))==13 | ST(length(mass))==14), MCOlow(i)=mass(length(mass)); end;
end;

%Redo figures 3,4 of Mandel & Farmer using old BSE file
M=zeros(500,1); R0=zeros(500,1); RMS=zeros(500,1); RHG=zeros(500,1); RHeB=zeros(500,1); RAGB=zeros(500,1); MCO=zeros(500,1);
%solar, Z=0.0142
prefix='/Users/ilyam/Work/COMPAS/COMPAS/examples/methods_paper_plots/fig_6_max_R/COMPAS_Output/Detailed_Output/BSE_Detailed_Output_';
for(i=1:500),
    file=[prefix,num2str(i-1),'.h5'];
    mass=h5read(file,'/Mass(1)');
    radius=h5read(file,'/Radius(1)');
    time=h5read(file,'/Time');
    ST=h5read(file,'/Stellar_Type(1)');
    M(i)=mass(1);
    R0(i)=radius(1);
    if(~isempty(max(radius(ST==1 | ST==0)))),  RMS(i)=max(radius(ST==1 | ST==0)); end;
    if(~isempty(max(radius(ST==2)))), RHG(i)=max(radius(ST==2)); end;
    if(~isempty(max(radius(ST==3 | ST==4)))), RHeB(i)=max(radius(ST==3 | ST==4)); end;
    if(~isempty(max(radius(ST==5 | ST==6)))), RAGB(i)=max(radius(ST==5 | ST==6)); end;
    if(ST(length(mass))==13 | ST(length(mass))==14), MCO(i)=mass(length(mass)); end;
end;
%Z=0.001
Mlow=zeros(500,1); MCOlow=zeros(500,1);
for(i=1:500),
    file=[prefix,num2str(i+499),'.h5'];
    mass=h5read(file,'/Mass(1)');
    ST=h5read(file,'/Stellar_Type(1)');
    Mlow(i)=mass(1);
    if(ST(length(mass))==13 | ST(length(mass))==14), MCOlow(i)=mass(length(mass)); end;
end;
    
figure(111)
set(gca,'FontSize',14), 
semilogy(M,R0, ...
    M,RMS,...
    M(RHG~=0),RHG(RHG~=0),...
    M(RHeB~=0),RHeB(RHeB~=0),...
    M(RAGB~=0),RAGB(RAGB~=0),...
    'LineWidth',2), 
axis([0 100 0.8*min(R0sol),1.5*max(RCHeBsol)]),
set(gca,'FontSize',20), 
xlabel('$M / M_\odot$', 'Interpreter','latex'), 
ylabel ('$R / R_\odot$', 'Interpreter','latex'),
legend('ZAMS','MS','HG','CHeB','AGB','Location','SouthEast');
    
figure(142)
set(gca,'FontSize',14),
loglog(M,MCO,'*',...
	Mlow,MCOlow,'*', 'LineWidth',3),
axis([0.8*min(Mlow(MCOlow~=0)) 1.2*max(M) 1.0 1.2*max(MCOlow)])
set(gca,'FontSize',22), grid on;
xlabel('$M_\textrm{ZAMS} / M_\odot$', 'Interpreter','latex'), 
ylabel ('$M_\textrm{CO} / M_\odot$', 'Interpreter','latex'),
%legend('$Z=Z_\odot$','$Z=0.1\,Z_\odot$','Interpreter','latex','Location','NorthWest');
legend('Z=0.0142','Z=0.0001','Location','NorthWest');

figure(143)
set(gca,'FontSize',14),
loglog(M,MHe,'*',...
	Mlow,MHelow,'*', 'LineWidth',3),
axis([0.8*min(Mlow(MHelow~=0)) 1.2*max(M) 1.0 1.2*max(MHelow)])
set(gca,'FontSize',22), grid on;
xlabel('$M_\textrm{ZAMS} / M_\odot$', 'Interpreter','latex'), 
ylabel ('$M_\textrm{He core} / M_\odot$', 'Interpreter','latex'),
%legend('$Z=Z_\odot$','$Z=0.1\,Z_\odot$','Interpreter','latex','Location','NorthWest');
legend('Z=0.0142','Z=0.0001','Location','NorthWest');

figure(160)
RLOF=0.49/(0.6+log(2));
Msunkg=1.98892e30;	c=299792458;		G=6.67428e-11;		Rsun = 695500000; 
T0=14*1e9*3.15e7;
beta=64/5*G^3*MCO.^2.*(MCO+MCO)*Msunkg^3/c^5;
amaxinit=(T0*4*beta).^0.25;
RLOFinit=RLOF*amaxinit/Rsun;
set(gca,'FontSize',18), 
loglog(M,R0, ...
    M,RMS,...
    M(RHG~=0),RHG(RHG~=0),...
    M(RHeB~=0),RHeB(RHeB~=0),...
    M(RAGB~=0 & RAGB>RHeB),RAGB(RAGB~=0 & RAGB>RHeB),...
    M(MCO~=0), RLOFinit(MCO~=0),'k-.',...
    'LineWidth',4), 
axis([min(M(MCO~=0)) 100 0.8*min(RLOFinit(MCO~=0)),1.5*max(max(RHeB),max(RAGB))]),
set(gca,'FontSize',26), grid on; 
xlabel('{\bf $M / M_\odot$}', 'Interpreter','latex','FontSize',30), 
ylabel ('$R / R_\odot$', 'Interpreter','latex', 'FontSize',30),
legend('ZAMS','MS','HG','HeB','AGB','Roche Lobe','Location','East');