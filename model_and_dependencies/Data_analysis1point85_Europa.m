%Calculating volume at each time
    set(groot,'defaultAxesFontName','Times')
    set(groot,'defaultAxesFontSize',20)
    set(groot,'defaulttextinterpreter','latex')
    set(groot,'defaultAxesTickLabelInterpreter','latex')
    set(groot,'defaultLegendInterpreter','latex')
    set(groot, 'DefaultFigureVisible', 'on');

blue = [57 106 177]./255;
red = [204 37 41]./255;
black = [83 81 84]./255;
green = [62 150 81]./255;
brown = [146 36 40]./255;
purple = [107 76 154]./255;

d = 10e3;   %Thickness of ice shell [in m]   




load("../Output/Europa03321_eta0_14kc1.85e-08_Ea_50_output_1C.mat"); %loading file
interface = (find(Grid.p.yc>0)); interface = interface(1,1); %First cell of no ocean
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean
VolumeDS = [sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)];
Time_arrDS = [tVec];
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+10);  %Just staying away from ocean
pene_depthDS= [Grid.p.yc(pene_ind(1))];

i = 10;

%kc1.85e-8m2
while i<=22200
%Individual file
load(sprintf("../Output/Europa03321_eta0_14kc1.85e-08_Ea_50_output_%sC.mat",num2str(i))); %loading file
Time_arrDS = [Time_arrDS;tVec(i)]; %Time array [years]
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean

Vol = sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3);

if VolumeDS(end)<Vol
    Vol = VolumeDS(end); 
end 

VolumeDS = [VolumeDS;Vol]; %Calculating volume [in m^3]

%Penetration depth
%pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
%pene_ind(pene_ind>interface+10);  %Just staying away from ocean
%pene_depthDS = [pene_depthDS;Grid.p.yc(pene_ind(1))];

[Xc,Yc] = meshgrid(Grid.p.xc,Grid.p.yc); Yc = Yc(:);
Yc_wet = Yc(find(phi>1e-5)); pene_depthdummy = min(Yc_wet(Yc_wet>0));
%pene_depthS = [pene_depthS;Grid.p.yc(pene_ind(1))];
pene_depthDS = [pene_depthDS;pene_depthdummy];

if i<200
    i = i+10;
else
    i = i+100;
end

end

%[Xc,Yc] = meshgrid(Grid.p.xc,Grid.p.yc); Yc = Yc(:);
%Yc_wet = Yc(find(phi>1e-4)); min(Yc_wet(Yc_wet>0))

%Stokes
load("../Output/Europa03321_eta0_14kc1.85e-16_Ea_50_output_1C.mat"); %loading file
interface = (find(Grid.p.yc>0)); interface = interface(1,1); %First cell of no ocean
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean
VolumeS = [sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)];
Time_arrS = [tVec];
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+10);  %Just staying away from ocean
pene_depthS= [Grid.p.yc(pene_ind(1))];
i = 100;
while i<=36000
%Individual file
load(sprintf("../Output/Europa03321_eta0_14kc1.85e-16_Ea_50_output_%sC.mat",num2str(i))); %loading file

Time_arrS = [Time_arrS;tVec(i)]; %Time array [years]
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean

%Calculating volume
VolumeS = [VolumeS;sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)]; %Calculating volume [in m^3]

%Penetration depth
%pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
%pene_ind(pene_ind>interface+1);  %Just staying away from ocean
[Xc,Yc] = meshgrid(Grid.p.xc,Grid.p.yc); Yc = Yc(:);
Yc_wet = Yc(find(phi>1e-5)); pene_depthdummy = min(Yc_wet(Yc_wet>0));
%pene_depthS = [pene_depthS;Grid.p.yc(pene_ind(1))];
pene_depthS = [pene_depthS;pene_depthdummy];

i = i+100;
end

%Mid
load("../Output/Europa03321_eta0_14kc1.85e-10_Ea_50_output_1C.mat"); %loading file
interface = (find(Grid.p.yc>0)); interface = interface(1,1); %First cell of no ocean
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean
VolumeM = [sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)];
Time_arrM = [tVec(1)]
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+10);  %Just staying away from ocean
pene_depthM= [Grid.p.yc(pene_ind(1))];
i = 2;
while i<=13600
%Individual file
load(sprintf("../Output/Europa03321_eta0_14kc1.85e-10_Ea_50_output_%sC.mat",num2str(i))); %loading file

Time_arrM = [Time_arrM;tVec(i)]; %Time array [years]
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean

%Calculating volume
Vol = sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3);

if VolumeM(end)<Vol
    Vol = VolumeM(end); 
end

VolumeM = [VolumeM;Vol]; %Calculating volume [in m^3]

%Penetration depth
%pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
%pene_ind(pene_ind>interface+1);  %Just staying away from ocean
%pene_depthM = [pene_depthM;Grid.p.yc(pene_ind(1))];

[Xc,Yc] = meshgrid(Grid.p.xc,Grid.p.yc); Yc = Yc(:);
Yc_wet = Yc(find(phi>1e-5)); pene_depthdummy = min(Yc_wet(Yc_wet>0));
%pene_depthS = [pene_depthS;Grid.p.yc(pene_ind(1))];
pene_depthM = [pene_depthM;pene_depthdummy];

if i<400
    i = i+1;
else
    i = i+100;
end

end

%Plotting
%Volume plot
figure()
semilogx(Time_arrDS,VolumeDS/1e9,'-','Linewidth',5,'color',red)
hold on
semilogx(Time_arrS,VolumeS/1e9,'-','Linewidth',5,'color',blue)
xlabel('Time [years]');
ylabel('Melt volume [km$$^3$$]');
legend('Darcy-Stokes','Stokes');

%Volume plot
figure()
semilogx(Time_arrDS,pene_depthDS/1e3,'-','Linewidth',5,'color',red)
hold on
semilogx(Time_arrS,pene_depthS/1e3,'-','Linewidth',5,'color',blue)
xlabel('Time [years]');
ylabel('Penetration depth [km]');
legend('Darcy-Stokes','Stokes');


%Combined
hh=figure
hh.Units = 'centimeters'
hh.Position = [1,1,20,15];
set(gcf,'visible','on')
yyaxis left
semilogx(Time_arrDS,VolumeDS/1e9,'k-','Linewidth',3)
hold on
semilogx(Time_arrM,VolumeM/1e9,'k-.','Linewidth',3)
semilogx(Time_arrS,VolumeS/1e9,'k--','Linewidth',3)
semilogx(Time_arrDS,VolumeDS/1e9,'-','Linewidth',3,'color',blue)
semilogx(Time_arrM,VolumeM/1e9,'-.','Linewidth',3,'color',blue)
hold on
semilogx(Time_arrS,VolumeS/1e9,'--','Linewidth',3,'color',blue)
xlabel('Time [years]');
ylabel('Melt volume [km$$^3$$]','Color',blue);
ax = gca
ax.YColor = blue;
yyaxis right
semilogx(Time_arrDS,10 - pene_depthDS*d/1e3,'-','Linewidth',3,'color',red)
hold on
semilogx(Time_arrM,10 - pene_depthM*d/1e3,'-.','Linewidth',3,'color',red)
semilogx(Time_arrS,10 - pene_depthS*d/1e3,'--','Linewidth',3,'color',red)
xlim([1e-4,1e5])
xticks([1e-4,1e-2,1,1e2,1e4])
set(gca, 'YDir','reverse')
ax = gca
ax.YColor = red;
ylabel('Penetration depth [km]','Color',red);
%lh = legend('$$\textrm{k}_0=1.85\cdot10^{-8}\textrm{m}^2$$','$$\textrm{k}_0=1.85\cdot 10^{-10}\textrm{m}^2$$','$$\textrm{k}_0=1.85\cdot 10^{-16}\textrm{m}^2$$','location','southwest');
%set( lh, 'Box', 'off' ) ;
saveas(hh,sprintf('../figures/CombinedVolandpene.png')); 
saveas(hh,sprintf('../figures/CombinedVolandpene.pdf')); 


%Combined Stokes
hh=figure()
yyaxis left
semilogx(Time_arrS,VolumeS/1e9,'-','Linewidth',5,'color',blue)
hold on
xlabel('Time [years]');
ylabel('Melt volume [km$$^3$$]','Color',blue);
ax = gca
ax.YColor = blue;
yyaxis right
semilogx(Time_arrS,10 - pene_depthS*d/1e3,'-','Linewidth',5,'color',red)
hold on
%ylim([0,45])
set(gca, 'YDir','reverse')
ax = gca
ax.YColor = red;
ylabel('Penetration depth [km]','Color',red);
saveas(hh,sprintf('../figures/CombinedStokes.png')); 
saveas(hh,sprintf('../figures/CombinedStokes.pdf')); 

%stop

%Combined Darcy-Stokes
hh=figure()
yyaxis left
semilogx(Time_arrDS,VolumeDS/1e9,'-','Linewidth',5,'color',blue)
hold on
xlabel('Time [years]');
ylabel('Melt volume [km$$^3$$]');
ax = gca
ax.YColor = blue;
yyaxis right
semilogx(Time_arrDS,10 - pene_depthDS*d/1e3,'-','Linewidth',5,'color',red)
hold on
set(gca, 'YDir','reverse')
ax = gca
ax.YColor = red;
ylabel('Penetration depth [km]');
saveas(hh,sprintf('../figures/CombinedDarcyStokes.png')); 
saveas(hh,sprintf('../figures/CombinedDarcyStokes.pdf')); 


%Combined Vol
hh=figure()
semilogx(Time_arrDS,VolumeDS/1e9,'-','Linewidth',5,'color',blue)
hold on
semilogx(Time_arrS,VolumeS/1e9,'--','Linewidth',5,'color',blue)
xlabel('Time [years]');
ylabel('Melt volume [km$$^3$$]');
legend('$$\textrm{k}_0=5.6\times10^{-11}\textrm{m}^2$$','$$\textrm{k}_0=5.6\times10^{-16}\textrm{m}^2$$','location','southwest');
saveas(hh,sprintf('../figures/CombinedVol.png')); 
saveas(hh,sprintf('../figures/CombinedVol.pdf')); 



%%Combined all k0

%Testing purposes only
d = 10e3
load("../Output/Europa03321_eta0_14kc1.85e-08_Ea_50_output_1C.mat"); %loading file
interface = (find(Grid.p.yc>0)); interface = interface(1,1); %First cell of no ocean
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean
Volumetest = [sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)];
Time_arrtest = [tVec];
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+10);  %Just staying away from ocean
pene_depthtest= [Grid.p.yc(pene_ind(1))];
i = 100;

%kc1.85e-8m2
while i<=22200
%Individual file
load(sprintf("../Output/Europa03321_eta0_14kc1.85e-08_Ea_50_output_%sC.mat",num2str(i))); %loading file
Time_arrtest = [Time_arrtest;tVec(i)]; %Time array [years]
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean

Volumetest = [Volumetest;sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)]; %Calculating volume [in m^3]

%Penetration depth
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+10);  %Just staying away from ocean
pene_depthtest = [pene_depthtest;Grid.p.yc(pene_ind(1))];

i = i+100;
end

hh=figure()
yline(phiOrig/1e9,'-','Linewidth',5,'color',red);
hold on
xlim([0, max(tVec)])
plot(Time_arrtest,Volumetest/1e9,'-','Linewidth',5,'color',blue);
plot(tVec,phiDrain1Vec/1e9,'-','Linewidth',5,'color',green);
xlabel('Time [years]');
ylabel('Melt volume [km$$^3$$]');
legend('Initial total','Left','Delivered','location','northeast');
saveas(hh,sprintf('../figures/1.85e-8_vol.png')); 
saveas(hh,sprintf('../figures/1.85e-8_temperature_vol.pdf'));

drain_V_array = [phiDrain1Vec(end)/1e9]; %in km3
kc_array    = [1.85e-8]; %in m2

%Testing purposes only kc=1.85e-10 m2
d = 10e3
load("../Output/Europa03321_eta0_14kc1.85e-10_Ea_50_output_1C.mat"); %loading file
interface = (find(Grid.p.yc>0)); interface = interface(1,1); %First cell of no ocean
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean
Volumetest = [sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)];
Time_arrtest = [tVec];
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+10);  %Just staying away from ocean
pene_depthtest= [Grid.p.yc(pene_ind(1))];
i = 100;

%kc1.85e-10m2
while i<=13600
%Individual file
load(sprintf("../Output/Europa03321_eta0_14kc1.85e-10_Ea_50_output_%sC.mat",num2str(i))); %loading file
Time_arrtest = [Time_arrtest;tVec(i)]; %Time array [years]
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean

Volumetest = [Volumetest;sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)]; %Calculating volume [in m^3]

%Penetration depth
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+10);  %Just staying away from ocean
pene_depthtest = [pene_depthtest;Grid.p.yc(pene_ind(1))];

i = i+100;
end

hh=figure()
yline(phiOrig/1e9,'-','Linewidth',5,'color',red);
hold on
xlim([0, max(tVec)])
plot(Time_arrtest,Volumetest/1e9,'-','Linewidth',5,'color',blue);
plot(tVec,phiDrain1Vec/1e9,'-','Linewidth',5,'color',green);
xlabel('Time [years]');
ylabel('Melt volume [km$$^3$$]');
legend('Initial total','Left','Delivered','location','northeast');
saveas(hh,sprintf('../figures/1.85e-10_vol.png')); 
saveas(hh,sprintf('../figures/1.85e-10_temperature_vol.pdf'));

drain_V_array = [drain_V_array; phiDrain1Vec(end)/1e9]; %in km3
kc_array    = [kc_array; 1.85e-10]; %in m2


%kc1.85e-12m2
d = 10e3
load("../Output/Europa03321_eta0_14kc1.85e-12_Ea_50_output_1C.mat"); %loading file
interface = (find(Grid.p.yc>0)); interface = interface(1,1); %First cell of no ocean
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean
Volumetest = [sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)];
Time_arrtest = [tVec];
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+10);  %Just staying away from ocean
pene_depthtest= [Grid.p.yc(pene_ind(1))];
i = 100;
while i<=23100
%Individual file
load(sprintf("../Output/Europa03321_eta0_14kc1.85e-12_Ea_50_output_%sC.mat",num2str(i))); %loading file
Time_arrtest = [Time_arrtest;tVec(i)]; %Time array [years]
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean

Volumetest = [Volumetest;sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)]; %Calculating volume [in m^3]

%Penetration depth
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+10);  %Just staying away from ocean
pene_depthtest = [pene_depthtest;Grid.p.yc(pene_ind(1))];

i = i+100;
end

hh=figure()
yline(phiOrig/1e9,'-','Linewidth',5,'color',red);
hold on
xlim([0, max(tVec)])
plot(Time_arrtest,Volumetest/1e9,'-','Linewidth',5,'color',blue);
plot(tVec,phiDrain1Vec/1e9,'-','Linewidth',5,'color',green);
xlabel('Time [years]');
ylabel('Melt volume [km$$^3$$]');
legend('Initial total','Left','Delivered','location','northeast');
saveas(hh,sprintf('../figures/1.85e-12_vol.png')); 
saveas(hh,sprintf('../figures/1.85e-12_temperature_vol.pdf'));

drain_V_array = [drain_V_array; phiDrain1Vec(end)/1e9]; %in km3
kc_array    = [kc_array; 1.85e-12]; %in m2


%kc1.85e-14m2
d = 10e3
load("../Output/Europa03321_eta0_14kc1.85e-14_Ea_50_output_1C.mat"); %loading file
interface = (find(Grid.p.yc>0)); interface = interface(1,1); %First cell of no ocean
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean
phi_arr_delivered = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr_delivered(interface:end,:) = 0; %Zeroing out the melt in domain
Volumetest_delivered = [sum(sum(phi_arr_delivered(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)];
Volumetest = [sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)];
Time_arrtest = [tVec];
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+10);  %Just staying away from ocean
pene_depthtest= [Grid.p.yc(pene_ind(1))];
i = 100;
while i<=33000
%Individual file
load(sprintf("../Output/Europa03321_eta0_14kc1.85e-14_Ea_50_output_%sC.mat",num2str(i))); %loading file
Time_arrtest = [Time_arrtest;tVec(i)]; %Time array [years]
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr_delivered = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean
phi_arr_delivered(interface:end,:) = 0; %Zeroing out above ocean

Volumetest = [Volumetest;sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)]; %Calculating volume [in m^3]
Volumetest_delivered = [Volumetest_delivered;sum(sum(phi_arr_delivered(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)]; %Calculating volume [in m^3]

%Penetration depth
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+10);  %Just staying away from ocean
pene_depthtest = [pene_depthtest;Grid.p.yc(pene_ind(1))];

i = i+100;
end

hh=figure()
yline(phiOrig/1e9,'-','Linewidth',5,'color',red);
hold on
xlim([0, max(tVec)])
plot(Time_arrtest,Volumetest/1e9,'-','Linewidth',5,'color',blue);
plot(tVec,phiDrain2Vec/1e9,'-','Linewidth',5,'color',green);
plot(Time_arrtest,Volumetest_delivered/1e9,'--','Linewidth',5,'color',black);
xlabel('Time [years]');
ylabel('Melt volume [km$$^3$$]');
legend('Initial total','Left','Delivered','Delivered new','location','northeast');
saveas(hh,sprintf('../figures/1.85e-14_vol.png')); 
xlim([0,1000])
saveas(hh,sprintf('../figures/1.85e-14_temperature_vol.pdf'));

%drain_V_array = [drain_V_array; phiDrain2Vec(end)/1e9]; %in km3
drain_V_array = [drain_V_array; max(Volumetest_delivered)/1e9]; %in km3
kc_array    = [kc_array; 1.85e-14]; %in m2

%kc1.85e-16m2
d = 10e3
load("../Output/Europa03321_eta0_14kc1.85e-16_Ea_50_output_1C.mat"); %loading file
interface = (find(Grid.p.yc>0)); interface = interface(1,1); %First cell of no ocean
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean
Volumetest = [sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)];
phi_arr_delivered = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr_delivered(interface:end,:) = 0; %Zeroing out the melt in domain
Volumetest_delivered = [sum(sum(phi_arr_delivered(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)];
Time_arrtest = [tVec];
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+10);  %Just staying away from ocean
pene_depthtest= [Grid.p.yc(pene_ind(1))];
i = 100;
while i<=39300
%Individual file
load(sprintf("../Output/Europa03321_eta0_14kc1.85e-16_Ea_50_output_%sC.mat",num2str(i))); %loading file
Time_arrtest = [Time_arrtest;tVec(i)]; %Time array [years]
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean
phi_arr_delivered = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr_delivered(interface:end,:) = 0; %Zeroing out above ocean

Volumetest = [Volumetest;sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)]; %Calculating volume [in m^3]
Volumetest_delivered = [Volumetest_delivered;sum(sum(phi_arr_delivered(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)]; %Calculating volume [in m^3]

%Penetration depth
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+10);  %Just staying away from ocean
pene_depthtest = [pene_depthtest;Grid.p.yc(pene_ind(1))];

i = i+100;
end

hh=figure()
yline(phiOrig/1e9,'-','Linewidth',5,'color',red);
hold on
xlim([0, max(tVec)])
plot(Time_arrtest,Volumetest/1e9,'-','Linewidth',5,'color',blue);
plot(tVec,phiDrain2Vec/1e9,'-','Linewidth',5,'color',green);
plot(Time_arrtest,Volumetest_delivered/1e9,'--','Linewidth',5,'color',black);
xlabel('Time [years]');
ylabel('Melt volume [km$$^3$$]');
legend('Initial total','Left','Delivered','Delivered new','location','northeast');
saveas(hh,sprintf('../figures/1.85e-16_vol.png')); 
saveas(hh,sprintf('../figures/1.85e-16_temperature_vol.pdf'));

%drain_V_array = [drain_V_array; phiDrain2Vec(end)/1e9]; %in km3
drain_V_array = [drain_V_array; max(Volumetest_delivered)/1e9]; %in km3
kc_array    = [kc_array; 1.85e-16]; %in m2


hhh=figure
hhh.Units = 'centimeters'
hhh.Position = [1,1,20,15];
semilogx(kc_array, drain_V_array,'ro','Linewidth',2,'color',red,'MarkerSize',12);
xlabel('$$\textrm{k}_0 [\textrm{m}^2]$$');
ylabel('Melt delivered to ocean [km$$^3$$]');
xticks([1e-16, 1e-14, 1e-12, 1e-10, 1e-8])
ylim([0,30])
saveas(hhh,sprintf('../figures/melt_deliveredvsKc.png')); 
saveas(hhh,sprintf('../figures/melt_deliveredvsKc.pdf'));




%New part
%Combined
tstamp_DS = [1,20,40,60,80,100,120,500,1000,2000]; 
timestamp_DS = []; Volume_timestamp_DS = []; pene_depth_timestamp_DS = [];
tstamp_M = [1,20,40,60,80,100,125,150,200,620];
timestamp_M = []; Volume_timestamp_M = []; pene_depth_timestamp_M = [];
tstamp_S = [1,500,1000,1500,2000,3000,4000,5000,6000,12000];
timestamp_S = []; Volume_timestamp_S = []; pene_depth_timestamp_S = [];

for j = 1:size(tstamp_DS,2)
    load(sprintf("../Output/Europa03321_eta0_14kc1.85e-08_Ea_50_output_%sC.mat",num2str(tstamp_DS(j))))
    timestamp_DS = [timestamp_DS,tVec(end)];  
    
    if j ==1
        interface = (find(Grid.p.yc>0)); interface = interface(1,1); %First cell of no ocean
    end
    phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
    phi_arr(1:interface+1,:) = 0; %Zeroing out ocean
    
    Vol = sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3);

    if j>1
        if Vol > Volume_timestamp_DS(end)
            Vol = Volume_timestamp_DS(end);
        end
    end
    Volume_timestamp_DS = [Volume_timestamp_DS,Vol]; %Calculating volume [in m^3]
    
    %Penetration depth
    [Xc,Yc] = meshgrid(Grid.p.xc,Grid.p.yc); Yc = Yc(:);
    Yc_wet = Yc(find(phi>1e-5)); pene_depthdummy = min(Yc_wet(Yc_wet>0));
    %pene_depthS = [pene_depthS;Grid.p.yc(pene_ind(1))];
    pene_depth_timestamp_DS = [pene_depth_timestamp_DS; pene_depthdummy];

end

hh=figure
hh.Units = 'centimeters'
hh.Position = [1,1,20,15];
set(gcf,'visible','on')
yyaxis left
plot(Time_arrDS,VolumeDS/1e9,'-','Linewidth',3,'color',blue)
hold on
plot(timestamp_DS,Volume_timestamp_DS/1e9,'Linewidth',1.5,'Marker','o',MarkerFaceColor=blue,MarkerSize=8,LineStyle='none',MarkerEdgeColor=black)

%semilogx(Time_arrM,VolumeM/1e9,'-.','Linewidth',3,'color',blue)
hold on
%semilogx(Time_arrS,VolumeS/1e9,'--','Linewidth',3,'color',blue)
xlabel('Time [years]');
ylabel('Melt volume [km$$^3$$]','Color',blue);
ax = gca
ax.YColor = blue;
yyaxis right
plot(Time_arrDS,10 - pene_depthDS*d/1e3,'--','Linewidth',3,'color',red)
hold on
plot(timestamp_DS,10 - pene_depth_timestamp_DS*d/1e3,'Linewidth',1.5,'Marker','o',MarkerFaceColor=red,MarkerSize=8,LineStyle='none',MarkerEdgeColor=black)
%xlim([1e-4,1e5])
%xticks([1e-4,1e-2,1,1e2,1e4])
set(gca, 'YDir','reverse')
ax = gca
ax.YColor = red;
ylabel('Penetration depth [km]','Color',red);
%lh = legend('$$\textrm{k}_0=1.85\cdot10^{-8}\textrm{m}^2$$','$$\textrm{k}_0=1.85\cdot 10^{-10}\textrm{m}^2$$','$$\textrm{k}_0=1.85\cdot 10^{-16}\textrm{m}^2$$','location','southwest');
%set( lh, 'Box', 'off' ) ;
xlim([-0.001,0.15])
saveas(hh,sprintf('../figures/EuropaCombinedVolandpene_newDS.png')); 
saveas(hh,sprintf('../figures/EuropaCombinedVolandpene_newDS.pdf')); 

set(gca, 'XScale', 'log');
xlim([-0.001,1e3])
saveas(hh,sprintf('../figures/EuropaCombinedVolandpene_newDS_log.png')); 
saveas(hh,sprintf('../figures/EuropaCombinedVolandpene_newDS_log.pdf')); 

%%
%Middle Case
for j = 1:size(tstamp_M,2)
    load(sprintf("../Output/Europa03321_eta0_14kc1.85e-10_Ea_50_output_%sC.mat",num2str(tstamp_M(j))))
    timestamp_M = [timestamp_M,tVec(end)];  
    
    if j ==1
        interface = (find(Grid.p.yc>0)); interface = interface(1,1); %First cell of no ocean
    end
    phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
    phi_arr(1:interface+1,:) = 0; %Zeroing out ocean
    
    Vol = sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3);

    if j>1
        if Vol > Volume_timestamp_M(end)
            Vol = Volume_timestamp_M(end);
        end
    end
    Volume_timestamp_M = [Volume_timestamp_M,Vol]; %Calculating volume [in m^3]
    
    %Penetration depth
    [Xc,Yc] = meshgrid(Grid.p.xc,Grid.p.yc); Yc = Yc(:);
    Yc_wet = Yc(find(phi>1e-5)); pene_depthdummy = min(Yc_wet(Yc_wet>0));
    %pene_depthS = [pene_depthS;Grid.p.yc(pene_ind(1))];
    pene_depth_timestamp_M = [pene_depth_timestamp_M; pene_depthdummy];

end

hh=figure
hh.Units = 'centimeters'
hh.Position = [1,1,20,15];
set(gcf,'visible','on')
yyaxis left
plot(Time_arrM,VolumeM/1e9,'-','Linewidth',3,'color',blue)
hold on
plot(timestamp_M,Volume_timestamp_M/1e9,'Linewidth',1.5,'Marker','o',MarkerFaceColor=blue,MarkerSize=8,LineStyle='none',MarkerEdgeColor=black)
hold on
xlabel('Time [years]');
ylabel('Melt volume [km$$^3$$]','Color',blue);
ax = gca
ax.YColor = blue;
yyaxis right
plot(Time_arrM,10 - pene_depthM*d/1e3,'--','Linewidth',3,'color',red)
hold on
plot(timestamp_M,10 - pene_depth_timestamp_M*d/1e3,'Linewidth',1.5,'Marker','o',MarkerFaceColor=red,MarkerSize=8,LineStyle='none',MarkerEdgeColor=black)
set(gca, 'YDir','reverse')
ax = gca
ax.YColor = red;
ylabel('Penetration depth [km]','Color',red);
xlim([-0.01,5])
saveas(hh,sprintf('../figures/EuropaCombinedVolandpene_newM.png')); 
saveas(hh,sprintf('../figures/EuropaCombinedVolandpene_newM.pdf')); 

set(gca, 'XScale', 'log');
xlim([-0.001,1e4])
saveas(hh,sprintf('../figures/EuropaCombinedVolandpene_newM_log.png')); 
saveas(hh,sprintf('../figures/EuropaCombinedVolandpene_newM_log.pdf')); 



%%
%Stokes Case
for j = 1:size(tstamp_S,2)
    load(sprintf("../Output/Europa03321_eta0_14kc1.85e-16_Ea_50_output_%sC.mat",num2str(tstamp_S(j))))
    timestamp_S = [timestamp_S,tVec(end)];  
    
    if j ==1
        interface = (find(Grid.p.yc>0)); interface = interface(1,1); %First cell of no ocean
    end
    phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
    phi_arr(1:interface+1,:) = 0; %Zeroing out ocean
    
    Vol = sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3);

    if j>1
        if Vol > Volume_timestamp_S(end)
            Vol = Volume_timestamp_S(end);
        end
    end
    Volume_timestamp_S = [Volume_timestamp_S,Vol]; %Calculating volume [in m^3]
    
    %Penetration depth
    [Xc,Yc] = meshgrid(Grid.p.xc,Grid.p.yc); Yc = Yc(:);
    Yc_wet = Yc(find(phi>1e-5)); pene_depthdummy = min(Yc_wet(Yc_wet>0));
    %pene_depthS = [pene_depthS;Grid.p.yc(pene_ind(1))];
    pene_depth_timestamp_S = [pene_depth_timestamp_S; pene_depthdummy];

end

hh=figure
hh.Units = 'centimeters'
hh.Position = [1,1,20,15];
set(gcf,'visible','on')
yyaxis left
plot(Time_arrS,VolumeS/1e9,'-','Linewidth',3,'color',blue)
hold on
plot(timestamp_S,Volume_timestamp_S/1e9,'Linewidth',1.5,'Marker','o',MarkerFaceColor=blue,MarkerSize=8,LineStyle='none',MarkerEdgeColor=black)
hold on
xlabel('Time [years]');
ylabel('Melt volume [km$$^3$$]','Color',blue);
ax = gca
ax.YColor = blue;
yyaxis right

%Interp triangle
triang_x = [timestamp_S(6)-50,timestamp_S(4)-50];  % chose location
triang_y = [10 - pene_depth_timestamp_S(6)*d/1e3,10 - pene_depth_timestamp_S(4)*d/1e3];

aaa1 = plot(triang_x([1,2,2]), triang_y([1,1,2]), 'k-','Linewidth',1.5);
aaa2 = plot(triang_x, triang_y, 'k-','Linewidth',1.5);

plot(Time_arrS,10 - pene_depthS*d/1e3,'--','Linewidth',3,'color',red)
hold on
plot(timestamp_S,10 - pene_depth_timestamp_S*d/1e3,'Linewidth',1.5,'Marker','o',MarkerFaceColor=red,MarkerSize=8,LineStyle='none',MarkerEdgeColor=black)
set(gca, 'YDir','reverse')
ax = gca
ax.YColor = red;
ylabel('Penetration depth [km]','Color',red);
ylim([3,10.2])
xlim([-5,1500])
saveas(hh,sprintf('../figures/EuropaCombinedVolandpene_newS.png')); 
saveas(hh,sprintf('../figures/EuropaCombinedVolandpene_newS.pdf')); 

delete(aaa1)
delete(aaa2)
set(gca, 'XScale', 'log');
xlim([-0.001,1e4])
saveas(hh,sprintf('../figures/EuropaCombinedVolandpene_newS_log.png')); 
saveas(hh,sprintf('../figures/EuropaCombinedVolandpene_newScl_log.pdf')); 

