%Calculating volume at each time
    set(groot,'defaultAxesFontName','Times')
    set(groot,'defaultAxesFontSize',20)
    set(groot,'defaulttextinterpreter','latex')
    set(groot,'defaultAxesTickLabelInterpreter','latex')
    set(groot,'defaultLegendInterpreter','latex')
    set(groot, 'DefaultFigureVisible', 'on');


d = 20e3;   %Thickness of ice shell [in m]   


load("../Output/impact_Wakita_eta0_14_Ea_50_output_100Tlayer-17C.mat"); %loading file
interface = (find(Grid.p.yc>0)); interface = interface(1,1); %First cell of no ocean
VolumeDS = [];
Time_arrDS = [];
pene_depthDS= [];
i = 100;

%Darcy-Stokes
while i<=50000
%Individual file
load(sprintf("../Output/Darcy-Stokes-kc5.6e-11_impact_Wakita_eta0_14_Ea_50_output_100Tlayer-17C/impact_Wakita_eta0_14_Ea_50_output_%sTlayer-17C.mat",num2str(i))); %loading file
Time_arrDS = [Time_arrDS;tVec(i)]; %Time array [years]
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean

VolumeDS = [VolumeDS;sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)]; %Calculating volume [in m^3]

%Penetration depth
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+10);  %Just staying away from ocean
pene_depthDS = [pene_depthDS;Grid.p.yc(pene_ind(1))];

i = i+100;
end

%Stokes
VolumeS = [];
Time_arrS = [];
pene_depthS= [];
i = 100;
while i<=29800
%Individual file
load(sprintf("../Output/impact_Wakita-Stokes_eta0_14_Ea_50_output_%sTlayer-17C.mat",num2str(i))); %loading file

Time_arrS = [Time_arrS;tVec(i)]; %Time array [years]
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface+1,:) = 0; %Zeroing out ocean

%Calculating volume
VolumeS = [VolumeS;sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)]; %Calculating volume [in m^3]

%Penetration depth
pene_ind= find((sum(phi_arr,2)));   %Finding penetration depth index
pene_ind(pene_ind>interface+100);  %Just staying away from ocean
pene_depthS = [pene_depthS;Grid.p.yc(pene_ind(1))];

i = i+100;
end

%Plotting
%Volume plot
figure()
semilogx(Time_arrDS,VolumeDS/1e9,'r-','Linewidth',5)
hold on
semilogx(Time_arrS,VolumeS/1e9,'b-','Linewidth',5)
xlabel('Time [years]');
ylabel('Melt volume [km$$^3$$]');
legend('Darcy-Stokes','Stokes');

%Volume plot
figure()
semilogx(Time_arrDS,pene_depthDS/1e3,'r-','Linewidth',5)
hold on
semilogx(Time_arrS,pene_depthS/1e3,'b-','Linewidth',5)
xlabel('Time [years]');
ylabel('Penetration depth [km]');
legend('Darcy-Stokes','Stokes');


%Combined
hh=figure()
yyaxis left
semilogx(Time_arrDS,VolumeDS/1e9,'k-','Linewidth',5)
hold on
semilogx(Time_arrS,VolumeS/1e9,'k--','Linewidth',5)
semilogx(Time_arrDS,VolumeDS/1e9,'-','Linewidth',5,'color',[0,114/255,189/255])
hold on
semilogx(Time_arrS,VolumeS/1e9,'--','Linewidth',5,'color',[0,114/255,189/255])
xlabel('Time [years]');
ylabel('Melt volume [km$$^3$$]');
yyaxis right
semilogx(Time_arrDS,60 - pene_depthDS*d/1e3,'r-','Linewidth',5,'color',[217/255,83/255,25/255])
hold on
semilogx(Time_arrS,60 - pene_depthS*d/1e3,'r--','Linewidth',5,'color',[217/255,83/255,25/255])
%ylim([0,45])
set(gca, 'YDir','reverse')
ylabel('Penetration depth [km]');
legend('Darcy-Stokes','Stokes','location','southwest');
saveas(hh,sprintf('../figures/CombinedVolandpene.png')); 
saveas(hh,sprintf('../figures/CombinedVolandpene.png.pdf')); 


%Combined Stokes
hh=figure()
yyaxis left
semilogx(Time_arrS,VolumeS/1e9,'-','Linewidth',5,'color',[0,114/255,189/255])
hold on
xlabel('Time [years]');
ylabel('Melt volume [km$$^3$$]');
yyaxis right
semilogx(Time_arrS,60 - pene_depthS*d/1e3,'r-','Linewidth',5,'color',[217/255,83/255,25/255])
hold on
%ylim([0,45])
set(gca, 'YDir','reverse')
ylabel('Penetration depth [km]');
saveas(hh,sprintf('../figures/CombinedStokes.png')); 
saveas(hh,sprintf('../figures/CombinedStokes.pdf')); 



%Combined Darcy-Stokes
hh=figure()
yyaxis left
semilogx(Time_arrDS,VolumeDS/1e9,'-','Linewidth',5,'color',[0,114/255,189/255])
hold on
xlabel('Time [years]');
ylabel('Melt volume [km$$^3$$]');
yyaxis right
semilogx(Time_arrDS,60 - pene_depthDS*d/1e3,'r-','Linewidth',5,'color',[217/255,83/255,25/255])
hold on
set(gca, 'YDir','reverse')
ylabel('Penetration depth [km]');
saveas(hh,sprintf('../figures/CombinedDarcyStokes.png')); 
saveas(hh,sprintf('../figures/CombinedDarcyStokes.pdf')); 


%Combined Vol
hh=figure()
semilogx(Time_arrDS,VolumeDS/1e9,'-','Linewidth',5,'color',[0,114/255,189/255])
hold on
semilogx(Time_arrS,VolumeS/1e9,'--','Linewidth',5,'color',[0,114/255,189/255])
xlabel('Time [years]');
ylabel('Melt volume [km$$^3$$]');
legend('Darcy-Stokes','Stokes','location','southwest');
saveas(hh,sprintf('../figures/CombinedVol.png')); 
saveas(hh,sprintf('../figures/CombinedVol.pdf')); 
