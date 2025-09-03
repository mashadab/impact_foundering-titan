%% set up plot commands
set(groot,'defaulttextinterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex')
set(groot,'defaultTextFontSize',20)
set(groot, 'DefaultAxesFontSize', 20);
addpath ../model_and_dependencies/
clear all; close all;
fs = 18;
fs_small = 12;
%% set up grid
Xmax = 20; %maximum 
Ymax = 60; %maximum  
grRes = 50; % grid resolution in radial direction
grZ = 180; % grid resolution in radial direction
d = 20;
Gridp.xmin = 0; Gridp.xmax = Xmax/d; Gridp.Nx = grRes; %radial direction
Gridp.ymin = -3/5; Gridp.ymax = Ymax/d;  Gridp.Ny = grZ; %vertical direction
Gridp.geom = 'cylindrical_rz';
Grid = build_stokes_grid_cyl(Gridp);
[X,Y] = meshgrid(Grid.p.xc,Grid.p.yc);
angle = 0;

%% make plot dimensions
fp = '../Output/';

%Original
%inds1 = [1,2500,3000,4500,5500];%Stokes-ish
%inds2 = [1,5000,10000,30000,42000]; %Darcy-Stokes

inds1_arr = [1,1000,2500,3000,3500,4000,4500,5500,10000,50000];%Stokes-ish
inds2_arr = [1,5000,10000,20000,22500,24000,25000,30000,35000,42000]; %Darcy-Stokes
caption = 1; %0 = no caption, 1=caption

for kk=1:2
    if kk==1
        inds1 = inds1_arr(1:5);
        inds2 = inds2_arr(1:5);
    else
        inds1 = inds1_arr(6:end);
        inds2 = inds2_arr(6:end);    
    end 
    inds = [inds1; inds2];
    
    f = figure;
    f.Units = 'centimeters';
    % [left bottom width height]
    f.Position = [1,1,22,9.5*2+0.5]; % should be 19 wide
    axs = cell(2,10);
    
    lMarg = 1;
    bMarg = 1;
    wid = 7;
    hMarg = -4.65;
    vMarg2 = -13;
    vMarg1 = 1;
    tFrame = zeros(2,5);
    for i = 1:2
        for j = 1:5
            if i ==1 %Close to Stokes
                load(sprintf("../Output/Simple-test_eta0_14kc1.85e-10_Ea_50_output_%sC.mat",num2str(inds1(j))))
            end
            if i==2  %Close to Darcy-Stokes
                load(sprintf("../Output/Simple-test_eta0_14kc1.85e-08_Ea_50_output_%sC.mat",num2str(inds2(j))))
            end
            %load([fp fns{i} '/i' num2str(inds(i,j)) '.mat'],'T','phi','tVec');
            TPlot = reshape(T,Grid.p.Ny,Grid.p.Nx)*173+100;
            phiPlot = reshape(phi,Grid.p.Ny,Grid.p.Nx);
    
    %         TPlot(phiPlot > 1e-2) = nan;
            phiPlot(phiPlot < 1e-16) = nan;
            TPlot(TPlot>273)=273; %%%%%new
            TPlot(TPlot<94)=94; %%%%%new
            axs{i,j} = subplot(4,5,j+(i-1)*10);
            hold on
            contourf(X*d,Y*d,TPlot,40,'linestyle','none');
            caxis([94 273]);
            axis equal
            xlim([0 20]);
            ylim([0 60]);
            axs{i,j}.XTickLabel = [];
            axs{i,j}.YTickLabel = [];
            axs{i,j}.Units = 'centimeters';
            colormap(axs{i,j},'bone');
    
            if i == 1 && j == 5
                c = colorbar('EastOutside');
                c.Label.Interpreter = 'latex';
                c.TickLabelInterpreter = 'latex';
                c.Label.FontSize = fs;
                c.Label.String = 'Temperature, K';
                if i == 1
                    cs{1} = c;
                else
                    cs{2} = c;
                end
             end
    
            axs{i,5+j} = subplot(4,5,(i-1)*10+5+j);
            axs{i,5+j}.Colormap = colormap(flipud(parula));
            contourf(X*d,Y*d,phiPlot,40,'linestyle','none');
            axis equal
            xlim([0 20]);
            ylim([0 60]);
            axs{i,5+j}.XTickLabel = [];
            axs{i,5+j}.YTickLabel = [];
            axs{i,5+j}.Units = 'centimeters';
            caxis([0 1]);
           
            
            if i == 1 && j == 5
                c = colorbar('EastOutside');
                c.Label.Interpreter = 'latex';
                c.TickLabelInterpreter = 'latex';
                c.Label.FontSize = fs;%8;
                c.Label.String = 'Melt fraction';
                if i == 1
                    cs{3} = c;
                else
                    cs{4} = c;
                end
            end
            tFrame(i,j) = round(tVec(end),2);
            
        end
    end
    
    
    %% rearrange
    for i = 1:2
        for j = 1:5
            if i == 1
                axs{i,j}.Position = [lMarg+(j-1)*(wid+hMarg),bMarg+3*(wid+vMarg1)+vMarg2 , wid, wid];
            elseif i == 2
                axs{i,j}.Position = [lMarg+(j-1)*(wid+hMarg),bMarg+1*(wid+vMarg1), wid, wid];
            end
            
            if i == 1
                axs{i,5+j}.Position = [lMarg+(j-1)*(wid+hMarg),bMarg+2*(wid+vMarg1)+vMarg2 , wid, wid];
            elseif i == 2
                axs{i,5+j}.Position = [lMarg+(j-1)*(wid+hMarg),bMarg, wid, wid];
            end
        end
    end
            
    yTikLab = {'0','20','40','60'};
    xTikLab = {'0','5','10','15','20'};
    tik1 = [0,10,20];
    tik2 = [0,20,40,60];
    
    
    % l b w h
    i = 2;
    j = 6;
    % axs{i,j}.Position = [lMarg,bMarg,wid,wid];
    axs{i,j}.YTick = tik2;
    axs{i,j}.YTickLabel = yTikLab;
    axs{i,j}.XTick = tik1;
    axs{i,j}.XTickLabel = xTikLab;
    axs{i,j}.YLabel.String = 'z-dir [km]';
    axs{i,j}.XLabel.String = 'r [km]';
    xtickangle(axs{i,j}, angle); axs{i,j}.TickDir = 'in';   axs{i,j}.TickLength = [0.02 0.02];     % ticks pointing inward % make them longer axs{i,j}.TickDir = 'in';      % ticks pointing inward axs{i,j}.TickLength = [0.02 0.02];  % make them longer
    
    axs{i,j}.YLabel.FontSize = fs;
    axs{i,j}.XLabel.FontSize = fs;
    
    i = 2;
    j = 1;
    % axs{i,j}.Position = [lMarg,bMarg+wid+vMarg1,wid,wid];
    axs{i,j}.YTick = tik2;
    axs{i,j}.YTickLabel = yTikLab;
    axs{i,j}.YLabel.String = 'z-dir [km]';
    
    axs{i,j}.YLabel.FontSize = fs;
    xtickangle(axs{i,j}, angle); axs{i,j}.TickDir = 'in';   axs{i,j}.TickLength = [0.02 0.02];     % ticks pointing inward % make them longer
    
    i = 1;
    j = 1;
    % axs{i,j}.Position = [lMarg,bMarg+2*(wid+vMarg1)+vMarg2,wid,wid];
    axs{i,j}.YTick = tik2;
    axs{i,j}.YTickLabel = yTikLab;
    axs{i,j}.YLabel.String = 'z-dir [km]';
    axs{i,j}.XTickLabel = [];
    axs{i,j}.YLabel.FontSize = fs;
    xtickangle(axs{i,j}, angle); axs{i,j}.TickDir = 'in';   axs{i,j}.TickLength = [0.02 0.02];     % ticks pointing inward % make them longer
    
    i = 1;
    j = 6;
    % axs{i,j}.Position = [lMarg,bMarg+3*(wid+vMarg1)+vMarg2,wid,wid];
    axs{i,j}.YTick = tik2;
    axs{i,j}.YTickLabel = yTikLab;
    axs{i,j}.YLabel.String = 'z-dir [km]';
    axs{i,j}.XTickLabel = [];
    axs{i,j}.YLabel.FontSize = fs;
    xtickangle(axs{i,j}, angle); axs{i,j}.TickDir = 'in';   axs{i,j}.TickLength = [0.02 0.02];     % ticks pointing inward % make them longer
    
    %colorbar loc
    for i = 1:length(cs)
        cs{i}.Position(3) = 0.01;
        cs{i}.Position(1) = 0.85;
    end
    %%
    figLab = 'abcdefghij';
    for i = 1:5
        colormap(axs{1,i},'bone'); 
        colormap(axs{2,i},'bone'); 
        set([axs{2,i},axs{2,i+5}],'Position',axs{1,i+5}.Position)
        axs{1,i+5}.Visible = 'off';
        axs{2,i+5}.Visible = 'off';
        set([axs{1,i},axs{1,i+5}],'Position',axs{1,i}.Position)
        if caption == 1
        text(axs{2,i+5},19,0,[num2str(tFrame(2,i)) ' yrs'],...
            'FontSize',fs_small,'HorizontalAlignment','right','Color',[1,0,0])
        text(axs{1,i+5},19,0,[num2str(tFrame(1,i)) ' yrs'],...
            'FontSize',fs_small,'HorizontalAlignment','right','Color',[1,0,0])
        text(axs{1,i+5},19,9.9,figLab(i),'fontweight','bold','interpreter','tex',...
            'FontSize',fs_small,'HorizontalAlignment','right','Color',[1,0,0])
        text(axs{2,i+5},19,9.9,figLab(i+5),'fontweight','bold','interpreter','tex',...
            'FontSize',fs_small,'HorizontalAlignment','right','Color',[1,0,0])
        end
    end
    
    
    i = 2;
    j = 1;
    axs{i,j}.XTick = tik1;
    axs{i,j}.XTickLabel = xTikLab;
    axs{i,j}.XLabel.String = 'r [km]';
    axs{i,j}.XLabel.FontSize = fs;
    xtickangle(axs{i,j}, angle); axs{i,j}.TickDir = 'in';   axs{i,j}.TickLength = [0.02 0.02];     % ticks pointing inward % make them longer
    
    i = 2;
    j = 2;
    axs{i,j}.XTick = tik1;
    axs{i,j}.YTick = tik2;
    axs{i,j}.XTickLabel = xTikLab;
    %axs{i,j}.XLabel.String = 'r [km]';
    axs{i,j}.XLabel.FontSize = fs;
    xtickangle(axs{i,j}, angle); axs{i,j}.TickDir = 'in';   axs{i,j}.TickLength = [0.02 0.02];     % ticks pointing inward % make them longer
    axs{i,j}.XTickLabel = [];
    
    i = 2;
    j = 3;
    axs{i,j}.XTick = tik1;
    axs{i,j}.YTick = tik2;
    axs{i,j}.XTickLabel = xTikLab;
    %axs{i,j}.XLabel.String = 'r [km]';
    axs{i,j}.XLabel.FontSize = fs;
    xtickangle(axs{i,j}, angle); axs{i,j}.TickDir = 'in';   axs{i,j}.TickLength = [0.02 0.02];     % ticks pointing inward % make them longer
    axs{i,j}.XTickLabel = [];
    
    i = 2;
    j = 4;
    axs{i,j}.XTick = tik1;
    axs{i,j}.YTick = tik2;
    axs{i,j}.XTickLabel = xTikLab;
    %axs{i,j}.XLabel.String = 'r [km]';
    axs{i,j}.XLabel.FontSize = fs;
    xtickangle(axs{i,j}, angle); axs{i,j}.TickDir = 'in';   axs{i,j}.TickLength = [0.02 0.02];     % ticks pointing inward % make them longer
    axs{i,j}.XTickLabel = [];
    
    i = 2;
    j = 5;
    axs{i,j}.XTick = tik1;
    axs{i,j}.YTick = tik2;
    axs{i,j}.XTickLabel = xTikLab;
    %axs{i,j}.XLabel.String = 'r [km]';
    axs{i,j}.XLabel.FontSize = fs;
    xtickangle(axs{i,j}, angle); axs{i,j}.TickDir = 'in';   axs{i,j}.TickLength = [0.02 0.02];     % ticks pointing inward % make them longer
    axs{i,j}.XTickLabel = [];
    
    i = 1;
    j = 2;
    axs{i,j}.XTick = tik1;
    axs{i,j}.YTick = tik2;
    xtickangle(axs{i,j}, angle); axs{i,j}.TickDir = 'in';   axs{i,j}.TickLength = [0.02 0.02];     % ticks pointing inward % make them longer
    
    i = 1;
    j = 3;
    axs{i,j}.XTick = tik1;
    axs{i,j}.YTick = tik2;
    xtickangle(axs{i,j}, angle); axs{i,j}.TickDir = 'in';   axs{i,j}.TickLength = [0.02 0.02];     % ticks pointing inward % make them longer
    
    i = 1;
    j = 4;
    axs{i,j}.XTick = tik1;
    axs{i,j}.YTick = tik2;
    xtickangle(axs{i,j}, angle); axs{i,j}.TickDir = 'in';   axs{i,j}.TickLength = [0.02 0.02];     % ticks pointing inward % make them longer
    
    i = 1;
    j = 5;
    axs{i,j}.XTick = tik1;
    axs{i,j}.YTick = tik2;
    xtickangle(axs{i,j}, angle); axs{i,j}.TickDir = 'in';   axs{i,j}.TickLength = [0.02 0.02];     % ticks pointing inward % make them longer
    
    f.PaperPosition = [0.0177,4.6440,4.4646,5.7121];
    Filename = sprintf('FinalStokesvsDarcyStokesEuropa_vectorpart_%s',num2str(kk));
    set(gcf,'PaperPositionMode','auto'); 
    saveas(gcf,[Filename '.pdf']);
    
    Filename = sprintf('FinalStokesvsDarcyStokesEuropa_vectorpart_%s',num2str(kk));
    set(gcf,'PaperPositionMode','auto');
    print(gcf, Filename, '-dpdf', '-painters');
end



%% 
%Stokes only

%Stokes-ish
inds1_arr = [1,1500,3000,6000,11000,20000,40000,70000,85000,100000];%Stokes
inds2_arr = [1,1000,2500,3000,3500,4000,4500,5500,10000,50000];%DarcyStokes-sh e-10
caption = 1; %0 = no caption, 1=caption

for kk=1:2
    if kk==1
        inds1 = inds1_arr(1:5);
        inds2 = inds2_arr(1:5);
    else
        inds1 = inds1_arr(6:end);
        inds2 = inds2_arr(6:end);    
    end 
    inds = [inds1; inds2];
    
    f = figure;
    f.Units = 'centimeters';
    % [left bottom width height]
    f.Position = [1,1,22,9.5*2+0.5]; % should be 19 wide
    axs = cell(2,10);
    
    lMarg = 1;
    bMarg = 1;
    wid = 7;
    hMarg = -4.65;
    vMarg2 = -13;
    vMarg1 = 1;
    tFrame = zeros(2,5);
    for i = 1:2
        for j = 1:5
            if i ==1 %Close to Stokes
                load(sprintf("../Output/Simple-test_eta0_14kc1.85e-16_Ea_50_output_%sC.mat",num2str(inds1(j))))
            end
            if i==2  %Close to Darcy-Stokes
                load(sprintf("../Output/Simple-test_eta0_14kc1.85e-10_Ea_50_output_%sC.mat",num2str(inds2(j))))
            end
            %load([fp fns{i} '/i' num2str(inds(i,j)) '.mat'],'T','phi','tVec');
            TPlot = reshape(T,Grid.p.Ny,Grid.p.Nx)*173+100;
            phiPlot = reshape(phi,Grid.p.Ny,Grid.p.Nx);
    
    %         TPlot(phiPlot > 1e-2) = nan;
            phiPlot(phiPlot < 1e-16) = nan;
            TPlot(TPlot>273)=273; %%%%%new
            TPlot(TPlot<94)=94; %%%%%new
            axs{i,j} = subplot(4,5,j+(i-1)*10);
            hold on
            contourf(X*d,Y*d,TPlot,40,'linestyle','none');
            caxis([94 273]);
            axis equal
            xlim([0 20]);
            ylim([0 60]);
            axs{i,j}.XTickLabel = [];
            axs{i,j}.YTickLabel = [];
            axs{i,j}.Units = 'centimeters';
            colormap(axs{i,j},'bone');
    
            if i == 1 && j == 5
                c = colorbar('EastOutside');
                c.Label.Interpreter = 'latex';
                c.TickLabelInterpreter = 'latex';
                c.Label.FontSize = fs;
                c.Label.String = 'Temperature, K';
                if i == 1
                    cs{1} = c;
                else
                    cs{2} = c;
                end
             end
    
            axs{i,5+j} = subplot(4,5,(i-1)*10+5+j);
            axs{i,5+j}.Colormap = colormap(flipud(parula));
            contourf(X*d,Y*d,phiPlot,40,'linestyle','none');
            axis equal
            xlim([0 20]);
            ylim([0 60]);
            axs{i,5+j}.XTickLabel = [];
            axs{i,5+j}.YTickLabel = [];
            axs{i,5+j}.Units = 'centimeters';
            caxis([0 1]);
           
            
            if i == 1 && j == 5
                c = colorbar('EastOutside');
                c.Label.Interpreter = 'latex';
                c.TickLabelInterpreter = 'latex';
                c.Label.FontSize = fs;%8;
                c.Label.String = 'Melt fraction';
                if i == 1
                    cs{3} = c;
                else
                    cs{4} = c;
                end
            end
            tFrame(i,j) = round(tVec(end),2);
            
        end
    end
    
    
    %% rearrange
    for i = 1:2
        for j = 1:5
            if i == 1
                axs{i,j}.Position = [lMarg+(j-1)*(wid+hMarg),bMarg+3*(wid+vMarg1)+vMarg2 , wid, wid];
            elseif i == 2
                axs{i,j}.Position = [lMarg+(j-1)*(wid+hMarg),bMarg+1*(wid+vMarg1), wid, wid];
            end
            
            if i == 1
                axs{i,5+j}.Position = [lMarg+(j-1)*(wid+hMarg),bMarg+2*(wid+vMarg1)+vMarg2 , wid, wid];
            elseif i == 2
                axs{i,5+j}.Position = [lMarg+(j-1)*(wid+hMarg),bMarg, wid, wid];
            end
        end
    end
            
    yTikLab = {'0','20','40','60'};
    xTikLab = {'0','5','10','15','20'};
    tik1 = [0,10,20];
    tik2 = [0,20,40,60];
    
    
    % l b w h
    i = 2;
    j = 6;
    % axs{i,j}.Position = [lMarg,bMarg,wid,wid];
    axs{i,j}.YTick = tik2;
    axs{i,j}.YTickLabel = yTikLab;
    axs{i,j}.XTick = tik1;
    axs{i,j}.XTickLabel = xTikLab;
    axs{i,j}.YLabel.String = 'z-dir [km]';
    axs{i,j}.XLabel.String = 'r [km]';
    xtickangle(axs{i,j}, angle); axs{i,j}.TickDir = 'in';   axs{i,j}.TickLength = [0.02 0.02];     % ticks pointing inward % make them longer axs{i,j}.TickDir = 'in';      % ticks pointing inward axs{i,j}.TickLength = [0.02 0.02];  % make them longer
    
    axs{i,j}.YLabel.FontSize = fs;
    axs{i,j}.XLabel.FontSize = fs;
    
    i = 2;
    j = 1;
    % axs{i,j}.Position = [lMarg,bMarg+wid+vMarg1,wid,wid];
    axs{i,j}.YTick = tik2;
    axs{i,j}.YTickLabel = yTikLab;
    axs{i,j}.YLabel.String = 'z-dir [km]';
    
    axs{i,j}.YLabel.FontSize = fs;
    xtickangle(axs{i,j}, angle); axs{i,j}.TickDir = 'in';   axs{i,j}.TickLength = [0.02 0.02];     % ticks pointing inward % make them longer
    
    i = 1;
    j = 1;
    % axs{i,j}.Position = [lMarg,bMarg+2*(wid+vMarg1)+vMarg2,wid,wid];
    axs{i,j}.YTick = tik2;
    axs{i,j}.YTickLabel = yTikLab;
    axs{i,j}.YLabel.String = 'z-dir [km]';
    axs{i,j}.XTickLabel = [];
    axs{i,j}.YLabel.FontSize = fs;
    xtickangle(axs{i,j}, angle); axs{i,j}.TickDir = 'in';   axs{i,j}.TickLength = [0.02 0.02];     % ticks pointing inward % make them longer
    
    i = 1;
    j = 6;
    % axs{i,j}.Position = [lMarg,bMarg+3*(wid+vMarg1)+vMarg2,wid,wid];
    axs{i,j}.YTick = tik2;
    axs{i,j}.YTickLabel = yTikLab;
    axs{i,j}.YLabel.String = 'z-dir [km]';
    axs{i,j}.XTickLabel = [];
    axs{i,j}.YLabel.FontSize = fs;
    xtickangle(axs{i,j}, angle); axs{i,j}.TickDir = 'in';   axs{i,j}.TickLength = [0.02 0.02];     % ticks pointing inward % make them longer
    
    %colorbar loc
    for i = 1:length(cs)
        cs{i}.Position(3) = 0.01;
        cs{i}.Position(1) = 0.85;
    end
    %%
    figLab = 'abcdefghij';
    for i = 1:5
        colormap(axs{1,i},'bone'); 
        colormap(axs{2,i},'bone'); 
        set([axs{2,i},axs{2,i+5}],'Position',axs{1,i+5}.Position)
        axs{1,i+5}.Visible = 'off';
        axs{2,i+5}.Visible = 'off';
        set([axs{1,i},axs{1,i+5}],'Position',axs{1,i}.Position)
        if caption == 1
        text(axs{2,i+5},19,0,[num2str(tFrame(2,i)) ' yrs'],...
            'FontSize',fs_small,'HorizontalAlignment','right','Color',[1,0,0])
        text(axs{1,i+5},19,0,[num2str(tFrame(1,i)) ' yrs'],...
            'FontSize',fs_small,'HorizontalAlignment','right','Color',[1,0,0])
        text(axs{1,i+5},19,9.9,figLab(i),'fontweight','bold','interpreter','tex',...
            'FontSize',fs_small,'HorizontalAlignment','right','Color',[1,0,0])
        text(axs{2,i+5},19,9.9,figLab(i+5),'fontweight','bold','interpreter','tex',...
            'FontSize',fs_small,'HorizontalAlignment','right','Color',[1,0,0])
        end
    end
    
    
    i = 2;
    j = 1;
    axs{i,j}.XTick = tik1;
    axs{i,j}.XTickLabel = xTikLab;
    axs{i,j}.XLabel.String = 'r [km]';
    axs{i,j}.XLabel.FontSize = fs;
    xtickangle(axs{i,j}, angle); axs{i,j}.TickDir = 'in';   axs{i,j}.TickLength = [0.02 0.02];     % ticks pointing inward % make them longer
    
    i = 2;
    j = 2;
    axs{i,j}.XTick = tik1;
    axs{i,j}.YTick = tik2;
    axs{i,j}.XTickLabel = xTikLab;
    %axs{i,j}.XLabel.String = 'r [km]';
    axs{i,j}.XLabel.FontSize = fs;
    xtickangle(axs{i,j}, angle); axs{i,j}.TickDir = 'in';   axs{i,j}.TickLength = [0.02 0.02];     % ticks pointing inward % make them longer
    axs{i,j}.XTickLabel = [];
    
    i = 2;
    j = 3;
    axs{i,j}.XTick = tik1;
    axs{i,j}.YTick = tik2;
    axs{i,j}.XTickLabel = xTikLab;
    %axs{i,j}.XLabel.String = 'r [km]';
    axs{i,j}.XLabel.FontSize = fs;
    xtickangle(axs{i,j}, angle); axs{i,j}.TickDir = 'in';   axs{i,j}.TickLength = [0.02 0.02];     % ticks pointing inward % make them longer
    axs{i,j}.XTickLabel = [];
    
    i = 2;
    j = 4;
    axs{i,j}.XTick = tik1;
    axs{i,j}.YTick = tik2;
    axs{i,j}.XTickLabel = xTikLab;
    %axs{i,j}.XLabel.String = 'r [km]';
    axs{i,j}.XLabel.FontSize = fs;
    xtickangle(axs{i,j}, angle); axs{i,j}.TickDir = 'in';   axs{i,j}.TickLength = [0.02 0.02];     % ticks pointing inward % make them longer
    axs{i,j}.XTickLabel = [];
    
    i = 2;
    j = 5;
    axs{i,j}.XTick = tik1;
    axs{i,j}.YTick = tik2;
    axs{i,j}.XTickLabel = xTikLab;
    %axs{i,j}.XLabel.String = 'r [km]';
    axs{i,j}.XLabel.FontSize = fs;
    xtickangle(axs{i,j}, angle); axs{i,j}.TickDir = 'in';   axs{i,j}.TickLength = [0.02 0.02];     % ticks pointing inward % make them longer
    axs{i,j}.XTickLabel = [];
    
    i = 1;
    j = 2;
    axs{i,j}.XTick = tik1;
    axs{i,j}.YTick = tik2;
    xtickangle(axs{i,j}, angle); axs{i,j}.TickDir = 'in';   axs{i,j}.TickLength = [0.02 0.02];     % ticks pointing inward % make them longer
    
    i = 1;
    j = 3;
    axs{i,j}.XTick = tik1;
    axs{i,j}.YTick = tik2;
    xtickangle(axs{i,j}, angle); axs{i,j}.TickDir = 'in';   axs{i,j}.TickLength = [0.02 0.02];     % ticks pointing inward % make them longer
    
    i = 1;
    j = 4;
    axs{i,j}.XTick = tik1;
    axs{i,j}.YTick = tik2;
    xtickangle(axs{i,j}, angle); axs{i,j}.TickDir = 'in';   axs{i,j}.TickLength = [0.02 0.02];     % ticks pointing inward % make them longer
    
    i = 1;
    j = 5;
    axs{i,j}.XTick = tik1;
    axs{i,j}.YTick = tik2;
    xtickangle(axs{i,j}, angle); axs{i,j}.TickDir = 'in';   axs{i,j}.TickLength = [0.02 0.02];     % ticks pointing inward % make them longer
    
    f.PaperPosition = [0.0177,4.6440,4.4646,5.7121];
    Filename = sprintf('FinalStokesvsStokesEuropa_vectorpart_%s',num2str(kk));
    set(gcf,'PaperPositionMode','auto'); 
    saveas(gcf,[Filename '.pdf']);
    
    Filename = sprintf('FinalStokesvsStokesEuropa_vectorpart_%s',num2str(kk));
    set(gcf,'PaperPositionMode','auto');
    print(gcf, Filename, '-dpdf', '-painters');
end

