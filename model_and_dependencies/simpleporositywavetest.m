function simpleporositywavetest(fn,eta_0,E_a)
    %{
    Function to evolve impact melt chambers on Europa. Simulations end
    when there is a negligible amount of melt left from the impact. The
    ouputs are saved at the breakpoint in line 294. This model is meant for
    short term simulations, on the order of the sinking of impact melts.
    
    Variables:
        fn (string): is the directory to find the initial conditions (the
        outputs from the impact simulation)i
        eta_0 (int): basal viscosity of ice in Pa s, 1e14 is a common value
        E_a (int): visocisty activation energy in Arrhenihus relationship
        given in J/mol, 50e3 is a common value
        no_of_tracers: Number of tracers = 0, 1, 2

    Author: Mohammad Afzal Shadab, mashadab@utexas.edu, December 25th, 2023 (search %%%%); 
            Evan Carnahan, evan.carnahan@utexas.edu, 11/20/2022
    %}

    set(groot,'defaultAxesFontName','Times')
    set(groot,'defaultAxesFontSize',20)
    set(groot,'defaulttextinterpreter','latex')
    set(groot,'defaultAxesTickLabelInterpreter','latex')
    set(groot,'defaultLegendInterpreter','latex')
    set(groot, 'DefaultFigureVisible', 'on');
    warning off; % matrix is close to singular due to viscosity contrast
    %% Load initial condition to be evolved
    % make ice shell thickness based on impact code passed from iSALE
    if fn == 'Wakita-Stokes'
        d = 20*1e3; % ice shell thickness, m
    end

    % threshold of initial fluid left in ice shell to stop simulation at
    termFrac = 0.005;
   
    %Initial volume
    
    %% Set physical parameters and make dimensionless scales
    % physical parameters for ice
    T_t = 100.0; % surface temperature, K
    a = 185; % ice specific heat parameter
    b = 7.037; % ice specific heat parameter
    T_b = 273.16; % melting temperature, K
    rho_i = 917; % ice density at melting temperature, kg/m^3
    grav = 9.81; % gravity: Europa =1.315, m/s^2 , Titan = 1.352, m/s^2
    DT = T_b - T_t; % difference in temperature
    alpha = 1.6e-4; % thermal expansion coefficient, 1/K
    R = 8.314; % universal gas constant, J K^-1 mol^-1
    Apar = E_a/R/T_b; % viscosity exponenet

    % physical properties of water
    latHeat = 334e3; % latent heat of fusion, J/kg
    rho_w = 1e3; % density of water, kg/m^3
    c_pw = 4.186e3; % specific ehat of water, J/kg
    kappa_w = 0.56; %thermal diffusivity of water, W/(m K)
    porViscPar = 45;%How much a certain amt of melt reduce viscosity of ice  %%%% effective solid visc
    
    % temperature and melt fraction dependent viscosity, Pa s
    %barrViscPhi = @(nonT,phi) max(exp(Apar*(T_b./(DT.*nonT+T_t)-1)).*exp(-porViscPar*phi),1e-2);
    %barrViscPhi = @(nonT,phi) max(exp(Apar*(T_b./(DT.*nonT+T_t)-1)).*exp(-porViscPar*phi).*(1-phi),1e-5); %%%%New viscosity of solid (1-phi)*mu_s
    %barrViscPhi = @(nonT,phi) max(exp(Apar*(T_b./(DT.*nonT+T_t)-1)).*exp(-porViscPar*phi),1e-2); %Old viscosity
    barrViscPhi = @(nonT,phi) max(exp(Apar*(T_b./(DT.*nonT+T_t)-1)).*exp(-porViscPar*phi),1e-2).*(1-phi); %Old viscosity with water softening
    % viscosity is max of Temp dependence on viscosity, melt dependence, threshold
    % threshold 1e-2 means two orders of magnitude reduction is essentially inviscid
    c_fun = @(nonT) a+b*(DT*nonT+T_t); %specific heat function, J/(kg K)
    kappa_b = 3.3; %thermal conductivity of ice, set to be consistent with Cox and Bauer, 2015
    D_T = kappa_b/(rho_i *c_fun(1)); % thermal diffusicvity of ice, m^2/s
    c_pi = c_fun(1); %constant specific heat, J/(kg K)
    
    % thershold for melting: dimensionless boundary between partial and
    % fully melted regions
    mixZone = (rho_w*latHeat)/(rho_i*c_pi*DT); %fully ice to fully water mix spectrum
                                               %similar to Stefan number = latent heat / sensible heat
    
    phi_fun = @(nonH) nonH * (rho_i*c_pi*DT)/(rho_w*latHeat); %porosity: fully ice = 0 to fully water = 1
    TWater_fun = @(nonH) nonH * (rho_i*c_pi)/(rho_w*c_pw) - latHeat/(DT*c_pw) + 1; %Dim Temperature, K
    compBouy_fun = @(phi)  phi*rho_i*(rho_w/rho_i-1)*(grav*d^3/(eta_0*D_T)); %Dimless RHS of mom balance
                                                                             %Ice is lighter than water
    
    % conversion in dimensionless units with constant specific heat
    nonH_fun = @(nonT) nonT - 1; 
    nonT_fun = @(nonH) nonH + 1;
    % condictivity is weighted average of mixture components
    porKappaPrime_fun = @(phi,nonT) (phi*kappa_w + (1-phi).*kappa_b)/kappa_b; %thermal cond. W/(m-K)
    porNonH_fun = @(phi,nonT) (1-phi).*(nonH_fun(nonT)) + ...
        (phi*rho_w)./(rho_i*c_pi*DT).*(latHeat+c_pw*DT*(nonT-1)); %Dimless enthalpy
    
    % characteristic scales for general convection
    t_c = d^2/D_T; % convection time scale, s
    Ra = rho_i*grav*alpha*d^3*DT/(eta_0*D_T); % basal Rayleigh number
    
    %%%%%%%%%
    kc = 5.6e-11; %Absolute permeability [in m^2] (From Meyer and Hewitt (2017)) 
    mu_f = 1e-3;  %Viscosity of water phase [in Pa.s] (Duh)
    rho_f = 1e3;  %Density of water phase [in kg/m^3] (Duh) 
    
    
    %Compaction length
    GG = 1; %G from compaction viscosity relation
    mm = 1; %m from compaction viscosity relation
    nn  = 2; %porosity-permeability relation
    phic = 1e-2; mu_max = eta_0 * exp(-porViscPar*phic); %Char. porosity and viscosity
    delta0 = sqrt(kc*phic^nn*mu_max/(phic^mm*mu_f)) %Compaction length
    Delta_rho = rho_w - rho_i;
    Kc = kc*Delta_rho*grav*phic^nn/mu_f; %Char. hydraulic conductivity
    HD = 1; %Dimensionless depth to compaction length ratio
    
    
    Pi_1 = eta_0 * kc/(mu_f * d^2);      %For updated hydraulic condcutivity of total mass balance
    Pi_3 = kc*d*rho_f*grav/(mu_f * D_T); %For RHS of mass balance (fluid diffusion over heat diffusion)
    Pi_5 = d^3 * rho_i * grav / (eta_0 * D_T);      %Pi 5 from the notes
    Pi_6 = mixZone;                      %Pi 6 from the notes
    %%%%%%%%%

    % thermal condutvity in convective ocean, set to maintain vertical
    % isotherm in ocean
    kappa_c = 100; %thermal cond. W/(m-K) is set high to make ocean isothermal

    %% build cylindirical grid for numerical solution
    % build grid
    Gridp.xmin = 0; Gridp.xmax = HD*delta0/d; Gridp.Nx = 5; %radial direction
    Gridp.ymin = 0; Gridp.ymax = HD*delta0/d; Gridp.Ny = 100; %vertical direction
    Gridp.geom = 'cartesian';       %geometry type: cylinderical r-z coordinates
    Grid = build_stokes_grid(Gridp); %build grid for Stokes equation in 
    [X,Y]= meshgrid(Grid.p.xc,Grid.p.yc);
    
    % convert inital condition to grid
    TGr = zeros(size(X));%reshape(T,grRes,Grid.p.Nx);     %Temp on the grid, K
    phiGr = phic*ones(size(X));%reshape(phi,grRes,Grid.p.Nx); %porosity on the grid
    
    T = TGr(:);
    phi = phiGr(:);
    
    
    
    %Analytic solution
    zDa = linspace(0,HD,1e3);
    % coefficients
    c1 = @(H) (exp(-H)-1)./(exp(H)-exp(-H));
    c2 = @(H) (exp(H)-1)./(exp(H)-exp(-H));
    % potentials
    hDa = @(z,H) z + c1(H).*exp(z) + c2(H).*exp(-z);
    uDa = @(z,H) -z - c1(H).*exp(z) - c2(H).*exp(-z);
    % overpressure
    pDa = @(z,H) c1(H).*exp(z) + c2(H).*exp(-z);
    % flux & velocity
    qDa = @(z,H) -1 - c1(H).*exp(z) + c2(H).*exp(-z);
    vDa = @(z,H) 1 + c1(H).*exp(z) - c2(H).*exp(-z);
    
    H = porNonH_fun(phi,T); %Dimensionless porosity array

    figure();
    contourf(X,Y,phiGr);
    
    figure();
    contourf(X,Y,TGr);   
    
    
    %% build operators
    Zp = sparse(Grid.p.N,Grid.p.N); %N by N Zero matrix
    Ip = speye(Grid.p.N); %N by N indentity matrix
    [D,Edot,Dp,Gp,I,Gyy]=build_stokes_ops(Grid); %Making stokes operators
    
    linInds = find(Gyy > 0);
    [row,~] = ind2sub(size(Gyy),linInds);
    
    %% Build boundary conditions for temperature and flow equation
    % Fix temperature at top of ice shell with Dirchlet BC
    pc = (mu_max * D_T)/d^2 ; %Characteristic pressure
    
    % Free slip boundary condition for Stokes equation
    %Dirichlet BC below; Natural BCs are automatically defined
    Param(1).dof_dir = [Grid.dof_ymax_vt(2:end-1);...  % tangential velocity on the top
              Grid.dof_pene;...     % no penetration on all bnd's
              Grid.dof_ymin_vt(2:end-1);...   
              Grid.dof_pc_comp_col];         % pressure constraint
    Param(1).dof_f_dir = [];         
    Param(1).g = [zeros(length(Grid.dof_ymax_vt(2:end-1)),1);...       % tangential velocity on the top
              zeros(Grid.N_pene,1);...      % no penetration on all bnd's
              zeros(length(Grid.dof_ymin_vt(2:end-1)),1);...  
              hDa(Grid.p.dy*d/(2*delta0),HD)*delta0 * Delta_rho * grav/pc];                         % pressure constraint
    [B,N,fn] = build_bnd(Param(1),Grid,I);
    fs_T = nan(size(T));
    %{
    Param(1).dof_dir =  [...
                      Grid.x.dof_xmax;...           %set x_max x-vel
                      Grid.x.dof_xmin;...           %set x_min x-vel
                      Grid.x.N+Grid.y.dof_ymin;...  %set y_min y-vel
                      Grid.x.N+Grid.y.dof_ymax;...  %set y_max y-vel
                      Grid.p.Nf+1];                 %set pressure
    Param(1).g = 0*Param.dof_dir;
    Param.g(end) = 0;
    B = I([Param.dof_dir],:);
    N = I;
    N(:,[Param.dof_dir]) = [];
    fs_T = nan(size(T));
    %}
    %%%% 
    
    % create arrays for time evolution storage
    it = 1e9;
    %% temporal evolution
    for i = 1:1e9
        
        % calculate porosity from enthalpy
        %[T,phi] = enthMelt(H,mixZone,nonT_fun,phi_fun,TWater_fun);
        T = ones(size(T));
        phi = phic*ones(size(T));

        %% Bouyancy force as the body force to Stokes flow
        %temp bouyancy
        Tplot= reshape(T,Grid.p.Ny,Grid.p.Nx);
        Tdiag = comp_mean(Tplot,1,1,Grid.p);
        Tvec = diag(Tdiag);
        Ty = Tvec(Grid.p.Nfx+1:Grid.p.Nf);
        fs_T(Ty>1) = -Ra*1;         %All melt only has maximum RHS of Rayleigh number (T>Tm)
        fs_T(Ty<=1) = -Ra*Ty(Ty<=1);%Cold ice without melt has a lower magnitude Ra than 1 (T<=Tm)
        
        %compositional bouyancy
        phiPlot= reshape(phi,Grid.p.Ny,Grid.p.Nx);
        phiDiag = diag(comp_mean(phiPlot,1,1,Grid.p));
        phiY = phiDiag(Grid.p.Nfx+1:Grid.p.Nf);
        fs_por = compBouy_fun(phiY); %calculating RHS of momentum balance
        
        %Gxx variable viscosity matrix
        nxxVec = zeros(Grid.x.Nfx,1);
        nxxVec(Grid.x.Ny+Grid.p.dof) = Tplot;
        
        %Gyy variable viscosity matrix
        nyyVec = zeros(Grid.y.Nfy,1);
        nyyVec(row) = Tplot;
        ncVec = comp_mean_corners(Tplot,-1,Grid.p);
        tempVec = [nxxVec;nyyVec;ncVec];
        
        %% Porosity and temperature dependent viscosity
        % porosity dependent viscosity
        phiPlot = reshape(phi,Grid.p.Ny,Grid.p.Nx);
        nxxVecPhi = zeros(Grid.x.Nfx,1);
        nxxVecPhi(Grid.x.Ny+Grid.p.dof) = phiPlot;
        
        %Gyy variable viscosity matrix
        nyyVecPhi = zeros(Grid.y.Nfy,1);
        nyyVecPhi(row) = phiPlot;
        
        ncVecPhi = comp_mean_corners(phiPlot,-1,Grid.p);
        phiVec = [nxxVecPhi;nyyVecPhi;ncVecPhi];
        
        % merge temperature and melt fraction viscosities
        viscVec = barrViscPhase(tempVec,phiVec,barrViscPhi); %%%% effective solid visc;
        viscVec(isnan(viscVec)) = 0; %Due to mean, 1/0 appear at corners so making them 0
        viscMat = spdiags(viscVec,0,length(viscVec),length(viscVec));
        
        %%%%
        Zd_vec =  GG * (1-phi)./(phi+1e-5).^mm - 2/3*(1-phi); %dim-less compaction viscosity

        %size(Zd_vec) 
        %size(comp_mean_corners(Zd_vec,-1,Grid.p))
        Zd = spdiags(Zd_vec,0,Grid.p.N,Grid.p.N);   %transforming into a matrix
        
        Kdprime_vec = comp_mean(phiPlot.^nn,-1,1,Grid.p);
        Kdprime = spdiags(Kdprime_vec,0,length(Kdprime_vec),length(Kdprime_vec)); %dimensional permeability

        fs_fluid_pressure =  - Pi_5 * (rho_f/rho_i) * ones(Grid.p.Nfy,1);  %Extra source term from Pf
        
        % higher porosity acts against Ra bouyancy force 
        fsVec = fs_T + fs_por + fs_fluid_pressure; %adding diffusion and convection equation RHS
        fsVecalter = -(Delta_rho*grav*(diag(comp_mean(1-phiGr,1,1,Grid.p)))*(Pi_5/(rho_i*grav))).* [zeros(Grid.p.Nfx,1); ones(Grid.p.Nfy,1)] ;
 
        
        min(fs_fluid_pressure*(mu_max*D_T)/d^3)
        min(fsVec*(mu_max*D_T)/d^3)
        min(fsVecalter*(mu_max*D_T)/d^3)
        min((fs_T + fs_por)*(mu_max*D_T)/d^3)
        
        %fs_mass = Dp * Kdprime * Pi_3 * [zeros(Grid.p.Nfx,1); ones(Grid.p.Nfy,1)];  %%%%source term of total mass balance; no Gamma term
        fs_mass  = zeros(Grid.p.N,1);   %No source term in mass balance %%%%
        fs = [zeros(Grid.p.Nfx,1); fsVec; fs_mass];  %%%%Added RHS term for total mass balance
        %fs = [fsVec; fs_mass];  %%%%Added RHS term for total mass balance new
        
        %%%%

        %% Stokes Flow calcualtion
        % make linear operators
        tau = D*2*viscMat*Edot;
         %L = [tau,Gp;                             %Two-phase slurry model
         %    Dp,Zp];
        L = [tau+ Gp * Zd * Dp,    -Gp;                 %Darcy-Stokes
             Dp,                   -Dp * Pi_1* Kdprime * Gp];
        % solve for flow velocities
        u = solve_lbvp(L,fs,B,Param.g,N);
        vx = u(1:Grid.p.Nfx);
        vy = u(Grid.p.Nfx+1:(Grid.p.Nfx+Grid.p.Nfy));
        vm = [vx;vy];       %velocity of the solid phase, m/s
        vmax= max(vm); %largest solid velocity
        p  = u((Grid.p.Nfx+Grid.p.Nfy+1):end); %Dimless total fluid pressure coupled with gravitational head
        %vf = vm - Pi_1 * comp_mean(phiPlot.^(nn-1),1,1,Grid.p) * ( Gp * p + (rho_w/rho_i) * Pi_5 * [zeros(Grid.p.Nfx,1); ones(Grid.p.Nfy,1)]); %calculate the dimless fluid velocity %%%%
        vf = vm - Pi_1 * comp_mean(phiPlot.^(nn-1),1,1,Grid.p) * ( Gp * p ); %calculate the dimless fluid velocity %%%%
        p_overall = p;
        p  = p - (rho_f/rho_i) * Pi_5 * Y(:);  %Calculating fluid from overall pressure  %%%%
        vfmax= max(vf); %largest solid velocity        
        % Adaptive time stepping based on competition b/w CFL and Neumann
        % conditions in each direction; CFL number set to 0.8
        
        pL = -rho_i * grav *Grid.p.yc*d; %Lithostatic or solid pressure [Pa]
        
        
        hhh=figure()
        set(gcf,'units','points','position',[0,0,3125,1250])
        % Enlarge figure to full screen.
        set(gcf, 'Position', [50 50 1500 600])
        t=sgtitle(sprintf('Dim-less depth = %.1f, Porosity = %0.3f',HD, phic)); t.FontSize = 20;
        %% Plotting and post-processing
        subplot 131
        plot(hDa(zDa,HD)*delta0,zDa*delta0,'linewidth',2), hold on
        plot((p_overall(1:Grid.p.Ny)*pc)/(Delta_rho*grav),Grid.p.yc*d,'--','linewidth',2)
        xlabel('h [m]','fontsize',22)
        ylabel('z [m]','fontsize',22)
        legend('analytic','numerical','location','northwest')
        set(gca,'fontsize',18)

        subplot 132
        plot(pDa(zDa,HD)*Delta_rho*grav*delta0,zDa*delta0,'linewidth',2), hold on
        plot(p_overall(1:Grid.p.Ny)*pc - Delta_rho * grav * Grid.p.yc*d,Grid.p.yc*d,'--','linewidth',2)
        xlabel('Overpressure [Pa]','fontsize',22)
        ylabel('z [m]','fontsize',22)
        legend('analytic','numerical','location','northeast')
        set(gca,'fontsize',18)

        subplot 133
        plot(vDa(zDa,HD)*Kc,zDa*delta0,'linewidth',2), hold on
        plot(vm(Grid.p.Nfx+1:Grid.p.Nfx+1+Grid.p.Ny),Grid.p.yf*d,'--','linewidth',2)
        xlabel('Solid velocity [m/s]','fontsize',22)
        ylabel('z[m]','fontsize',22)
        legend('analytic','numerical','location','northeast')
        set(gca,'fontsize',18)
        saveas(hhh,sprintf('instant_comp_column_phic%d_HD%d.png',phic,HD));

        
        
        
        
        break %end simulation
        
  
    end
%{
%%%%
%% Making a video out of frames
 % create the video writer with fps of the original video
 Data_result= sprintf('../figures/case%s_t%syrs_kc%sTlayer%sC..avi',num2str(fn),num2str(tTot),num2str(kc),num2str(Tlayer));
 writerObj = VideoWriter(Data_result);
 writerObj.FrameRate = 20; % set the seconds per image
 open(writerObj); % open the video writer
% write the frames to the video
for i=1:frameno
    %'Frame number'; i
    % convert the image to a frame
    frameimg = FF(i) ;
    writeVideo(writerObj, frameimg);
end
% close the writer object
close(writerObj);    
    %%%%
%}  
end