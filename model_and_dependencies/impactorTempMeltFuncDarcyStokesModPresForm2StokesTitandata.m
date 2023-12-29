function impactorTempMeltFuncDarcyStokesModPresForm2StokesTitandata(fn,eta_0,E_a)
    %{
    Function to evolve impact melt chambers on Europa. Simulations end
    when there is a negligible amount of melt left from the impact. The
    ouputs are saved at the breakpoint in line 294. This model is meant for
    short term simulations, on the order of the sinking of impact melts.
    
    Variables:
        fn (string): is the directory to find the initial conditions (the
        outputs from the impact simulation)
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
    if fn == 'Wakita'
        d = 20*1e3; % ice shell thickness, m
    end

    % threshold of initial fluid left in ice shell to stop simulation at
    termFrac = 0.005;
    
    % load simulations from initial conditions folder (Digitized output via
    % Python from iSALE simulations)
    fp = '../initial_conditions/'; %loading the initial conditions
    load([fp 'Wakita_data.mat']); %loading porosity phi and temp T

    T   = flipud(Tdata);
    phi = flipud(phidata);
    
    %% Set physical parameters and make dimensionless scales
    % physical parameters for ice
    T_t = 100; % surface temperature, K
    a = 185; % ice specific heat parameter
    b = 7.037; % ice specific heat parameter
    T_b = 273; % melting temperature, K
    rho_i = 920; % ice density at melting temperature, kg/m^3
    grav = 1.352; % gravity: Europa =1.315, m/s^2 , Titan = 1.352, m/s^2
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
    barrViscPhi = @(nonT,phi) max(exp(Apar*(T_b./(DT.*nonT+T_t)-1)).*exp(-porViscPar*phi),1e-2); %Old viscosity
    %barrViscPhi = @(nonT,phi) max(exp(Apar*(T_b./(DT.*nonT+T_t)-1)).*exp(-porViscPar*phi).*(1-phi),1e-2); %Old viscosity with water softening
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
    
    
    Pi_1 = eta_0 * kc/(mu_f * d^2);      %For updated hydraulic condcutivity of total mass balance
    Pi_3 = kc*d*rho_f*grav/(mu_f * D_T); %For RHS of mass balance (fluid diffusion over heat diffusion)
    Pi_5 = d^3 * rho_i * grav / (eta_0 * D_T);      %Pi 5 from the notes
    Pi_6 = mixZone;                      %Pi 6 from the notes
    %%%%%%%%%
    
    % non-dimensionalize temperature
    T = (T - T_t)/DT;

    % thermal condutvity in convective ocean, set to maintain vertical
    % isotherm in ocean
    kappa_c = 100; %thermal cond. W/(m-K) is set high to make ocean isothermal

    %% build cylindirical grid for numerical solution
    % build grid
    grRes = size(Ydata,2); % grid resolution in radial direction
    grZ = size(Ydata,1); % grid resolution in radial direction
    ocTh = floor(3*grZ/5); % make ocean below the ice shell: 1/5 factor being 20% of the ice shell
    Gridp.xmin = 0; Gridp.xmax = 1.5; Gridp.Nx = grRes; %radial direction
    Gridp.ymin = -3/5; Gridp.ymax = 3; Gridp.Ny = 3*grZ+ocTh; %vertical direction
    Gridp.geom = 'cylindrical_rz';       %geometry type: cylinderical r-z coordinates
    Grid = build_stokes_grid_cyl(Gridp); %build grid for Stokes equation in 
    [X,Y]= meshgrid(Grid.p.xc,Grid.p.yc);
    
    % convert inital condition to grid
    TGr = T;%reshape(T,grRes,Grid.p.Nx);     %Temp on the grid, K
    phiGr = phi;%reshape(phi,grRes,Grid.p.Nx); %porosity on the grid

    %%%%
    %vertical variation
    %TGr_layer = kron(linspace(1,T(2,end),2*grZ)',ones(1,Grid.p.Nx));
    
    % Adding a layer of temperate ice (fixed temperature)
    TGr_layer   = ones(2*grZ,Grid.p.Nx);%T(2,end)*ones(2*grZ,Grid.p.Nx);  %melting temperature in the ice %%%%
    phiGr_layer = zeros(2*grZ,Grid.p.Nx);
    %%%%
    
    % get initial melt volumes: phiOrig, m^3
    phiOrig = sum(sum(phiGr(10:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)
    % First sum is for porosity in the z direction since volume is same
    % Second sum is after multiplying the grid volume at bottom cells
    % Lastly d^3 comes from the redimensionalization to calculate melt volumes
    % First 10 rows are ignored since ocean there is neglected so 
    % phiGr(10:end,:) and not phiGr(1:end,:)
   
    % build ocean layer
    TOc = ones(ocTh,Grid.p.Nx);    %Temperature of ocean, K
    phiOc = ones(ocTh,Grid.p.Nx);  %Porosity of ocean, K 
    
    % combine ice shell and ocean to get the entire domain fields
    TGr = [TOc; TGr_layer;TGr]; %Temperature, K
    phiGr = [phiOc;phiGr_layer; phiGr]; %Porosity or melt fraction
    T = TGr(:); %Making Temperature array from grid, N by 1
    phi = phiGr(:); %Making porosity array from grid, N by 1
    H = porNonH_fun(phi,T); %Dimensionless porosity array

    figure();
    contourf(X,Y,phiGr);
    
    figure();
    contourf(X,Y,TGr);   
    
    %%%%
    %initializing tracer with max conc. unity in a certain region
    %1. Organic tracer
    %dimensionless z = 0.8 to 1.0
    trc1 = ones(Grid.p.N,1); %initializing tracer with max conc. unity
    trc1(Y(:)<2.85) = 0;    %How thick is the organics layer? 1 - 0.975 dimensional units
    
    trc1(phiGr>0.05) = 1;     %Melted region being initialized
    trc1(Y(:)<0.05) = 0;      %Taking out the ocean
    
    %Clathrates tracer
    %dimensionless z = 0.8 to 1.0
    trc2 = ones(Grid.p.N,1); %initializing tracer with max conc. unity
    trc2(Y(:)<2.5) = 0;    %How thick is the clathrate layer? 1 - 0.975 dimensional units
    
    trc2(phiGr>0.1) = 1;   %Melted region being initialized
    trc2(Y(:)<0.5)  = 1; %second layer
    trc2(Y(:)<0.0)  = 0;   %Taking out the ocean
    
    %%%%
    
    %% build operators
    Zp = zeros(Grid.p.N); %N by N Zero matrix
    Ip = speye(Grid.p.N); %N by N indentity matrix
    [D,Edot,Dp,Gp,I,Gyy]=build_stokes_ops_cyl(Grid); %Making stokes operators
    % Refer to https://mhesse.github.io/numerical_modeling/spring2022/StokesOps.pdf
    % D: Stokes divergence, Edot = operator for strain rate = \nabla (v) + \nabla (v))^T
    % when operates on a velocity vector v, Dp: N by Nf divergence operator on pressure grid
    % Gp: Nf by N gradient operator on pressure grid
    % I is (Nf + N) by (Nf + N) identity matrix to apply BCs
    % Gyy: Gradient of y directional grid in y-direction
    
    linInds = find(Gyy > 0);
    [row,~] = ind2sub(size(Gyy),linInds);
    
    %% Build boundary conditions for temperature and flow equation
    % Fix temperature at top of ice shell with Dirchlet BC
    T0 = 0;
    H0 = nonH_fun(T0);
    Param = struct('H',{},'g',{},'dof_dir',{});
    Param(1).H = struct('dof_dir',{},'dof_f_dir',{},'g',{},'dof_neu',{},'dof_f_neu',{},'qb',{},'dof_out',{});
     
    % fix bottom heat flux, Neumann BC, to maintain linear geotherm in ice shell
    qPrime = 1;
    %%%%
    %{
    Param.H(1).dof_dir = [Grid.p.dof_ymax];
    Param.H(1).dof_f_dir = [Grid.p.dof_f_ymax];
    Param.H(1).g = [H0*ones(length(Grid.p.dof_ymax),1)];
    
    Param.H(1).dof_neu = [Grid.p.dof_xmin;Grid.p.dof_xmax;Grid.p.dof_ymin];
    Param.H(1).dof_f_neu = [Grid.p.dof_f_xmin;Grid.p.dof_f_xmax;Grid.p.dof_f_ymin];
    Param.H(1).qb = [0*Grid.p.dof_f_xmin;0*Grid.p.dof_f_xmax;qPrime*ones(size(Grid.p.dof_f_ymin))];
        %}
    
    Param.H(1).dof_dir = [Grid.p.dof_ymax,Grid.p.dof_ymin];
    Param.H(1).dof_f_dir = [Grid.p.dof_f_ymax,Grid.p.dof_f_ymin];
    Param.H(1).g = [H0*ones(length(Grid.p.dof_ymax),1);H(Grid.p.dof_ymin)];
    
    Param.H(1).dof_neu = [Grid.p.dof_xmin;Grid.p.dof_xmax];
    Param.H(1).dof_f_neu = [Grid.p.dof_f_xmin;Grid.p.dof_f_xmax];
    Param.H(1).qb = [0*Grid.p.dof_f_xmin;0*Grid.p.dof_f_xmax];    

    %%%%
    
    [BH,NH,fn_H] = build_bnd(Param.H,Grid.p,Ip);
    Param.H(1).dof_out = [Grid.p.dof_ymin];
    
    % Free slip boundary condition for Stokes equation
    %Dirichlet BC below; Natural BCs are automatically defined
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
    
    %%%% 
    %% Build boundary conditions for conservative tracer equation, c   
    Param(1).c = struct('dof_dir',{},'dof_f_dir',{},'g',{},'dof_neu',{},'dof_f_neu',{},'qb',{},'dof_out',{}); %%%%
    
    Param.c(1).dof_dir   =  [];
    Param.c(1).dof_f_dir =  [];
    Param.c(1).g         =  [];
    Param.c(1).dof_neu   =  [];
    Param.c(1).dof_f_neu =  [];
    Param.c(1).qb        =  [];
    [B_c,N_c,fn_c]    =  build_bnd(Param.c(1),Grid.p,Ip);
    %%%%
    
    % create arrays for time evolution storage
    it = 1e9;
    netMelt = [];
    phiDrain1 = 0;
    phiDrain2 = 0;
    tTot = 0;
    tVec = [];
    phiDrain1Vec = []; 
    phiDrain2Vec = [];
    phiFracRem = [];
    frameno = 0; %%%%Initializing frame number for plotting
    
    %% temporal evolution
    for i = 1:1e9
        
        % calculate porosity from enthalpy
        [T,phi] = enthMelt(H,mixZone,nonT_fun,phi_fun,TWater_fun);
        H = porNonH_fun(phi,T);
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
        GG = 1; %G from compaction viscosity relation
        mm = 1; %m from compaction viscosity relation
        Zd_vec =  GG * (1-phi)./(phi+1e-5).^mm - 2/3*(1-phi); %dim-less compaction viscosity

        %size(Zd_vec) 
        %size(comp_mean_corners(Zd_vec,-1,Grid.p))
        Zd = spdiags(Zd_vec,0,Grid.p.N,Grid.p.N);   %transforming into a matrix
        nn  = 3; %porosity-permeability relation
        
        Kdprime_vec = comp_mean(phiPlot.^nn,1,1,Grid.p);
        Kdprime = spdiags(Kdprime_vec,0,length(Kdprime_vec),length(Kdprime_vec)); %dimensional permeability

        fs_fluid_pressure =  - Pi_5 * (rho_f/rho_i) * ones(Grid.p.Nfy,1);  %Extra source term from Pf
        
        % higher porosity acts against Ra bouyancy force 
        fsVec = fs_T + fs_por + fs_fluid_pressure; %adding diffusion and convection equation RHS
        
        %fs_mass = Dp * Kdprime * Pi_3 * [zeros(Grid.p.Nfx,1); ones(Grid.p.Nfy,1)];  %%%%source term of total mass balance; no Gamma term
        fs_mass  = zeros(Grid.p.N,1);   %No source term in mass balance %%%%
        fs = [zeros(Grid.p.Nfx,1); fsVec; fs_mass];  %%%%Added RHS term for total mass balance
        
        %%%%

        %% Stokes Flow calcualtion
        % make linear operators
        tau = D*2*viscMat*Edot;
         %L = [tau,Gp;                             %Two-phase slurry model
         %    Dp,Zp];
        L = [tau+ Gp * Zd * Dp,    -Gp;                 %Darcy-Stokes
             Dp,                    -Dp * Pi_1* Kdprime * Gp];
        % solve for flow velocities
        u = solve_lbvp(L,fs,B,Param.g,N);
        vx = u(1:Grid.p.Nfx);
        vy = u(Grid.p.Nfx+1:(Grid.p.Nfx+Grid.p.Nfy));
        vm = [vx;vy];       %velocity of the solid phase, m/s
        vmax= max(abs(vm)); %largest solid velocity
        p  = u((Grid.p.Nfx+Grid.p.Nfy+1):end); %Dimless fluid pressure coupled with gravitational head
        %vf = vm - Pi_1 * comp_mean(phiPlot.^(nn-1),1,1,Grid.p) * ( Gp * p + (rho_w/rho_i) * Pi_5 * [zeros(Grid.p.Nfx,1); ones(Grid.p.Nfy,1)]); %calculate the dimless fluid velocity %%%%
        vf = vm - Pi_1 * comp_mean(phiPlot.^(nn-1),1,1,Grid.p) * ( Gp * p ); %calculate the dimless fluid velocity %%%%

        p  = p - (rho_f/rho_i) * Pi_5 * Y(:);  %Calculating fluid from overall pressure  %%%%
        vfmax= max(abs(vf)); %largest solid velocity        
        
        % Adaptive time stepping based on competition b/w CFL and Neumann
        % conditions in each direction; CFL number set to 0.8
        dt = min([0.5*Grid.p.dx^2/kappa_c, Grid.p.dx/vmax,Grid.p.dx/vfmax,Grid.p.dy/vfmax,0.5*Grid.p.dy^2/kappa_c, Grid.p.dy/vmax])*0.1;
        
        %% non-linear thermal conducitivity matricies
        kappaPrime = porKappaPrime_fun(phi,T); %thermal conductivity, K
        % select near boundary ocean cells
        ocLog = Y(:) < 2/grRes & phi > 0.5; %Porosity greater than 50% and z (2)
        kappaPrime(ocLog) = kappa_c;
        kappaPrimePlot = reshape(kappaPrime,Grid.p.Ny,Grid.p.Nx);
        kappaFace = comp_mean(kappaPrimePlot,1,1,Grid.p);
                
        %% Advection of enthalpy, diffusion of temperature
        %AH = build_adv_op(vm,H,dt,Gp,Grid.p,Param.H,'mc')*H; %Upwinding the enthalpy from center to faces
        AHs = build_adv_op(vm,(1-phi).*(T-1),dt,Gp,Grid.p,Param.H,'gd')*((1-phi).*(T-1)); %Upwinding the solid enthalpy from center to faces %%%%%%%
        AHf = build_adv_op(vf,phi*Pi_6,dt,Gp,Grid.p,Param.H,'gd')*(phi*Pi_6); %Upwinding the fluid enthalpy from center to faces %%%%%%%
        AH  = AHs + AHf; %Combine the advection operator
        
        %AHnew2 = build_adv_op(vm,((1-phi).*(T-1) + phi*Pi_6),dt,Gp,Grid.p,Param.H,'mc')*((1-phi).*(T-1) + phi*Pi_6);
        %AHnew3 = (build_adv_op(vm,(1-phi).*(T-1),dt,Gp,Grid.p,Param.H,'mc')*((1-phi).*(T-1))+build_adv_op(vm,phi*Pi_6,dt,Gp,Grid.p,Param.H,'mc')*(phi*Pi_6))
        %norm(H - ((1-phi).*(T-1) + phi*Pi_6))
        
        L_T_I = Ip;  % Implicit operator (Unity as the method is explicit) 
        L_T_E_T = - dt*(-Dp*kappaFace*Gp); %Linear operator of heat diffusion
        L_T_E_H = H - dt*(Dp*AH);         %Linear operator of heat advection
        RHS_T = L_T_E_H + L_T_E_T*T + (fn_H)*dt; %RHS of enthalpy balance
        H = solve_lbvp(L_T_I,RHS_T,BH,Param.H.g,NH); %time marching the enthalpy equation
        
        %% Advection of organic tracer fn_c is the same as natural BCs
        %Organics
        %Ac1 = build_adv_op(vm,trc1,dt,Gp,Grid.p,Param.c,'mc');%Upwinding the tracer conc. from center to faces
        %Ac1s = build_adv_op(vm,trc1.*(1-phi),dt,Gp,Grid.p,Param.c,'mc');%Upwinding the tracer conc. from center to faces (solid)  %%%%
        %Ac1f = build_adv_op(vf,trc1.*phi,dt,Gp,Grid.p,Param.c,'mc');%Upwinding the tracer conc. from center to faces (fluid) %%%%       
        %Ac1  = Ac1s + Ac1f; %Total transport flux of tracer 1 %%%%
        phi_face_vec = comp_mean(phiPlot,1,1,Grid.p)*ones(Grid.p.Nf,1);
        v_tot = vm.*(1-phi_face_vec) + vf.*phi_face_vec;  %%%% This is where we need to add more variation in velocity
        Ac1 = build_adv_op(v_tot,trc1,dt,Gp,Grid.p,Param.c,'gd');%Upwinding the tracer conc. from center to faces
        L_c_I = Ip;  % Implicit operator (Unity as the method is explicit) 
        L_c_E1 = Ip - dt*(Dp * Ac1);% Explicit operator of tracer advection
        RHS_c1 = L_c_E1 * trc1 + (fn_c) * dt; %Forming the vector B of Ax = B
        trc1 = solve_lbvp(L_c_I,RHS_c1, B_c, Param.c.g,N_c); %time marching the tracer equation
        
        %Clathrates
        %Ac2 = build_adv_op(vm,trc2,dt,Gp,Grid.p,Param.c,'mc');%Upwinding the tracer conc. from center to faces
        %Ac2s = build_adv_op(vm,trc2.*(1-phi),dt,Gp,Grid.p,Param.c,'mc');%Upwinding the tracer conc. from center to faces (solid)  %%%%
        %Ac2f = build_adv_op(vf,trc2.*phi,dt,Gp,Grid.p,Param.c,'mc');%Upwinding the tracer conc. from center to faces (fluid) %%%%       
        %Ac2  = Ac2s + Ac2f; %Total transport flux of tracer 1 %%%%
        Ac2 = build_adv_op(v_tot,trc2,dt,Gp,Grid.p,Param.c,'gd');%Upwinding the tracer conc. from center to faces
        L_c_E2 = Ip - dt*(Dp * Ac2);% Explicit operator of tracer advection
        RHS_c2 = L_c_E2 * trc2 + (fn_c) * dt; %Forming the vector B of Ax = B
        trc2 = solve_lbvp(L_c_I,RHS_c2, B_c, Param.c.g,N_c); %time marching the tracer equation
        
        %% calculate net melt and melt transported to "ocean"
        % make two planes to measure the melt transported through
        % plane 1
        abInd = 5;
        phiGr = reshape(phi,Grid.p.Ny,Grid.p.Nx);
        bdPhi = mean(phiGr(ocTh+abInd:ocTh+abInd+1,:),1).*Grid.p.V(Grid.p.dof_ymin)';
        vyGr = reshape(vy,Grid.y.Ny,Grid.y.Nx);
        bdVy = vyGr(ocTh+abInd+1,:);
        meltTrans = bdPhi .* bdVy / Grid.p.dy * dt * d^3;
        phiDrain2 = phiDrain2 - sum(meltTrans(bdVy < 0)); % m^3
        phiDrain2Vec = [phiDrain2Vec phiDrain2];
        
        % plane 2
        phiGr = reshape(phi,Grid.p.Ny,Grid.p.Nx);
        bdPhi = mean(phiGr(ocTh+(abInd-1):ocTh+(abInd-1)+1,:),1).*Grid.p.V(Grid.p.dof_ymin)';
        vyGr = reshape(vy,Grid.y.Ny,Grid.y.Nx);
        bdVy = vyGr(ocTh+(abInd-1)+1,:);
        meltTrans = bdPhi .* bdVy / Grid.p.dy * dt * d^3;
        phiDrain1 = phiDrain1 - sum(meltTrans(bdVy < 0)); % m^3
        phiDrain1Vec = [phiDrain1Vec phiDrain1];
        
        % total melt left
        netMelt = [netMelt sum(phi .* Grid.p.V * d^3)];
        % select near bouundary ocean cells
        % time stepping
        tTot = tTot + dt*t_c/(3.154e7); %Yrs
        tVec = [tVec tTot];
    
        % calculate remaining phi above near ocean, bottom 20% of ice shell
        phiRem = sum(sum(phiGr(ocTh+grRes/5:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3);
        phiFracRem = [phiFracRem phiRem/phiOrig];

        % condition for ending simulation
        if phiFracRem(end) < termFrac || (i > 1000 && phiFracRem(end) > phiFracRem(end-1)) || i >5000 %1500 to 5000
            % save point
            save(['impact_' fn '_eta0_' num2str(log10(eta_0)) '_Ea_' num2str(E_a/1e3) '_output.mat'],...
                'Tplot','phi','Grid','phiDrain1Vec','phiDrain2Vec','phiOrig','tVec',...
                'phiFracRem','T','phi','tVec','phiDrain1Vec','phiDrain2Vec','phiOrig')
            break
        end
    
        %% PLOTTING
         if mod(i,20) == 0
            greens = interp1([0;1],[1 1 1; 0.45, 0.65, 0.38],linspace(0,1,256));
            reds = interp1([0;1],[1 1 1;  190/255  30/255  45/255],linspace(0,1,256));
            blues = interp1([0;1],[1 1 1; 39/255  170/255  225/255],linspace(0,1,256));
            
             i
            %streamfunction plot
            h=figure('Visible', 'off'); %For visibility: h=figure(4);
            set(gcf,'units','points','position',[0,0,3125,1250])
            % Enlarge figure to full screen.
            [PSI,psi_min,psi_max] = comp_streamfun(vm,Grid.p);
            set(gcf, 'Position', [50 50 1500 600])
            ax1 = subplot(1,4,1);
            cla;
            hold on
            axis equal
            Tplot_dummy=Tplot;
            contourf(X*d/1e3,Y*d/1e3-Grid.p.dy,Tplot*DT+T_t,40,'linestyle','none'),view(2),hold on
            c1 = colorbar('NorthOutside');
            caxis([min(Tplot(:)) max(Tplot(:))]*DT+T_t);
            cStrVal = linspace(min(PSI(:)),max(PSI(:)),10);
            contour(Grid.p.xf*d/1e3,Grid.p.yf*d/1e3,PSI,'k','LevelList',cStrVal);
            c1.Label.String = 'Temperature, K';
            xlabel('x-dir, km');
            ylabel('z-dir, km');
            colormap(ax1,reds);

            %melt fraction plot
            ax2 = subplot(1,4,2);
            t = sgtitle(sprintf('time=%.3f years',tTot)); t.FontSize = 25;
            cla;
            axis equal
            hold on
            contourf(X,Y,reshape(phi,Grid.p.Ny,Grid.p.Nx),40,'linestyle','none'),view(2),hold on
            phi_dummy=phi;       
            contourf(X*d/1e3,Y*d/1e3-Grid.p.dy,reshape(phi_dummy,Grid.p.Ny,Grid.p.Nx),40,'linestyle','none'),view(2),hold on
            c2 = colorbar('NorthOutside');
            colormap(ax2,blues);
            xlabel('x-dir, km');
            ylabel('z-dir, km');
            c2.Label.String = 'Melt fraction, 1';
            
            %Organics tracer location concentration plot  
            % Create green-to-red colormap
            ax3 = subplot(1,4,3);
            cla;
            axis equal
            hold on
            c3 = colorbar('NorthOutside');
            %contour(X,Y,reshape(phi,Grid.p.Ny,Grid.p.Nx),'r','LevelList',5e-2),hold on
            trc1_plot = trc1; %trc1_plot(trc1<=1e-8) = nan;
            contourf(X*d/1e3,Y*d/1e3-Grid.p.dy,reshape(trc1_plot,Grid.p.Ny,Grid.p.Nx),40,'linestyle','none'),view(2)
            xlabel('x-dir, km');
            ylabel('z-dir, km');
            c3.Label.String = 'Organics Conc., 1';
            colormap(ax3,greens);
            
            %Clathrates tracer location concentration plot
            % Create green-to-red colormap
            ax4 = subplot(1,4,4);
            cla;
            axis equal
            hold on
            c4 = colorbar('NorthOutside');
            %contour(X,Y,reshape(phi,Grid.p.Ny,Grid.p.Nx),'r','LevelList',5e-2),hold on
            trc2_plot = trc2; %trc1_plot(trc1<=1e-8) = nan;
            contourf(X*d/1e3,Y*d/1e3-Grid.p.dy,reshape(trc2_plot,Grid.p.Ny,Grid.p.Nx),40,'linestyle','none'),view(2)
            xlabel('x-dir, km');
            ylabel('z-dir, km');
            c4.Label.String = 'Clathrates conc., 1';
            colormap(ax4,flipud(gray));
            
%             if i<1500
            saveas(h,sprintf('../figures/res_fig%d.png',i));            
%             end
            %{
            [PSI,psi_min,psi_max] = comp_streamfun(vm,Grid.p);
            set(gcf, 'Position', [50 50 1500 600])
            subplot(3,3,1)
            cla;
            hold on
            axis equal
            contourf(X*d/1e3,Y*d/1e3-Grid.p.dy,Tplot*DT+T_t,40,'linestyle','none'),view(2),hold on
            c = colorbar('NorthOutside');
            caxis([min(Tplot(:)) max(Tplot(:))]*DT+T_t);
            cStrVal = linspace(min(PSI(:)),max(PSI(:)),10);
            contour(Grid.p.xf*d/1e3,Grid.p.yf*d/1e3,PSI,'k','LevelList',cStrVal);
            c.Label.String = 'Temperature, K';
            xlabel('x-dir, km')
            ylabel('z-dir, km')

            %melt fraction plot
            subplot(3,3,2)
            sgtitle(sprintf('time=%.3f years',tTot));
            cla;
            axis equal
            hold on
            contourf(X,Y,reshape(phi,Grid.p.Ny,Grid.p.Nx),40,'linestyle','none'),view(2),hold on
            c = colorbar('NorthOutside');
            contour(X,Y,reshape(phi,Grid.p.Ny,Grid.p.Nx),'r','LevelList',5e-2)
            xlabel('x-dir, 1')
            ylabel('z-dir, 1')
            c.Label.String = 'Melt fraction, 1';
            
            %Tracer location concentration plot          
            ax3 = subplot(3,3,3);
            cla;
            axis equal
            hold on
            c = colorbar('NorthOutside');
            contour(X,Y,reshape(phi,Grid.p.Ny,Grid.p.Nx),'r','LevelList',5e-2),hold on
            trc1_plot = trc1; %trc1_plot(trc1<=1e-8) = nan;
            contourf(X,Y,reshape(trc1_plot,Grid.p.Ny,Grid.p.Nx),40,'linestyle','none'),view(2)
            xlabel('x-dir, 1')
            ylabel('z-dir, 1')
            c.Label.String = 'Tracer Conc.';
            colormap(ax3,flipud(gray));
            
            
            %average in radial  direction temperature in fraction plot along z
            subplot(3,3,4)
            cla;
            plot(mean(reshape(T,Grid.p.Ny,Grid.p.Nx),2),Grid.p.yc)
            xlabel('Avg. temp');
            ylabel('z-dir, 1');
            
            %total melt remaining as a function of time          
            subplot(3,3,5)
            cla;
            hold on
            plot(phiFracRem);
            ylabel('Total melt remaining, \%');
            xlabel('Time, yrs')

            %dimensionless thermal conductivity at each time
            subplot(3,3,6)
            cla;
            hold on
            axis equal
            contourf(X,Y,reshape(kappaPrimePlot,Grid.p.Ny,Grid.p.Nx),40,'linestyle','none'),view(2),hold on
            c = colorbar('NorthOutside');
            xlabel('x-dir, 1')
            ylabel('z-dir, 1')
            c.Label.Interpreter = 'latex';
            c.TickLabelInterpreter = 'latex';
            c.Label.String = 'Non-dim thermal conductivity';
            
            %Melt drained through the two planes to calculate the total
            %melt drained from the ocean
            [Xc,Yf] = meshgrid(Grid.p.xc,Grid.p.yf);
            subplot(3,3,7)
            cla;
            plot(tVec,phiDrain1Vec);
            hold on
            plot(tVec,phiDrain2Vec);
            legend('intefrace','1 above interface','Location','NorthWest');
            ylabel('Melt drained, m$^3$');
            xlabel('Time, yrs')
            
            %Percentage of melt drained through the two planes to 
            %calculate the total melt drained from the ocean
            subplot(3,3,8)
            cla;
            plot(tVec,phiDrain1Vec/phiOrig*100);
            hold on
            plot(tVec,phiDrain2Vec/phiOrig*100);
            legend('interface','1 above interface','Location','NorthWest');
            ylabel('Percentage of melt drained, \%');
            xlabel('Time, yrs')
            
            %Total dimensionless buoyancy
            subplot(3,3,9)
            cla;
            axis equal
            contourf(Xc,Yf,reshape(-fsVec,Grid.p.Ny+1,Grid.p.Nx),40,'linestyle','none'),view(2),hold on
            c = colorbar('NorthOutside');
            c.Label.Interpreter = 'latex';
            c.TickLabelInterpreter = 'latex';
            c.Label.String = 'Total Bouyancy, 1';
            xlabel('x-dir, 1')
            ylabel('z-dir, 1')
            axis equal
            saveas(h,sprintf('../figures/fig%d.png',i));
            
            %}
            
            % convert the image to a frame
            frameno = frameno + 1;
            FF(frameno) = getframe(gcf) ;
         end
    end
    
%%%%
%% Making a video out of frames
 % create the video writer with fps of the original video
 Data_result= sprintf('../figures/case%s_t%syrs.avi',num2str(fn),num2str(tTot));
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
    
end