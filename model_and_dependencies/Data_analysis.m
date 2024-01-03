%Calculating volume at each time

d = 20e3;   %Thickness of ice shell [in m]   


load("../Output/impact_Wakita_eta0_14_Ea_50_output_100Tlayer-17C.mat"); %loading file
interface = (find(Grid.p.yc>0)); interface = interface(1,1); %First cell of no ocean
Volume = [];
Time_arr = [];

%Individual file
load("../Output/impact_Wakita_eta0_14_Ea_50_output_100Tlayer-17C.mat"); %loading file
phi_arr = (reshape(phi,Grid.p.Ny,Grid.p.Nx)); %reshaping to find phi
phi_arr(1:interface-1,:) = 0;

Volume = [Volume;sum(sum(phi_arr(1:end,:),1).*Grid.p.V(Grid.p.dof_ymin)' * d^3)] %Calculating volume [in m^3]
Time_arr = [Time_arr;tVec(100)]