% FEM 2-D Scalar Analysis
% 2D Acoustics
% Plot: geometry, mesh, boundary nodes, and field

% Load ASCII data: mesh, boundary nodes, domains (COMSOL is used as a mesh
% generator)
load_data;

% to do:
% 2D: Präsi: hohe f, niedrige f vergleichen (extreme druckunterschiede), hohe und niedrige Impedanz Wand, andere
%     Quellen position, rechteckram ohne boundaries, raummoden zeigen, Z->0
%     (freifeld, sommerfeld ansprechen),Time componten einbauen (moden
%     schwingung zeigen), optimize performance ansprechen
%     Probleme: wieso läuft durch boundary, el_mat, domains, H-Q ersetzen
%     einheitlich, Solver integrieren (+ander TEchnik Highlights)
% 1D: pressure über x zeigen und pressure über f, Z variieren

%% 
%==================================================================
% Parameters definition
% **************************************************************
rho0 = 1.2; % Air density in kg/m^3
c0 = 340; % Speed of sound in m/s
air_damp = 0.00; % Damping by the cavity air
Z = 100; % % Wall impedance defined as pressure at wall divided by particle velocity Z=p/v
freq = 500; % Frequency of piston excitation in Hz (should be maximum 500 Hz otherwise less then 6 elements per wavelength)
solver = 1; % 1 - sparse matrix solver; 2 - GMRES iterative solver
piston_xy = [0.55,2]; %  Piston position x and y (original position source [0.55,2])

% *************************************************************************
%% Acoustical properties calucation

% Boundary
Z0 = rho0 * c0; % Calculate impedance of air
if (Z~=0), beta= 1/Z; else beta=1e6; end % Calculate wall admittance

% Excitation
w=2*pi*freq; % Calculate angular frequency of loudspeaker excitation frequency
lamda = c0/freq; % Calculate wave length of loudspeaker excitation frequency
lamda_min= c0/max(freq); % Smallest wavelength in the room
k0 = w /c0; % wave number

% Room modes
[x_min_bc,y_min_bc,y_max_bc,x_max_bc,S]=find_limits(Nn,x_no,y_no); % get room wall positions
L=x_max_bc-x_min_bc; % Length of the room in meter
W=y_max_bc-y_min_bc; % Width of the room in meter

n_modes=100;
f_ax_mode=zeros(10,1);
f_ay_mode=zeros(10,1);
f_tan_mode=zeros(10,10);

for i=linspace(1,n_modes)
    f_ax_mode(i)=c0/2*sqrt((i/L)^2);
    f_ay_mode(i)=c0/2*sqrt((i/W)^2);
    for j = linspace(1,n_modes)
        f_tan_mode(i,j)=c0/2*sqrt((i/L)^2+(j/W)^2);
    end
end

% Approximated number of elements per wavelength 
S_e_average=L*W/Ne; % Average area of an element (rough approximation)
approx_average_length_element = sqrt(S_e_average*2);
elements_per_lamda_min= lamda_min/ approx_average_length_element;
if(elements_per_lamda_min<6)
    fprintf('frequency to high for element number');
end

%% Step 1: Set up mesh

% Mesh paramater calculation
[~, Ne_Nn] = size(elements); % number of nodes at each element
[~, N_dim] = size(nodes); % dimension (here x,y = 2)
Px1=100;Px2=600;Py1=100;Py2=620; % Size of figures
K = sparse( Nn , Nn ); % Gobal stiffness matrix
M = sparse( Nn , Nn ); % Global mass matrix

% Get index of all non boundary nodes
nodes_no_bc_index=zeros(Nn,1);
nodes_no_bc_index(:,1)=linspace(1,Nn,Nn);
nodes_no_bc_index(unique(bc_elements),:) = [];

% Get node index of piston
piston_node = nodes_no_bc_index(dsearchn(nodes(nodes_no_bc_index,:),delaunayn(nodes(nodes_no_bc_index,:)),piston_xy)); % calculate closest node to piston point
% piston_node = 50; % dsearchn(nodes,delaunayn(nodes),piston_xy); 

%% Modify mesh

% make room without interior objects

% bc_elements_compelete=bc_elements;
% bc_inside_index =unique([find(nodes(:,1)==0); find(nodes(:,2)==0)]); % find alle node indexes which have x or y == 0
% bc_elements=bc_elements(unique([find(ismember(bc_elements(:,2),bc_inside_index)==1); find(ismember(bc_elements(:,1),bc_inside_index)==1)]),:); % only store nodes with this indexes in bc_elements

%% Build stiffness and mass matrix
fprintf('Assemble stiffness and mass matrix: \n');
tic
for e = 1:Ne
    %% Compute elementary stiffness and mass matrices
    % shape functions for linear triangular with 3 nodes 
    % Get global x y position of all nodes of element
    e_n_xy = zeros(Ne_Nn,N_dim);
    b = zeros(Ne_Nn,1);
    c= zeros(Ne_Nn,1);
    for e_n = 1:Ne_Nn
       e_n_xy(e_n,:)= nodes(el_no(e,e_n),:);
    end
    
    %  y local coordinates
    b(1) = e_n_xy(2,2) - e_n_xy(3,2); % y_23
    b(2) = e_n_xy(3,2) - e_n_xy(1,2); % y_31
    b(3) = e_n_xy(1,2) - e_n_xy(2,2); % y_12

    % x local coordinates
    c(1) = e_n_xy(3,1) - e_n_xy(2,1); % x_32
    c(2) = e_n_xy(1,1) - e_n_xy(3,1); % x_13
    c(3) = e_n_xy(2,1) - e_n_xy(1,1); % x_21

    % area of triangular element
    area2 = abs(b(1)*c(2) - b(2)*c(1)); % y_23*x_13 - y_31*x_32
    area  = area2 / 2;

    % elementary stiffness matrix
    B = [b(1),b(2), b(3);
          c(1),c(2),c(3)];
    B = B / area2;
    Ke = sparse(transpose(B)*B*area);

    % elementary mass matrix
    Me = [ 2 , 1 , 1 ;
         1 , 2 , 1 ;
         1 , 1 , 2 ];
    Me = sparse(Me * area / 12);
    
    % Me = [ 1 , 0 , 0 ;
           % 0 , 1 , 0 ;
           % 0 , 0 , 1 ];
    % Me = Me * area / 3 ;
    % Find the equation number list for the i-th element
    % eqnum = el_no(e,:); % the 3 node numbers at element e
    % for i = 1 : Ne_Nn
      % for j = 1 : Ne_Nn
        % K(eqnum(i),eqnum(j)) = K(eqnum(i),eqnum(j)) + Ke(i,j);
        % M(eqnum(i),eqnum(j)) = M(eqnum(i),eqnum(j)) + Me(i,j);
      % end
    % end

    
    %% Assemble global stiffness and mass matrix 
    
    % positions where to assemble
    LM=sparse(Ne_Nn,Nn);
    for e_n= 1:Ne_Nn
        node_number = el_no(e,e_n);
        LM(e_n,node_number) = 1;
    end
    
    % assembling
    K = K + LM'*Ke*LM; % K Stiffness Matrix
    M = M + LM'*Me*LM; % M Mass Matrix
end
toc

%% Solving the FEM system to get sound pressure in room
fprintf('Integrate boundaries and force vector and solve weak form for every frequency: \n');
tic

V = zeros(Nn,length(w));  % normalized squared sound pressure in dB
for nw=1:length(w)
    
    % Get current angular frequency
    w_n = w(nw);
    
    % integrate impedance conditions for frequency into damping matrix
    C = sparse(Nn,Nn);
    Boundary_vector_non_unique = [bc_elements(:,2);bc_elements(:,1)];
    Boundary_vector = unique([bc_elements(:,2);bc_elements(:,1)]);
    [Bn, Bm] = size(C);
    diagIdx = 1:Bn+1:Bn*Bm;
    C(diagIdx(Boundary_vector) ) = beta;

    % Build force vector for frequency
    b = sparse(Nn,1); 
    b(piston_node) = w_n^2;

    % Add global stiffness, mass and boundary matrix to one matrix 
    % representing the left side of the weak form of the Helmholtz equation
    A = K/rho0 - w_n^2*M/((rho0*c0^2)*(1+1i*air_damp))+1i*w_n*C/(rho0*c0);
    
    switch solver
            % ***************************
    case 1  % direct sparse matrix solver
            % ***************************
        spparms('autoamd',0);
        permut=colamd(A); % Matrix reordering (bandwidth reduction)
        Vp = A (permut,permut) \ b(permut); % Direct solution
        Vc(permut)=Vp; % Mapping to the original numbering 
            % **********************
    case 2  % GMRES iterative solver
            % **********************
        M1=sparse(Nn,Nn);
        diag_A=diag(A);
        for p=1:Nn
           M1(p,p)=diag_A(p); 
        end
        [Vc,flag,relres,iter]=gmres(A,b,[],1e-7,Nn,M1);
    end    
    P_magnitude = abs(Vc); % sound pressure magnitude
    L_P_dB = 20*log10(P_magnitude/max(P_magnitude)); % normalized squared sound pressure level in dB
    V(:,nw)=L_P_dB;
end
toc
%% Plot field
fprintf('Plot field: \n');
tic

if 1 % Field Plot
% V2=Field;
fig1=figure('Position',[Px1 Py1 Px2 Py2],'Color',[1 1 1]);
set_figure_1;
find_min_max;
field_plot;
geometry_plot;
end

toc
%% strg r, t
% if 1 % Geometry Plot
% figure('Position',[Px1 Py1 Px2 Py2],'Color',[1 1 1]);
% set_figure_1;
% geometry_plot;
% end
% 
% if 1 % Mesh Plot
% figure('Position',[Px1 Py1 Px2 Py2],'Color',[1 1 1]);
% set_figure_1;
% mesh_plot;
% geometry_plot;
% end
% 
% if 1 % Boundary Nodes Plot
% figure('Position',[Px1 Py1 Px2 Py2],'Color',[1 1 1]);
% set_figure_1;
% geometry_plot;
% plot_boundary_nodes;
% end

% Set the field values for testing the field visulaization function
% xc=0;yc=0;
% Field=zeros(Nn,1);
% for i=1:Nn
    % Field(i)=sqrt((x_no(i)-xc)^2+(y_no(i)-yc)^2);
% end
