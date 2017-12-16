% FEM 2-D Scalar Analysis
% 2D Acoustics
% Plot: geometry, mesh, boundary nodes, and field

% Load ASCII data: mesh, boundary nodes, domains (COMSOL is used as a mesh
% generator)
load_data;

% to do:
% el_mat , domains
% A wieso complex
%%
%==================================================================
% Parameters definition
% **************************************************************
rho0 = 1.2; % Air density in kg/m^3
c0 = 340; % Speed of sound in m/s
air_damp = 0.00; % Damping by the cavity air
Z = 2055-4678i; % % Wall impedance defined as pressure at wall divided by particle velocity Z=p/v
freq = 100; % Frequency of loudspeaker excitation
solver = 1; % 1 - sparse matrix solver; 2 - GMRES iterative solver
piston_xy = [3,3]; %  Piston source position node
solver=1; % 1 - sparse matrix solver; 2 - GMRES iterative solver

% *************************************************************************

%% Step 1: Set up Mesh

% Mesh paramater calculation
[~, Ne_Nn] = size(elements);
[~, N_dim] = size(nodes);
Px1=100;Px2=600;Py1=100;Py2=620;    % Size of figures
[x_min_bc,y_min_bc,y_max_bc,x_max_bc,S]=find_limits(Nn,x_no,y_no); % get room wall positions
L=x_max_bc-x_min_bc; % Length of the room in meter
W=y_max_bc-y_min_bc; % Width of the room in meter

% make room without interior objects
% bc_elements_compelete=bc_elements;
% bc_inside_index =unique([find(nodes(:,1)==0); find(nodes(:,2)==0)]); % find alle node indexes which have x or y == 0
% bc_elements=bc_elements(unique([find(ismember(bc_elements(:,2),bc_inside_index)==1); find(ismember(bc_elements(:,1),bc_inside_index)==1)]),:); % only store nodes with this indexes in bc_elements

nodes_no_bc_index=zeros(Nn,1); % index vector with positions of all non boundary nodes
nodes_no_bc_index(:,1)=linspace(1,Nn,Nn);
nodes_no_bc_index(unique(bc_elements),:) = [];

piston_node = nodes_no_bc_index(dsearchn(nodes(nodes_no_bc_index,:),delaunayn(nodes(nodes_no_bc_index,:)),piston_xy)); % calculate closest node to piston point
% piston_node = 50; % dsearchn(nodes,delaunayn(nodes),piston_xy); 


% S_e_average=L*W/Ne; % Average area of an element
% approx_average_length_element = sqrt(S_e_average*2);
% elements_per_lamda_min= lamda_min/ approx_average_length_element;
% if(elements_per_lamda_min<6)
    % fprintf('frequency to high for element number');
%end

K = sparse( Nn , Nn );     % Gobal stiffness matrix
% C = sparse( Nn , Nn );     % Global damping matrix
M = sparse( Nn , Nn );     % Global mass matrix
force = complex( zeros(Nn,1) , zeros(Nn,1) ); % Force vector
% H = sparse(Nn,Nn); % Stiffness Matrix
% Q = sparse(Nn,Nn); % Mass Matrix
H = zeros(Nn,Nn); % Stiffness Matrix (rigidity)
Q = zeros(Nn,Nn); % Mass Matrix



%% Acoustical parameter calucation
% Boundary
Z0 = rho0 * c0; % Calculate impedance of air
if (Z~=0), beta= 1/Z; else beta=1e6; end % Calculate wall admittance

% Excitation
omega=2*pi*freq; % Calculate angular frequency of loudspeaker excitation frequency
lamda = c0/freq; % Calculate wave length of loudspeaker excitation frequency
lamda_min= c0/max(freq); % Smallest wavelength in the room
k0 = omega /c0; % wave number

%% Step 2: Compute Elementary matrices and %% Step 3: Assembling (not optimized sparse matrix)

for e = 1:Ne
   % Step 2: Compute shape functions for elementary matrices
   
   elementnodes= zeros(Ne_Nn,N_dim);
   b = zeros(Ne_Nn,1);
   c= zeros(Ne_Nn,1);
   for e_n = 1:Ne_Nn
       elementnodes(e_n,:)= nodes(el_no(e,e_n),:);
   end

    % calculate stiffness matrix of triangular element
    % y value
    b(1) = elementnodes(2,2) - elementnodes(3,2); % y_23
    b(2) = elementnodes(3,2) - elementnodes(1,2); % y_31
    b(3) = elementnodes(1,2) - elementnodes(2,2); % y_12

    % x value
    c(1) = elementnodes(3,1) - elementnodes(2,1); % x_32
    c(2) = elementnodes(1,1) - elementnodes(3,1); % x_13
    c(3) = elementnodes(2,1) - elementnodes(1,1); % x_21
    
    %Elemental values
    % Se=(c(3)*b(2)-c(2)*b(3))/2; %  det Jacobi Matrix /2
    % Hey=1/(imj*omega*mu0)*((y2-y3)*V1+(y3-y1)*V2+(y1-y2)*V3)/(2*Se);
    % Hex=-1/(imj*omega*mu0)*((x3-x2)*V1+(x1-x3)*V2+(x2-x1)*V3)/(2*Se);
    
    J_e_det = abs(b(1)*c(2) - b(2)*c(1)); % y_23*x_13 - y_31*x_32=a(1) % ander Schreibweise det Jacobi Matrix
    area  = J_e_det / 2;
    %=======

    J_e = [b(1),b(2), b(3); % correct Jacobi Matrix
          c(1),c(2),c(3)];

    J_e = J_e / J_e_det;

    Ke = transpose(J_e)*J_e*area; % area?
    

    %=======
    N_e = [ b(1)+c(1),b(2)+c(2),b(3)+c(3)]; % shape functions
    
    Me = [ 2 , 1 , 1 ;
         1 , 2 , 1 ;
         1 , 1 , 2 ];

    Me = k0^2 * Me * area / 12 ; % k0^2=1/c^2
    % Me = Ke;
    Me = [ 1 , 0 , 0 ;
         0 , 1 , 0 ;
         0 , 0 , 1 ];

    Me = k0^2 * Me * area / 3 ;
    
    % Shape function 

    
    % Step 3: Assembling (not optimized sparse matrix)

    % Find the equation number list for the i-th element
    lnods=el_no(e,:); % the 3 node numbers at element e
    eqnum = lnods;

    for i = 1 : Ne_Nn
      for j = 1 : Ne_Nn
        K(eqnum(i),eqnum(j)) = K(eqnum(i),eqnum(j)) + Ke(i,j);
        M(eqnum(i),eqnum(j)) = M(eqnum(i),eqnum(j)) + Me(i,j);
      end
    end
    LM=zeros(Ne_Nn,Nn);
    for e_n= 1:Ne_Nn
        % LM(e_n,2*ie-1:2*ie+1)=1;
        node_number = el_no(e,e_n);
        LM(e_n,node_number) = 1;
    end
    
    % H = H + sparse(LM'*Ke*LM); % K Stiffness Matrix
    % Q = Q + sparse(LM'*Me*LM); % M Mass Matrix
    H = H + LM'*Ke*LM; % K Stiffness Matrix
    Q = Q + LM'*Me*LM; % M Mass Matrix
 end

%% Step 4: Impedance condition at x=boundaries : dPdn(boundarie)+j*k*beta*P(boundarie) = 0 -> Here always the same impedance
% B = sparse(Nn,Nn);
C = zeros(Nn,Nn);
Boundary_vector_non_unique = [bc_elements(:,2);bc_elements(:,1)];
Boundary_vector = unique([bc_elements(:,2);bc_elements(:,1)]);
[Bn, Bm] = size(C);
diagIdx = 1:Bn+1:Bn*Bm;
C(diagIdx(Boundary_vector) ) = (1252 -2209i); % beta/(rho0*c0); %% NN fehlt, *rho, kein c0

%% Step 5: Solving the system
Pquad = zeros(Nn,length(omega));
% for n=1:length(omega)
n = 1;
w = omega(n);
% FE solution
S = zeros(Nn,1); % Initialize force vector
S(piston_node) = 100;% w^2; % U0 displacement of the piston fixed to 1 %% N fehlt, 1i fehlt

T_1= H - w^2*Q;
histfit(T_1(:),10)
min(T_1(:))
max(T_1(:))

A = (H - w^2*Q/(1+1i*air_damp)+1i*w*C);
b = S;
test=inv(w^2*Q);
max(abs(test(:,50)).^2)
min(abs(test(:,50)).^2)
size_A=size(A);
eye_matrix=eye(size_A(1));
P = A\S; % Solve the system direct


% P=Vc;
% Pquad_direct(n) = (rho0*c0^2)* real(P'*Q*P)/(2*L); % Compute the space avergaed quadratic pressure
Pquad = abs(P); % Compute the  average quadratic pressure at node n
% end

if 1 % Field Plot
V= log10(Pquad); % round((Pquad/min(Pquad)-1)/max((Pquad/min(Pquad)-1))*10000000);

% V2=Field;
fig1=figure('Position',[Px1 Py1 Px2 Py2],'Color',[1 1 1]);
set_figure_1;
find_min_max;
field_plot;
geometry_plot;
end



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
        M=diag(A);
        for p=1:Nn
           M1(p,p)=M(p); 
        end
        [Vc,flag,relres,iter]=gmres(A,b,[],1e-7,Nn,M1);
end









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
