% FEM 2-D Scalar Analysis
% 2D Acoustics: Model sound pressure level in a room excited by a piston
% Galerkin method is used, linear triangular shape functions
% Plot: geometry, mesh, boundary nodes, and field

% Load ASCII data: mesh, boundary nodes, domains (COMSOL is used as a mesh
% generator)
load_data;
set(gcf,'color','w');

% to do:
%     save datein löschen, dateien paramter initialisieren, zeilen korrigiern, delete git

% interesting plots:
% # freq=34 (x-axial mode (n_x=1), very strong resonance)
% # freq=68 (x-axial mode (n_x=2), even with piston sound pressure does not increase significantly at left wave node )
% # freq=110 (tangential (2,2) mode, very strong nodes and increase of sound pressure in corners)
% # By increasing the frequency the wave equation allows more resonances and so the sound field gets more homogenous (diffuse)
% # By increasing the absorption degrees or reducing the impedance of the boundary the standing waves get less distinct (the Q factors of the resonances decrease)


%% 
%==================================================================
% Parameters definition
% **************************************************************

%% Acoustical constants
rho0 = 1.2; % Air density in kg/m^3
c0 = 340; % Speed of sound in m/s
air_damp = 0.00; % Damping by the cavity air
p_0 = 2e-5; % reference sound pressure (hearing threshold 1kHz)
Z0 = rho0 * c0; % Specific impedance of fluid

%% Source parameters
freq = [68,120]; % Sound source (piston) frequency in Hz (can be a frequency vector, should be maximum 500 Hz otherwise less then 6 elements per wavelength at given model)
u_n = repelem(1e-7,length(freq)); % "Strenght of sound source": expressed by particle displacement in m (can be a vector corresponding to every frequency of excitation)
piston_xy = [0.6,2]; % x and y position of source in m (original source location [0.6,2])

%% Boundaries paramters
% Model boundaries either by practial measured absorption coefficients or by physical more correct impedances under assumption of local reaction
model_Z_by_alpha = true ; % true: model acoustical impedances by absorption coefficients with Mommertz's (1996) method assuming phase difference p/v at boundaries is zero
% Mommertz, E. (1996): Untersuchung akustischer Wandeigenschaften und Modellierung der Schallrückwürfe in der binauralen Raumsimulation. Dissertation. RWTH Aachen, Fakultät für Elektrotechnik und Informationstechnik, S. 122. 

% either:
alpha = [0.9,0.1,0.1]; % absorption coefficent at wall, scattering object, piston casing (expresses: sound energy absorbed per reflection in percent, value range: ]0,1[)

% or (for experts):
Z = [Z0*100,Z0*100,Z0*100]; % acoustical impedance at wall, scattering object, piston casing (expresses: sound pressure divided by particle velocity at boundary, value range: ]0,inf[ + i*]0,inf[ )

%% Solver parameters
solver = 1; % 1 - sparse matrix solver; 2 - GMRES iterative solver
video = true; % false: show steady state sound pressure field to analyse resonances, true: show sound propagation as video (plotting each frame is computationally slow with the standard plot function)
video_time = 100; % set movie time in milliseconds
video_framerate = 2; % set the number of frames per millisecond to calculate
% *************************************************************************
%% Non-FEM calculations
%% Acoustical properties calculation

% Piston source
w=2*pi*freq; % Calculate angular frequency of loudspeaker excitation frequency
Nfreq = length(freq); % number of frequencies
lamda_min= c0/max(freq); % Smallest wavelength in the room
k0 = w /c0; % wave number

% Boundaries
beta = 1./Z; % Calculate admittances from impedances
Nmat = 3; % number of boundary materials
% Model acoustical admittances with Mommertz (1996) method
% Mommertz, E. (1996): Untersuchung akustischer Wandeigenschaften und Modellierung der Schallrückwürfe in der binauralen Raumsimulation. Dissertation. RWTH Aachen, Fakultät für Elektrotechnik und Informationstechnik, S. 122. 

if(model_Z_by_alpha)
    fprintf('Modelled specific acoustical impedance (Mommertz, 1996) for \n');
    for mat = 1:Nmat
        switch mat
            case 1
                fprintf('wall with absorption coefficient = %.2f:\n', alpha(mat))
            case 2
                fprintf('scattering object with absorption coefficient = %.2f:\n', alpha(mat))
            case 3
                fprintf('piston casing with absorption coefficient = %.2f:\n', alpha(mat))
        end
        beta(mat)=1/(Mommertz(alpha(mat))*Z0);
        fprintf('Z/Z0 = %.0f \n',(1/beta(mat))/Z0);
    end
end

%% Analytical room modes calculation for comparision of FEM results and theoretical results for rectangular rooms:
% Vorländer, M. (2008). Auralization: Fundamentals of Acoustics, Modelling, Simulation,
% Algorithms and Acoustic Virtual Reality. Berlin Heidelberg: Springer Verlag, p. 55.
[x_min_bc,y_min_bc,y_max_bc,x_max_bc,S]=find_limits(Nn,x_no,y_no); % get room wall positions
L=x_max_bc-x_min_bc; % Length of the room in meter
W=y_max_bc-y_min_bc; % Width of the room in meter

n_modes=100;
f_ax_mode=zeros(10,1);
f_ay_mode=zeros(10,1);
f_tan_mode=zeros(10,10);

for nf=linspace(1,n_modes)
    f_ax_mode(nf)=c0/2*sqrt((nf/L)^2);
    f_ay_mode(nf)=c0/2*sqrt((nf/W)^2);
    for j = linspace(1,n_modes)
        f_tan_mode(nf,j)=c0/2*sqrt((nf/L)^2+(j/W)^2);
    end
end

for nf=1:length(w)
    fprintf('Frequency %.0f Hz: \n',w(nf));
    % print closest modes
    nr_closest_modes=3;
    min_ax=abs(f_ax_mode-freq(nf));
    min_ay=abs(f_ay_mode-freq(nf));
    min_tan=abs(f_tan_mode-freq(nf));
    for ax_i=1:nr_closest_modes
        % x-axial
        [~,ax_mode_nearest_idx] = min(min_ax);
        ax_mode_nearest_f=f_ax_mode(ax_mode_nearest_idx,1);
        fprintf('%.0f',ax_i);
        fprintf('. closest x-axial mode (n_x=%.0f) at %.0f Hz\n',[ax_mode_nearest_idx,ax_mode_nearest_f]);
        min_ax(ax_mode_nearest_idx,1)=inf;
    end
    for ax_i=1:3
        % y-axial
        [~,ay_mode_nearest_idx] = min(min_ay);
        ay_mode_nearest_f=f_ay_mode(ay_mode_nearest_idx,1);
        fprintf('%.0f',ax_i);
        fprintf('. closest y-axial y mode (n_y=%.0f) at %.0f Hz\n',[ay_mode_nearest_idx,ay_mode_nearest_f]);
        min_ay(ay_mode_nearest_idx,1)=inf;
    end
    for ax_i=1:3
        % xy-tangential
        [~,tan_mode_nearest_idx]=min(min_tan(:));
        [row,col]=ind2sub(size(min_tan),tan_mode_nearest_idx);
        tan_mode_nearest_idx=[row,col];
        tan_mode_nearest_f=f_tan_mode(tan_mode_nearest_idx(1),tan_mode_nearest_idx(2));
        fprintf('%.0f',ax_i);
        fprintf('. closest xy-tangential mode (n_x=%.0f,n_y=%.0f) at %.0fHz\n',[tan_mode_nearest_idx,tan_mode_nearest_f]);
        min_tan(tan_mode_nearest_idx(1),tan_mode_nearest_idx(2))=inf;
    end
end

%% Calcuation of the approximated number of elements per wavelength
S_e_average=L*W/Ne; % Average area of an element (rough approximation)
approx_average_length_element = sqrt(S_e_average*2);
elements_per_lamda_min= lamda_min/ approx_average_length_element;
if(elements_per_lamda_min<6)
    fprintf('frequency to high for given number of elements');
end

%% FEM calculations
%% Step 1: Set up the mesh

% Mesh paramater calculation
[~, Ne_Nn] = size(elements); % number of nodes at each element
[~, N_dim] = size(nodes); % dimension (here x,y = 2)
Px1=100;Px2=600;Py1=100;Py2=620; % Size of figures
K = sparse( Nn , Nn ); % Gobal stiffness matrix
M = sparse( Nn , Nn ); % Global mass matrix

% Get index of all non boundary nodes
nodes_no_bc_index=zeros(Nn,1);
nodes_no_bc_index(:,1)=linspace(1,Nn,Nn);
nodes_no_bc_index(unique(bc_elements)) = [];

% Get node index of piston
piston_node = nodes_no_bc_index(dsearchn(nodes(nodes_no_bc_index,:),delaunayn(nodes(nodes_no_bc_index,:)),piston_xy)); % calculate closest node to piston point
% piston_node = 50; % dsearchn(nodes,delaunayn(nodes),piston_xy); 

%% Step 2: Assign boundary material to the corresponding mesh nodes

% assign specific materials to boundary domains wall, object, piston casing
nodes_beta=zeros(Nn,1);

nodes_bc_index = unique(bc_elements);
nodes_bc_xy = nodes(nodes_bc_index,:);

nodes_bc_material = {};
% nodes representing wall
nodes_bc_material=[nodes_bc_material; nodes_bc_index(unique([find(nodes_bc_xy(:,1)==0); find(nodes_bc_xy(:,1)==L); find(nodes_bc_xy(:,2)==0); find(nodes_bc_xy(:,2)==W)]))];
% nodes representing object
nodes_bc_material=[nodes_bc_material; nodes_bc_index(find(nodes_bc_xy(:,1)>=L/2 & nodes_bc_xy(:,1)<L  & nodes_bc_xy(:,2)>0 & nodes_bc_xy(:,2)<W))];
% nodes representing piston casing
nodes_bc_material=[nodes_bc_material; nodes_bc_index(find(nodes_bc_xy(:,1)>0 & nodes_bc_xy(:,1)<L/2  & nodes_bc_xy(:,2)>0 & nodes_bc_xy(:,2)<W))];

for nf = 1:Nfreq
    nodes_beta(nodes_bc_material{nf})= beta(nf);
end

%% Step 3: Build and assemble stiffness and mass matrix
fprintf('Build and assemble stiffness and mass matrix: \n');
tic
for e = 1:Ne
    %% Step 3a: Build elementary stiffness and mass matrices
    % Shape functions: linear triangular with 3 nodes
    % Galerkin method
    
    % Get global x y position of all nodes of element
    e_n_xy = zeros(Ne_Nn,N_dim);
    f = zeros(Ne_Nn,1);
    c= zeros(Ne_Nn,1);
    for e_n = 1:Ne_Nn
       e_n_xy(e_n,:)= nodes(el_no(e,e_n),:);
    end
    
    %  y local coordinates nodes of element
    f(1) = e_n_xy(2,2) - e_n_xy(3,2); % y_23
    f(2) = e_n_xy(3,2) - e_n_xy(1,2); % y_31
    f(3) = e_n_xy(1,2) - e_n_xy(2,2); % y_12

    % x local coordinates nodes of element
    c(1) = e_n_xy(3,1) - e_n_xy(2,1); % x_32
    c(2) = e_n_xy(1,1) - e_n_xy(3,1); % x_13
    c(3) = e_n_xy(2,1) - e_n_xy(1,1); % x_21

    % area of triangular element
    area  = abs(f(1)*c(2) - f(2)*c(1)) / 2; % y_23*x_13 - y_31*x_32

    % elementary strain matrix
    B = [f(1),f(2), f(3);
          c(1),c(2),c(3)]/area/2;
    
    % elementary stiffness matrix
    K_e = sparse(transpose(B)*B*area);

    % elementary mass matrix
    M_e = [ 2 , 1 , 1 ;
         1 , 2 , 1 ;
         1 , 1 , 2 ];
    M_e = sparse(M_e * area / 12 / c0^2);
    
    % Find the equation number list for the i-th element
    % eqnum = el_no(e,:); % the 3 node numbers at element e
    % for i = 1 : Ne_Nn
      % for j = 1 : Ne_Nn
        % K(eqnum(i),eqnum(j)) = K(eqnum(i),eqnum(j)) + K_e(i,j);
        % M(eqnum(i),eqnum(j)) = M(eqnum(i),eqnum(j)) + M_e(i,j);
      % end
    % end

    %% Step 3b: Assemble to global stiffness and mass matrix 
    
    % positions where to assemble
    pos_e=sparse(Ne_Nn,Nn);
    for e_n= 1:Ne_Nn
        node_number = el_no(e,e_n);
        pos_e(e_n,node_number) = 1;
    end
    
    % assembling
    K = K + pos_e'*K_e*pos_e; % K Stiffness Matrix
    M = M + pos_e'*M_e*pos_e; % M Mass Matrix
end
toc

%% Step 4: Solve FEM matrices to get sound pressure at nodes

for nf=1:Nfreq % for every frequency
    % Get current angular frequency
    w_n = w(nf);
    fprintf('Frequency %.0f Hz: \n',w(nf));
    fprintf('Integrate boundaries and force vector and solve weak form: \n');
    tic

%         t_video=linspace(0,1,100);
%         
%         for t_video_i=t_video
%             if video
%                 t_factor=real(exp(1i*t_video_i*w(nf)));
%             else
%                 t_factor=1



    % integrate impedance conditions for current frequency into damping matrix
    A = sparse(Nn,Nn);
    Boundary_vector_non_unique = [bc_elements(:,2);bc_elements(:,1)];
    Boundary_vector = unique([bc_elements(:,2);bc_elements(:,1)]);
    [Bn, Bm] = size(A);
    diagIdx = 1:Bn+1:Bn*Bm;

    A(diagIdx(Boundary_vector)) = nodes_beta(Boundary_vector)*rho0;

    % Build force vector for angular frequency
    f = sparse(Nn,1); 
    f(piston_node) = w_n^2*rho0*u_n(nf);

    % Add global stiffness, mass and boundary matrix to one matrix 
    Matrix = K - w_n^2*M/(1+1i*air_damp)+1i*w_n*A;

    % solve the system of equations numerical
    switch solver
            % ***************************
    case 1  % direct sparse matrix solver
            % ***************************
        spparms('autoamd',0);
        permut=colamd(Matrix); % Matrix reordering (bandwidth reduction)
        Vp = Matrix (permut,permut) \ f(permut); % Direct solution
        Vc(permut)=Vp; % Mapping to the original numbering 
            % **********************
    case 2  % GMRES iterative solver
            % **********************
        M1=sparse(Nn,Nn);
        diag_A=diag(Matrix);
        for p=1:Nn
           M1(p,p)=diag_A(p); 
        end
        [Vc,flag,relres,iter]=gmres(Matrix,f,[],1e-7,Nn,M1);
    end    

    toc
    %% Step 5: Plot sound pressure field
    fprintf('Plot sound pressure field: \n');
    tic

    if video == true
        frames_nr = video_time*video_framerate; % number of frames to plot
        t_video = linspace(0,video_time/1000,frames_nr); % Time vector for frames in seconds
        Vc_t = (Vc'*exp(-1i*t_video*w(nf)))'; % Complex sound pressure pointer at time t
        P_real_t = real(Vc_t); % Sound pressure at time t
        P_out=P_real_t; % Sound pressure at time t
    else
        P_magnitude = abs(Vc); % sound pressure magnitude
        P_out = P_magnitude/sqrt(2); % Root mean squared (RMS) sound pressure
        frames_nr = 1; % number of frames to plot
    end
    
    for plot_n=1:frames_nr
        fprintf('Frame number %.0f: \n',plot_n);
        % calculate sound pressure level of steady state field
        L_P_absolut_dB = 10*log10((P_out(plot_n,:).^2/p_0.^2))'; % absolut sound pressure level in dB
        L_P_normalized_dB = 10*log10((P_out(plot_n,:).^2)/max(P_out(plot_n,:).^2))'; % normalized sound pressure level in dB
        L_P_average_dB = 10*log10(mean(P_out(plot_n,:).^2)/p_0.^2)'; % average sound pressure level in room in dB

        % assign which field value to plot
        V=L_P_absolut_dB;
        fig1=figure('Position',[Px1 Py1 Px2 Py2],'Color',[1 1 1]);
        set_figure_1;
        find_min_max;
        field_plot;
        geometry_plot;
        
        if video == true
        % save figure for video
            newtitle = get(title_plot,'String');
            newtitle{end+1}=['Time: ',num2str(t_video(plot_n)*1000,	'%.0f'),' ms'];
            set(title_plot,'String',newtitle)
            video_matrix(plot_n) = getframe(gcf);
            if (plot_n ~= frames_nr)
                close % close figure of current frame
            end
        end
    end
    if video == true
        video_frame_rate=1;
        % play video
        movie(gcf,video_matrix,video_frame_rate,3)
        % save video
        v = VideoWriter(strcat('soundwave_',int2str(round(w(nf)/2/pi)),'Hz','.avi'));
        v.FrameRate = video_frame_rate;
        v.Quality = 99;
        open(v);
        writeVideo(v,video_matrix);
        close(v);
    end
    toc
end



% %% strg r, t
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
% 
