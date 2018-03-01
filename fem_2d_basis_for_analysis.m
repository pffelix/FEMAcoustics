% FEM 2-D Scalar Analysis
% 2D Acoustics: Model sound pressure level in a given rectangular room geometry with
% scattering object excited by a radial sound source in a casing
% Input Parameters: Frequency dependent source conditions and boundary conditions
% FEM approach: Galerkin method is used, linear triangular shape functions
% Plot: geometry, mesh, boundary nodes, and sound pressure at all nodes expressed
% as logarithmic unit sound pressure level (SPL) with 2e-5 Pa as reference (hearing threshold)

% interesting plots:
% # freq=34, wall absorption coefficient low: first x-axial room resonance (also called x-axial room mode)
% # freq=68, wall absorption coefficient low: second x-axial resonance, strong local pressure node in front of sound source where sound field can hardly be excited
% # freq=110, wall absorption coefficient low: tangential (2,2) mode, very strong anti-nodes and increase of sound pressure in corners
% # By increasing the wall absorption coefficient at the above examples or reducing the corresponding acoustical impedance the resonances get less distinct (Q factors of the resonances decrease)
% # By further increasing the source frequency more and more resonances are locally excited and the sound field gets more homogeneous (diffuse)


% Load ASCII data: mesh, boundary nodes, domains (COMSOL is used as a mesh
% generator)
load_data;
set(gcf,'color','w');

%% 
%==================================================================
% Parameters definition (adaptable)
% **************************************************************

%% Acoustical constants of medium
rho0 = 1.2; % Air density in kg/m^3
c0 = 340; % Speed of sound in m/s
air_damp = 0.00; % Damping by the cavity air
p_0 = 2e-5; % reference sound pressure (hearing threshold 1kHz) in Pa
Z0 = rho0 * c0; % Specific impedance of fluid in Pa·s/m

%% Source parameters
freq = [400]; % Radial sound source (piston) frequency in Hz (can be a frequency vector, should be maximum 500 Hz otherwise less than 6 elements per wavelength)
u_n = [1e-5]; % Particle displacement in m by source excitation (vector corresponding to every frequency of excitation, adapt to change strength of sound)
source_xy = [0.6,2]; % x and y position of source in m (original source location [0.6,2])

%% Boundaries parameters
% Model boundaries either by practical measured absorption coefficients (simpler) or by physical more correct acoustical impedances under assumption of local reaction
model_Z_by_alpha = true ; % true: model acoustical impedances by absorption coefficients with Mommertz's (1996) method assuming the phase difference between sound pressure and normal particle velocity at boundaries is zero
% Literature: Mommertz, E. (1996): Untersuchung akustischer Wandeigenschaften und Modellierung der Schallrückwürfe in der binauralen Raumsimulation. Dissertation. RWTH Aachen, Fakultät für Elektrotechnik und Informationstechnik, S. 122. 

% A parameter matrix is needed with dimension nf x 3 (nf = number of source frequencies) 
% specify for every frequency 3 boundary condition material parameters (wall, scattering object, piston casing)
% either by:
alpha = [0.9,0.9,0.1]; % absorption coefficient at wall, scattering object, piston casing (expresses: sound energy absorbed per reflection in percent at boundary, value range: ]0,1[)

% or (if model_Z_by_alpha = false):
Z = [Z0*100,Z0*100,Z0*100]; % complex acoustical impedance at wall, scattering object, piston casing (expresses: sound pressure divided by normal particle velocity at boundary, value range: ]0,inf[ + 1i*]0,inf[ )

%% Solver parameters
simplify_damping_matrix = false; % use only diagonal elements of damping matrix
solver = 2; % 1 - sparse matrix solver; 2 - GMRES iterative solver

%% Plotting parameters
plot_propagation_video = false; % false: show steady state sound pressure field to analyse room resonances (faster), true: show sound propagation as video and make avi file (plotting for every frame is slow with standard plot function)
propagation_video_time = 2; % set movie time in milliseconds (computationally recommended 5-20)
propagation_video_framerate = 2; % set the number of frames per millisecond to calculate (computationally recommended 2)
plot_geometry_mesh_boundary = false; % plot additional figures about geometry, mesh and boundary
% *************************************************************************

%% 
%==================================================================
% Non-FEM calculations
% **************************************************************
%% Acoustical properties calculation

% Piston sound source
w = 2*pi*freq; % calculate angular frequency of loudspeaker excitation frequency
Nfreq = length(freq); % number of frequencies
lamda_min = c0/max(freq); % smallest wavelength in the room
k0 = w/c0; % wave number

% Boundaries
beta = 1./Z; % calculate admittances from impedances
Nmat = 3; % number of boundary materials

% Model acoustical admittances with Mommertz (1996) method
if(model_Z_by_alpha)
    for nf = 1:Nfreq
        fprintf(['Modelled specific acoustical impedance (Mommertz, 1996) at ', num2str(freq(nf),'%.0f'),'Hz for \n'])
        for mat = 1:Nmat
            switch mat
                case 1
                    fprintf('wall with absorption coefficient = %.2f:\n', alpha(nf,mat))
                case 2
                    fprintf('scattering object with absorption coefficient = %.2f:\n', alpha(nf,mat))
                case 3
                    fprintf('piston casing with absorption coefficient = %.2f:\n', alpha(nf,mat))
            end
            beta(nf,mat) = 1/(conversion_alpha_to_Z(alpha(nf,mat))*Z0);
            fprintf('Z/Z0 = %.0f \n',(1/beta(nf,mat))/Z0);
        end
    end
end

%% Get a theoretical calculation of modes that are closest to source frequency for comparison with the FEM result
analytical_modes; % print 3 closest axial and tangential modes to source frequency

%% Approximated calcuation of the number of elements per wavelength
S_e_average = L*W/Ne; % Average area of an element (rough approximation)
approx_average_length_element = sqrt(S_e_average*2); % Average diameter of an element (rough approximation)
elements_per_lamda_min = lamda_min/ approx_average_length_element;
if(elements_per_lamda_min<6)
    fprintf('Only %.2f elements per frequency (roughly approximated) \n',elements_per_lamda_min);
end
% *************************************************************************

%% 
%==================================================================
% FEM calculations
% **************************************************************
%% Step 1: Set up the mesh

% Mesh paramater calculation
[~, Ne_Nn] = size(elements); % number of nodes at each element
[~, N_dim] = size(nodes); % dimension (here x,y = 2)

% Size of figures
Px1 = 100;
Px2 = 600;
Py1 = 100;
Py2 = 620; 
K = sparse( Nn , Nn ); % Gobal stiffness matrix
M = sparse( Nn , Nn ); % Global mass matrix

% Get index of all non boundary nodes
nodes_no_bc_index = zeros(Nn,1);
nodes_no_bc_index(:,1) = linspace(1,Nn,Nn);
nodes_no_bc_index(unique(bc_elements)) = [];

% Get node index of source position
piston_node = nodes_no_bc_index(dsearchn(nodes(nodes_no_bc_index,:),delaunayn(nodes(nodes_no_bc_index,:)),source_xy)); % calculate closest node to piston point
% piston_node = 50; % dsearchn(nodes,delaunayn(nodes),piston_xy); 

%% Step 2: Assign boundary material to the corresponding mesh nodes

% assign selected material coefficients to boundary domains wall, object, piston casing

% create material vector with impedance entry for every node and every
% frequency
nodes_beta=zeros(Nn,Nfreq);

% get  index of every boundary node
nodes_bc_index = unique(bc_elements);
% get position of every boundary node
nodes_bc_xy = nodes(nodes_bc_index,:);

nodes_bc_material = {};
% get index of nodes representing wall
nodes_bc_material = [nodes_bc_material; nodes_bc_index(unique([find(nodes_bc_xy(:,1)==0); find(nodes_bc_xy(:,1)==L); find(nodes_bc_xy(:,2)==0); find(nodes_bc_xy(:,2)==W)]))];
% get index of nodes representing object
nodes_bc_material = [nodes_bc_material; nodes_bc_index(find(nodes_bc_xy(:,1)>=L/2 & nodes_bc_xy(:,1)<L  & nodes_bc_xy(:,2)>0 & nodes_bc_xy(:,2)<W))];
% get index of nodes representing piston casing
nodes_bc_material = [nodes_bc_material; nodes_bc_index(find(nodes_bc_xy(:,1)>0 & nodes_bc_xy(:,1)<L/2  & nodes_bc_xy(:,2)>0 & nodes_bc_xy(:,2)<W))];

% for the nodes at the boundaries (wall, object, piston casing) and for every frequency
for nf = 1:Nfreq
    for nb = 1:Nmat
        % assign the corresponding admittance condition
        nodes_beta(nodes_bc_material{nb},nf) = beta(nf,nb);
    end
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
    c = zeros(Ne_Nn,1);
    for e_n = 1:Ne_Nn
       e_n_xy(e_n,:) = nodes(el_no(e,e_n),:);
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
    area = abs(f(1)*c(2) - f(2)*c(1)) / 2; % y_23*x_13 - y_31*x_32
    
    % Derivation of the formulas for elementary matrices is given in literature:
    % Ochmann M., Mechel F. (2002) Analytical and numerical methods in acoustics. In: Mechel F (ed) Formulas of Acoustics, Springer Verlag, Berlin / Heidelberg, p. 1077.
    % or in a more detailed explained online at:
    % http://what-when-how.com/the-finite-element-method/fem-for-two-dimensional-solids-finite-element-method-part-1/
    
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
    

    %% Step 3b: Assemble to global stiffness and mass matrix 
    
    % positions where to assemble
    pos_e = sparse(Ne_Nn,Nn);
    for e_n = 1:Ne_Nn
        node_number = el_no(e,e_n);
        pos_e(e_n,node_number) = 1;
    end
    
    % assembling
    K = K + pos_e'*K_e*pos_e; % K Stiffness Matrix
    M = M + pos_e'*M_e*pos_e; % M Mass Matrix
end
toc

%% Step 4: Solve FEM matrices to get sound pressure at nodes

% Solve separated for every given frequency
for nf=1:Nfreq 
    % Get current angular frequency
    w_n = w(nf);
    fprintf('Frequency %.0f Hz: \n',freq(nf));
    fprintf('Integrate boundaries and force vector and solve weak form: \n');
    tic

    if(simplify_damping_matrix == false)
        %% Global damping matrix
        A = sparse( Nn , Nn );
        for b = 1:Nb
            % elementary damping matrix
            nodes_at_b = bc_elements(b,:);
            nodes_at_b_xy = nodes(nodes_at_b,:);
            h = sqrt((nodes_at_b_xy(1,1)-nodes_at_b_xy(2,1))^2+(nodes_at_b_xy(1,2)-nodes_at_b_xy(2,2))^2);
            A_b = [2,1;1,2]*h*beta(nf)*rho0; % assume both nodes at boundary element have same impedance


            % Assemble to global damping matrix 

            % positions where to assemble
            pos_b = sparse(Ne_Nn-1,Nn);
            for b_n = 1:Ne_Nn-1
                node_number = bc_elements(b,b_n);
                pos_b(b_n,node_number) = 1;
            end

            % assembling

            A = A + pos_b'*A_b*pos_b; % A Damping Matrix
        end
    else
        % integrate boundary conditions into a damping matrix A
        A = zeros(Nn,Nn);
        Boundary_vector = unique(bc_elements);
        diagIdx = 1:Nn+1:Nn*Nn; %get index of all diagonal elements in matrix with same shape as A
        A(diagIdx(Boundary_vector)) = nodes_beta(Boundary_vector,nf)*rho0; % assign admittances to diagonal elements
    end
    
    % Build force vector
    f = sparse(Nn,1); 
    f(piston_node) = w_n^2*rho0*u_n(nf);

    % Add global stiffness, mass and damping matrix to one matrix 
    Matrix = K - w_n^2*M/(1+1i*air_damp)+1i*w_n*A;

    % solve the system of equations with selected solver
    switch solver
            % ***************************
    case 1  % direct sparse matrix solver
            % ***************************
        spparms('autoamd',0);
        permut = colamd(Matrix); % Matrix reordering (bandwidth reduction)
        Vp = Matrix (permut,permut) \ f(permut); % Direct solution
        Vc(permut) = Vp; % Mapping to the original numbering 
            % **********************
    case 2  % GMRES iterative solver
            % **********************
        M1 = sparse(Nn,Nn);
        diag_A = diag(Matrix);
        for p = 1:Nn
           M1(p,p) = diag_A(p);
        end
        [Vp,flag,relres,iter] = gmres(Matrix,f,[],1e-7,Nn,M1);
        Vc = Vp';
    end    

    toc
    %% Step 5: Plot sound pressure field figure (in case of propagation video for each frame)
    fprintf('Plot sound pressure field: \n');
    tic
    % calculate how much time frames should be plotted 
    % select field variable to plot (video: sound pressure at time t, no-video: root man squared (RMS) sound pressure) 
    if plot_propagation_video == true
        frames_nr = propagation_video_time*propagation_video_framerate; % number of frames to plot
        t_video = linspace(0,propagation_video_time/1000,frames_nr); % Time vector for frames in seconds
        Vc_t = (Vc'*exp(-1i*t_video*w(nf)))'; % Complex sound pressure pointer at time t
        P_real_t = real(Vc_t); % Sound pressure at time t
        P_out = P_real_t; % Sound pressure at time t
    else
        P_magnitude = abs(Vc); % sound pressure magnitude
        P_out = P_magnitude/sqrt(2); % Root mean squared (RMS) sound pressure
        frames_nr = 1; % no video, only one figure to plot
    end
    
    % plot for each frame
    for plot_n = 1:frames_nr
        % calculate sound pressure level (log scale) from pressure input
        % (absolut scale)
        L_P_absolut_dB = 10*log10((P_out(plot_n,:).^2/p_0.^2))'; % sound pressure level (SPL) in dB
        % L_P_normalized_dB = 10*log10((P_out(plot_n,:).^2)/max(P_out(plot_n,:).^2))'; % normalized sound pressure level (SPL) in dB
        % L_P_average_dB = 10*log10(mean(P_out(plot_n,:).^2)/p_0.^2)'; % average sound pressure level (SPL) in room in dB

        % Plot the field with given standard plot function
        V = L_P_absolut_dB;
        fig1 = figure('Position',[Px1 Py1 Px2 Py2],'Color',[1 1 1]);
        set_figure_1;
        find_min_max;
        field_plot;
        geometry_plot;
        
        % save obtained time frame figure to propagation video matrix
        if plot_propagation_video == true
            fprintf('Frame number %.0f: \n',plot_n);
            newtitle = get(title_plot,'String');
            newtitle{1} = ['Sound pressure level (SPL) at ',num2str(t_video(plot_n)*1000,	'%.0f'),' ms',' for ', num2str(freq(nf),'%0.0f\n'),' Hz', ];
            set(title_plot,'String',newtitle);
            video_matrix(plot_n) = getframe(gcf);
            if (plot_n ~= frames_nr)
                close % close figure of current frame
            end
        end
    end
%% Step 6 (optional): Play and save video of sound propagation

    % after all frames calculated: play propagation video 3 times and save it to avi file
    if plot_propagation_video == true
        video_frame_rate=1;
        % play video
        movie(gcf,video_matrix,video_frame_rate,3)
        % save video
        v = VideoWriter(strcat('soundfield_',int2str(round(freq(nf))),'Hz','.avi'));
        v.FrameRate = video_frame_rate;
        v.Quality = 90;
        open(v);
        writeVideo(v,video_matrix);
        close(v);
    end
    toc
end
% *************************************************************************

%% 
%==================================================================
% Plot geometry, mesh, boundary nodes
% **************************************************************
if plot_geometry_mesh_boundary == true
    fprintf('Plot geometry, mesh and boundary: \n');
    tic
    if 1 % Geometry Plot
    figure('Position',[Px1 Py1 Px2 Py2],'Color',[1 1 1]);
    set_figure_1;
    geometry_plot;
    end
    
    if 1 % Mesh Plot
    figure('Position',[Px1 Py1 Px2 Py2],'Color',[1 1 1]);
    set_figure_1;
    mesh_plot;
    geometry_plot;
    end
    
    if 1 % Boundary Nodes Plot
    figure('Position',[Px1 Py1 Px2 Py2],'Color',[1 1 1]);
    set_figure_1;
    geometry_plot;
    plot_boundary_nodes;
    end
    toc
end
% *************************************************************************