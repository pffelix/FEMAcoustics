% FEM 1-D Scalar Analysis
% 1D Acoustics waveguide with piston at one end and wall at other end
% Plot: geometry, mesh, boundary nodes, and field


%%
%==================================================================
% Parameters definition
% **************************************************************

% Generate mesh data: mesh, boundary nodes, domains 
L=1; % Length of the tube
shape_type = 2; % 1 - linear shape functions; 2 - quadratic shape functions 
boundaries_x = 0; %  boundaries x-position (vector with multiple positions possible)
pistons_x = L; %  Piston sources x-position (vector with multiple positions possible)

% Acoustical constants
rho0 = 1.2; % Air density in kg/m^3
c0 = 340; % Speed of sound in m/s
air_damp = 0.00; % Damping by the cavity air
p_0 = 2e-5; % reference root man squared sound pressure (hearing threshold 1kHz)

% Piston source parameters
freq = [10:2:1000]; % Frequencies piston excitation in Hz (can be a frequency vector, should be maximum 500 Hz otherwise less then 6 elements per wavelength)
u_n = repelem(1/10000000,length(freq)); % sound particle displacement in m at every frequency of piston excitation (adapt to change radiated sound power)

% Boundary material paramters
model_Z_by_alpha = true ; % model acoustical impedances by practical more commonly measured random incidence absorption coefficients with Mommertz's (1996) method assuming phase difference p/v at boundaries is zero
% Mommertz, E. (1996): Untersuchung akustischer Wandeigenschaften und Modellierung der Schallrückwürfe in der binauralen Raumsimulation. Dissertation. RWTH Aachen, Fakultät für Elektrotechnik und Informationstechnik, S. 122. 
alpha = repelem(0.1,length(freq)); % absorption coefficent wall (sound energy absorbed per reflection in percent, value range ]0,1[)
% For experts ( acoustical impedance Z=p/v respectively Admittance A=v/p in the praxis normally not measured/available):
% set model_Z_by_alpha==false -> 
Z = repelem(2055-4678i,length(freq)); % acoustical impedance at wall (=sound pressure divided by particle velocity at boundary, value range ]0,inf[)

% *************************************************************************

%% Acoustical properties calculation

% Piston source
w = 2*pi*freq; % Calculate angular frequency of loudspeaker excitation frequency
Nfreq = length(freq); % number of frequencies
lamda = c0./freq; % Calculate wave length of loudspeaker excitation frequency
lamda_min = c0/max(freq); % Smallest wavelength in the room
k0 = w /c0; % wave number

% Boundaries
Z0 = rho0 * c0; % Calculate impedance of air
beta= 1./Z; % Calculate admittances from impedances
% Model acoustical admittances with Mommertz (1996) method
% Mommertz, E. (1996): Untersuchung akustischer Wandeigenschaften und Modellierung der Schallrückwürfe in der binauralen Raumsimulation. Dissertation. RWTH Aachen, Fakultät für Elektrotechnik und Informationstechnik, S. 122. 

if(model_Z_by_alpha)
    fprintf('fodelled acoustical impedance (Mommertz, 1996) for ');
    for nf = 1:Nfreq
        fprintf('frequency %.0f Hz: \n',freq(nf));
        beta(nf)=1/(Mommertz(alpha(nf))*Z0);
        fprintf('%.0f Pa s/m3 \n',1/beta(nf));
    end
end




%% Generate Mesh
Ne_per_lamda_min = 6; % minimum number of elements per wavelength
Ne = ceil(Ne_per_lamda_min * L / lamda_min); % Total number of elements 
Nn = shape_type*Ne + 1; % Total number of nodes
Ne_Nn = shape_type + 1; % Number of nodes per element
h=L/Ne; % Element length
x = 0:h/shape_type:L; % Coordinates 
[~, piston_node] = min(abs(x - piston_x));
[~, boundary_node] = min(abs(x - boundary_x));
%% Compute elementary matrices
switch shape_type
    case 1  
        K_e = [1,-1;-1,1]*1/h;
        M_e = [2,1;1,2]*h/6/c0^2;
    case 2
        K_e = [7,-8,1;-8,16,-8;1,-8,7]/(3*h);
        M_e = [4,2,-1;2,16,2;-1,2,4]*h/30/c0^2;
end
%% Assemble stiffness and mass matrix

I=sparse(eye(Ne_Nn,Ne_Nn));
K=sparse(Nn,Nn); M = sparse(Nn,Nn);
for e=1:Ne
    pos_e = sparse(Ne_Nn,Nn); 
    pos_e(:,(Ne_Nn-1)*e-(Ne_Nn-2):(Ne_Nn-1)*e+1)=I;
    K = K + pos_e'*K_e*pos_e;
    M = M + pos_e'*M_e*pos_e;
end
    
%% Add boundaries
A = zeros(Nn,Nn);
A([boundary_node],[boundary_node]) = beta*rho0; % *rho0
%% Add force vector and solve
P_average_rms=zeros(1,Nn);
for nf=1:Nfreq
    w_n=w(nf);
    f=zeros(Nn,1); % Force vector
    f(piston_node) = w_n^2*rho0*u_n(nf); % at piston position
    Matrix = K - w_n^2*M/(1+1i*air_damp)+1i*w_n*A;
    P = Matrix\f; % Solve the system to get pressure in tube
    P_average_rms(nf) = (rho0*c0^2)* real(P'*M*P)/(2*L); % space avergaed quadratic pressure /(rho0*c0^2) because of M
end

% Step 6: Comparison with the exact solution
plot(freq,10*log10(P_average_rms(1:length(freq))))
hold on