% FEM 1-D Scalar Analysis
% 1D Acoustics: Model sound pressure level in a waveguide with piston at one end and wall at other end
% Plot: geometry, mesh, boundary nodes, and field


%%
%==================================================================
% Parameters definition
% **************************************************************

% Generate mesh data: mesh, boundary nodes, domains 
L=1; % Length of the tube
shape_type = 2; % 1 - linear shape functions; 2 - quadratic shape functions 
boundaries_x = [0,L]; %  wall boundaries x-position (vector can have multiple possible, 0<=x<=L)
pistons_x = L; %  piston sources x-position (vector can have multiple positions, 0<=x<=L)

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
alpha = [0.9,0.9]; % absorption coefficent at every boundary (sound energy absorbed per reflection in percent, value range ]0,1[)
% For experts ( acoustical impedance Z=p/v respectively Admittance A=v/p in the praxis normally not measured/available):
% set model_Z_by_alpha==false -> 
Z = [2055-4678i,10000]; % acoustical impedance at every boundary (=sound pressure divided by particle velocity at boundary, value range ]0,inf[)

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
Nmat = length(boundaries_x); % number of boundary materials

% Model acoustical admittances with Mommertz (1996) method
% Mommertz, E. (1996): Untersuchung akustischer Wandeigenschaften und Modellierung der Schallrückwürfe in der binauralen Raumsimulation. Dissertation. RWTH Aachen, Fakultät für Elektrotechnik und Informationstechnik, S. 122. 
if(model_Z_by_alpha)
    fprintf('Modelled acoustical impedance (Mommertz, 1996) for \n');
    for mat = 1:Nmat
        fprintf('boundary at x = %.3f with absorption coefficient = %.2f:\n', [boundaries_x(mat) alpha(mat)])
        beta(mat)=1/(Mommertz(alpha(mat))*Z0);
        fprintf('%.0f Pa s/m3 \n',1/beta(mat));
    end
end

%% Set up Mesh
Ne_per_lamda_min = 20; % minimum number of elements per wavelength
Ne = ceil(Ne_per_lamda_min * L / lamda_min); % Total number of elements 
Nn = shape_type*Ne + 1; % Total number of nodes
Ne_Nn = shape_type + 1; % Number of nodes per element
h=L/Ne; % Element length
x = 0:h/shape_type:L; % Coordinates of nodes
% get nodes at pistons x-position
boundary_nodes = 0;
piston_nodes = 0;
for i=1:length(pistons_x)
    [~, piston_nodes(i)] = min(abs(x - pistons_x(i)));
end
% get nodes at boundaries x-position
for i=1:length(boundaries_x)
    [~, boundary_nodes(i)] = min(abs(x - boundaries_x(i)));
end
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
for i = 1:length(boundary_nodes)
    A(boundary_nodes(i),boundary_nodes(i)) = beta(i)*rho0; % *rho0
end
%% Add force vector and solve
L_P_absolut_dB=zeros(Nn,Nfreq);
L_P_normalized_dB=zeros(Nn,Nfreq);
L_P_average_dB=zeros(Nn,Nfreq);
P_real=zeros(Nn,Nfreq);

for nf=1:Nfreq
    w_n=w(nf);
    f=zeros(Nn,1); % Force vector
    f(piston_nodes) = w_n^2*rho0*u_n(nf); % at piston position
    Matrix = K - w_n^2*M/(1+1i*air_damp)+1i*w_n*A;
    Vc = Matrix\f; % Solve the system to get pressure in tube
    P_average_rms = (rho0*c0^2)* real(Vc'*M*Vc)/(2*L); % space avergaed quadratic pressure /(rho0*c0^2) because of M
    
    P_magnitude = abs(Vc); % sound pressure magnitude
    P_phase = angle(Vc); % sound pressure phase
    P_real(:,nf)=real(P_magnitude.*exp(1i*P_phase*w_n));
    % L_P_absolut_dB(:,nf) = 20*log10((P_magnitude/sqrt(2))/p_0)'; % absolut sound pressure level in dB
    % L_P_normalized_dB(:,nf) = 20*log10((P_magnitude/sqrt(2))/max(P_magnitude/sqrt(2))); % normalized sound pressure level in dB
    L_P_average_dB(:,nf) = 20*log10(mean(P_magnitude/sqrt(2))/p_0); % average sound pressure level in room in dB
    
    
end

% Plot sound pressure level in waveguide for selected frequencies
freq_to_plot=[200]; % which frequencies should be plotted
for nf_plot= 1:length(freq_to_plot)
    freq_to_plot_i=ismember(freq,freq_to_plot(nf_plot));
    plot(x, P_real(:,freq_to_plot_i)) 
    hold on
    xlabel('x-position');
    ylabel('sound pressure');

end

title_1=['Sound pressure in waveguide at ', num2str(freq(nf),'%0.0f\n'),' Hz'];
title_2 = ['waveguide average SPL = ',num2str(L_P_average_dB,'%0.0f'),' dB'];
title({title_1;title_2})
figure
plot(freq,10*log10(P_average_rms(1:length(freq))))
hold on
xlabel('Frequency in Hz');
ylabel(' RMS sound pressure');
