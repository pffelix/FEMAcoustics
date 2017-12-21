% FEM 1-D Scalar Analysis
% 1D Acoustics: Model sound pressure level in a waveguide with piston at one end and wall at other end
% Plot: geometry, mesh, boundary nodes, and field

set(gcf,'color','w');

%%
%==================================================================
% Parameters definition
% **************************************************************

% Generate mesh data: mesh, boundary nodes, domains 
L=1; % Length of the tube
shape_type = 2; % 1 - linear shape functions; 2 - quadratic shape functions 
boundary_x = L; %  wall boundary x-position (0<=x<=L)
piston_x = 0; %  piston source x-position (0<=x<=L)

% Acoustical constants
rho0 = 1.2; % Air density in kg/m^3
c0 = 340; % Speed of sound in m/s
air_damp = 0.00; % Damping by the cavity air
p_0 = 2e-5; % reference root man squared sound pressure (hearing threshold 1kHz)
Z0 = rho0 * c0; % Specific impedance of fluid

% Piston source parameters
freq = [2:2:1000]; % Frequencies piston excitation in Hz (can be a frequency vector, should be maximum 500 Hz otherwise less then 6 elements per wavelength)
u_n = repelem(1,length(freq)); % sound particle displacement in m at every frequency of piston excitation (adapt to change radiated sound power)

% Boundary material paramters
model_Z_by_alpha = false ; % model acoustical impedances by practical more commonly measured random incidence absorption coefficients with Mommertz's (1996) method assuming phase difference p/v at boundaries is zero
% Mommertz, E. (1996): Untersuchung akustischer Wandeigenschaften und Modellierung der Schallr�ckw�rfe in der binauralen Raumsimulation. Dissertation. RWTH Aachen, Fakult�t f�r Elektrotechnik und Informationstechnik, S. 122. 
alpha = 0.9; % absorption coefficent at boundary (sound energy absorbed per reflection in percent, value range ]0,1[)
% For experts ( acoustical impedance Z=p/v respectively Admittance A=v/p in the praxis normally not measured/available):
% set model_Z_by_alpha==false -> 
Z = Z0*100; % acoustical impedance at boundary (=sound pressure divided by particle velocity at boundary, value range: ]0,inf[ + ]0,inf[*i )

% *************************************************************************

%% Acoustical properties calculation

% Piston source
w = 2*pi*freq; % Calculate angular frequency of loudspeaker excitation frequency
Nfreq = length(freq); % number of frequencies
lamda = c0./freq; % Calculate wave length of loudspeaker excitation frequency
lamda_min = c0/max(freq); % Smallest wavelength in the room
k0 = w /c0; % wave number

% Boundaries
beta= 1./Z; % Calculate admittances from impedances
Nmat = length(boundary_x); % number of boundary materials

% Model acoustical admittances with Mommertz (1996) method
% Mommertz, E. (1996): Untersuchung akustischer Wandeigenschaften und Modellierung der Schallr�ckw�rfe in der binauralen Raumsimulation. Dissertation. RWTH Aachen, Fakult�t f�r Elektrotechnik und Informationstechnik, S. 122. 
if(model_Z_by_alpha)
    fprintf('Modelled specific acoustical impedance (Mommertz, 1996) for \n');
    for mat = 1:Nmat
        fprintf('boundary at x = %.3f with absorption coefficient = %.2f:\n', [boundary_x(mat) alpha(mat)])
        beta(mat)=1/(Mommertz(alpha(mat))*Z0);
        fprintf('Z/Z0 = %.0f \n',(1/beta(mat))/Z0);
    end
end

%% Set up Mesh
Ne_per_lamda_min = 6; % minimum number of elements per wavelength
Ne = ceil(Ne_per_lamda_min * L / lamda_min); % Total number of elements 
Nn = shape_type*Ne + 1; % Total number of nodes
Ne_Nn = shape_type + 1; % Number of nodes per element
h=L/Ne; % Element length
x = 0:h/shape_type:L; % Coordinates of nodes
% get nodes at pistons x-position
boundary_nodes = 0;
piston_nodes = 0;
for i=1:length(piston_x)
    [~, piston_nodes(i)] = min(abs(x - piston_x(i)));
end
% get nodes at boundaries x-position
for i=1:length(boundary_x)
    [~, boundary_nodes(i)] = min(abs(x - boundary_x(i)));
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
    A(boundary_nodes(i),boundary_nodes(i)) = beta(i)*rho0; %
end
%% Add force vector and solve
L_P_average_dB=zeros(1,Nfreq);
P_real=zeros(Nn,Nfreq);

for nf=1:Nfreq
    w_n=w(nf);
    f=zeros(Nn,1); % Force vector
    f(piston_nodes) = w_n^2*rho0*u_n(nf); % at piston position
    Matrix = K - w_n^2*M/(1+1i*air_damp)+1i*w_n*A;
    Vc = Matrix\f; % Solve the system to get pressure in tube
    P_magnitude = abs(Vc); % sound pressure magnitude
    P_phase = angle(Vc); % sound pressure phase
    L_P_average_dB(:,nf) = 20*log10(mean(P_magnitude/sqrt(2))/p_0); % average sound pressure level in room in dB
    
end

% Plot sound pressure level in waveguide for selected frequencies
figure(100)
plot(freq,10*log10(L_P_average_dB(:,1:length(freq))),'LineWidth', 2, 'DisplayName',['Z/Z0 = ',num2str(1/beta(1)/Z0,'%0.0f')])
hold on
title_1=['Sound pressure level SPL in waveguide'];
% if (model_Z_by_alpha)
    % title_2 = ['wall absorption coefficient = ',num2str(alpha(1),'%0.2f'), '; modelled specific wall impedance Z/Z0 = ',num2str(1/beta(1)/Z0,'%0.0f')];
% else
    % title_2 = ['specific wall impedance Z/Z0 = ',num2str(1/beta(1)/Z0,'%0.0f')];
% end
% title({title_1;title_2})
xlabel('Frequency in Hz','fontweight','bold','fontsize',16);
ylabel('Sound pressure level SPL','fontweight','bold','fontsize',16);
% legend(gca,'show')
