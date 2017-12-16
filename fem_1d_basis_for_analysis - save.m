% FEM 1-D Scalar Analysis
% 1D Acoustics
% Plot: geometry, mesh, boundary nodes, and field

% Generate mesh data: mesh, boundary nodes, domains (COMSOL is used as a mesh
% generator)
L=1; % Length of the tube


%% Input data
rho0=1.2; % Fluid density
c0=342.2; % Speed of sound
neta_a=0.00; % Damping in the cavity
U0=1; % Piston displacement
%% Impedance properties
Z0=rho0 * c0;
Z= input('\n Input the normalized complex impedance at one end : ');
if (Z~=0), beta= 1/Z;else beta=1e6; end % convert to admittance
%% Frequency domain
freq=[10:2:1000]; % Chosen arbitrarily
omega=2*pi*freq;
lamda_min= c0/max(freq); % Smallest wavelength in the fluid
%% Step 1: Mesh
nel=input ('\n Number of quadratic element / wavelength:');
ne= ceil(nel* L/lamda_min);
fprintf('\n Number of quadratic elements : %d \n', ne')

nnt = 2*ne+1; % Total number of nodes
h=L/ne; % Length of the elements
x=[0:h/2:L]; % Coordinates table
%% Step 2: Compute Elementary matrices
He=[7,-8,1;-8,16,-8;1,-8,7]/(3*h)/rho0;
Qe=[4,2,-1;2,16,2;-1,2,4]*h/30/(rho0*c0^2);
%% Step 3: Assembling
I=eye(3,3);
H=zeros(nnt,nnt); Q = zeros(nnt,nnt);
for ie=1: ne
LM=zeros(3,nnt); LM(:,2*ie-1:2*ie+1)=I;
H = H + LM'*He*LM;
Q = Q + LM'*Qe*LM;
end
%% Step 4: Impedance condition at x=0 : dPdn(0)+j*k*beta*P(0) = 0
A = zeros(nnt,nnt);
A(1,1) = beta/(rho0*c0);
%% Step 5: Solving the system with the Force vector : Piston at a x=L
ndof = nnt;
Pquad_direct=zeros(1,ndof);
Pquad_modal=zeros(1,ndof);
Pquad_exact=zeros(1,ndof);
for n=1:length(omega)
    w=omega(n);
    % FE solution
    R=zeros(ndof,1); % Initialize force vector
    R(ndof) = w^2*U0; % displacement of the piston fixed to 1
    P = (H - w^2*Q/(1+1i*neta_a)+1i*w*A)\R; % Solve the system
    Pquad_direct(n) = (rho0*c0^2)* real(P'*Q*P)/(2*L); % Compute the space avergaed quadratic pressure
    % Analytical solution
    k0=w/c0/sqrt(1+1i*neta_a);
    r=(Z-1)/(Z+1);
    a=(1j*w*U0)*Z0/(exp(1j*k0*L)-r*exp(-1j*k0*L));
    x=[0:L/100:L]; p=a*(exp(-1i*k0*x)+r*exp(1i*k0*x));
    Pquad_exact(n) = real(norm(p)^2/length(x)/2);
end

% Step 6: Comparison with the exact solution
plot(freq,10*log10(Pquad_direct(1:length(freq))),'k','LineWidth',2)
hold on
plot(freq,10*log10(Pquad_exact(1:length(freq))),':','LineWidth',2,'Color',[0.5 0.5 0.5])
xlabel('Frequency (Hz)');
ylabel(' Quadratic pressure (dB ref.1)');
legend('Analytical', 'FEM ');
text(100 , 95,['Analytical vs. FEM using ', num2str(ne), ' quadratic elements'],'Color','k','FontSize',12)
text(100 , 90,['Specified normalized Impedance Z = ', num2str(Z)],'FontSize',12)