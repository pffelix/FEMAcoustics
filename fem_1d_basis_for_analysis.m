%
% Problem data
%
Lx=1; % Length of the domain; fixed arbitrary to 1
c0=342; % Speed of sound (m/s)
rho0=1.2; % Fluid density (kg/m^3))
neta_a= input(‘ loss factor (damping) in the fluid: ’); % Play with this parameter
freq=[10:2:1000]; % Frequency domain
omega=2*pi*freq;
%
% Step 1: Mesh
%
ne= input(‘ Number of quadratic element: ’); % number of elements
nnt = 2*ne+1; % Total number of nodes
h=Lx/ne; % Length of the elements
x=[0:h/2:Lx]; % Coordinates table
%
% Step 2: Compute Elementary matrices
%
Ke=c0^2*[7,-8,1;-8,16,-8;1,-8,7]/(3*h);
Me=[4,2,-1;2,16,2;-1,2,4]*h/30;
%
% Step 3: Assembling
%
I=eye(3,3);
K=zeros(nnt,nnt); M = zeros(nnt,nnt);
for ie=1: ne
    L=zeros(3,nnt); L(:,2*ie-1:2*ie+1)=I; % Location matrix for element ie
    K = K + L’*Ke*L;
    M = M + L’*Me*L;
end
%
% Step 4: Boundary conditions: p(1) = 0
%
K = K(2:nnt,2:nnt);
M = M(2:nnt,2:nnt);
ndof = nnt-1; % Final number of unknown (equations)
%
% Step 5: Solving the system with the Force vector : Piston at a x=L
%
R=zeros(1,ndof); % Initialize the force vector
Pquad=zeros(1,ndof);
Pquad_exact=zeros(1,ndof);
for n=1:length(omega)
    w=omega(n);

    R(ndof) = w^2*rho0*c0^2*(1+1i*neta_a); % Piston (fixed displacement) at x=L
    MK=K*(1+1i*neta_a)-w^2*M;
    P=inv(MK)*R’; % Solve the linear system
    Pquad(n) = real((P’*M*P)/2); % Compute the space avergaed quadratic pressure
    % Analytical solution (to be derived as an exercise)
    k0=w/c0/sqrt(1+1i*neta_a);
    amp=sin(k0*x);
    Pquad_exact(n) = real(rho0*c0*sqrt(1+1i*neta_a)*w/cos(k0*Lx))^2 .* norm(amp)^2/length(x)/2;
end
%
% Step 6: Comparison with the exact solution
%
plot(freq,10*log10(Pquad_exact),freq,10*log10(Pquad),‘r-’);
xlabel(‘Frequency (Hz)’);
ylabel(‘ Quadratic pressure (dB)’);
legend(‘Analytical’, ‘FEM ’);
sprintf(‘Analytical vs. FEM: %d quadratic elements; /eta=%f’, ne, neta_a);
text=sprintf(‘Analytical vs. FEM using %d quadratic elements; damping =%4.2f %%’, ne, neta_a);
Title(text)