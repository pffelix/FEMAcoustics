%% Analytical room modes calculation for comparision of FEM results and theoretical results for rectangular rooms:
% Vorländer, M. (2008). Auralization: Fundamentals of Acoustics, Modelling, Simulation,
% Algorithms and Acoustic Virtual Reality. Berlin Heidelberg: Springer Verlag, p. 55.
[x_min_bc,y_min_bc,y_max_bc,x_max_bc,S]=find_limits(Nn,x_no,y_no); % get room wall positions
L=x_max_bc-x_min_bc; % Length of the room in meter
W=y_max_bc-y_min_bc; % Width of the room in meter

% number of modes to predict
n_modes=100;
% initialize output vectors
f_ax_mode=zeros(10,1);
f_ay_mode=zeros(10,1);
f_tan_mode=zeros(10,10);

% calculate modes analytically
for nf=linspace(1,n_modes)
    f_ax_mode(nf)=c0/2*sqrt((nf/L)^2);
    f_ay_mode(nf)=c0/2*sqrt((nf/W)^2);
    for j = linspace(1,n_modes)
        f_tan_mode(nf,j)=c0/2*sqrt((nf/L)^2+(j/W)^2);
    end
end

% show the 3 modes as output that are closest to the source excitation frequency
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