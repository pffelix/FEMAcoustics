% Load and prepare data
clc;
clear;

load nodes.txt;
load elements.txt;
load domains.txt;
load bcs.txt;

% x= 3.93
% 470,512,513

[Nn m]=size(nodes);
[Ne m]=size(elements);
[Nb m]=size(bcs);

for i=1:Ne
   for j=1:3
       el_no(i,j)=elements(i,j)+1;  % which 3 nodes at element e
   end
   el_mat(i)=domains(i); % what material does element e has

end
for i=1:Nn
   x_no(i)=nodes(i,1); %  x position of node n
   y_no(i)=nodes(i,2); %  y position of node n
   st_no(i)=0; % is node n an boundary element: no
end

for i=1:Nb
   bc_elements(i,1)=bcs(i,1)+1; % node 1 and 2 of boundary line element b
   bc_elements(i,2)=bcs(i,2)+1;
   in=bcs(i,1)+1;
   st_no(in)=1;
   in=bcs(i,2)+1;
   st_no(in)=1; % is node n an boundary element: yes
end

