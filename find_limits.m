function [x_min_bc,y_min_bc,y_max_bc,x_max_bc,R]=find_limits(Nnodes,x_nodes,y_nodes)
% Find the node (0,0)
x_min_bc=x_nodes(1);
x_max_bc=x_nodes(1);
y_min_bc=y_nodes(1);
y_max_bc=y_nodes(1);
R=sqrt(x_nodes(1)^2+y_nodes(1)^2);
for i=1:Nnodes
   if(x_min_bc>x_nodes(i)) 
       x_min_bc=x_nodes(i);
   end
   if(x_max_bc<x_nodes(i)) 
       x_max_bc=x_nodes(i);
   end
   if(y_min_bc>y_nodes(i)) 
       y_min_bc=y_nodes(i);
   end
   if(y_max_bc<y_nodes(i)) 
       y_max_bc=y_nodes(i);
   end
   if(R<sqrt(x_nodes(i)^2+y_nodes(i)^2))
       R=sqrt(x_nodes(i)^2+y_nodes(i)^2);
   end
end
return;
