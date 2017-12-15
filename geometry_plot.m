% plot geometry

for i=1:Nb % Loop over boundary line elements
    in1=bc_elements(i,1);
    in2=bc_elements(i,2);
    x1=x_no(in1);x2=x_no(in2);
    y1=y_no(in1);y2=y_no(in2);
    line_plot(x1,y1,x2,y2,'k',2);
end
