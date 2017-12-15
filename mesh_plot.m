% mesh plot

for i=1:Ne
    % if(el_mat(i)==1)
        line_plot(x_no(el_no(i,1)),y_no(el_no(i,1)),x_no(el_no(i,2)),y_no(el_no(i,2)),'b',1);
        line_plot(x_no(el_no(i,2)),y_no(el_no(i,2)),x_no(el_no(i,3)),y_no(el_no(i,3)),'b',1);
        line_plot(x_no(el_no(i,3)),y_no(el_no(i,3)),x_no(el_no(i,1)),y_no(el_no(i,1)),'b',1);
    % end
end

hold on;
