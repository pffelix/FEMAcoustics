% draw field
colormap(jet);

cmap=colormap;
[n m]=size(cmap);

Nl=60; % number of field lines equal to the number of colors in palette
delta_l=int8(n/Nl-rem(n,Nl)/Nl);
for i=1:Nl
    cmap1(i,:)=cmap((i-1)*delta_l+1,:);
end

delta_field=(max_field-min_field)/(Nl-1);
for i=1:Nl
   field_lines(i)=min_field+(i-1)*delta_field; 
end

% field_lines=logspace(floor(log10(min_field)),ceil(log10(max_field)),Nl);

% Loop over elements

for i=1:Ne
    % Field values of nodes
    V1=V(el_no(i,1));
    V2=V(el_no(i,2));
    V3=V(el_no(i,3));
      
    if(V1<=V2 & V2<=V3)
        in1=1; % node of element with lowest field value
        in2=2;
        in3=3; % node of element with highest field values
    end
    
    if(V1<=V2 & V3<=V2 & V1<=V3)
        in1=1;
        in2=3;
        in3=2;
    end

    if(V1<=V2 & V3<=V2 & V3<=V1)
        in1=3;
        in2=1;
        in3=2;
    end

    if(V2<=V1 & V3<=V2 & V3<=V1)
        in1=3;
        in2=2;
        in3=1;
    end
    
    if(V2<=V1 & V2<=V3 & V1<=V3)
        in1=2;
        in2=1;
        in3=3;
    end
    
    if(V2<=V1 & V2<=V3 & V3<=V1)
        in1=2;
        in2=3;
        in3=1;
    end
    
    V1=V(el_no(i,in1)); % lowest field value
    V2=V(el_no(i,in2));
    V3=V(el_no(i,in3)); % highest field values
    
    x1=x_no(el_no(i,in1)); % x position node with lowest field value
    y1=y_no(el_no(i,in1)); % y position node with lowest field value
    x2=x_no(el_no(i,in2));
    y2=y_no(el_no(i,in2));
    x3=x_no(el_no(i,in3));
    y3=y_no(el_no(i,in3));
    
    for j=1:Nl
        if(field_lines(j)>=V1 & field_lines(j)>=V2 & field_lines(j)<=V3) % only the highest node field value is bigger than current field line j
            ratio=(field_lines(j)-V1)/(V3-V1);
            xt1=x1+ratio*(x3-x1);
            yt1=y1+ratio*(y3-y1);
            ratio=(field_lines(j)-V2)/(V3-V2);
            xt2=x2+ratio*(x3-x2);
            yt2=y2+ratio*(y3-y2);
            % line_plot(xt1,yt1,xt2,yt2,'k',1);
            line_plot(xt1,yt1,xt2,yt2,cmap1(j,:),1);
        end
        if(field_lines(j)>=V1 & field_lines(j)<=V2 & field_lines(j)<=V3)  % the highest and the second highest node field value is bigger than current field line j
            ratio=(field_lines(j)-V1)/(V3-V1);
            xt1=x1+ratio*(x3-x1);
            yt1=y1+ratio*(y3-y1);
            ratio=(field_lines(j)-V1)/(V2-V1);
            xt2=x1+ratio*(x2-x1);
            yt2=y1+ratio*(y2-y1);
            % line_plot(xt1,yt1,xt2,yt2,'k',1);
            line_plot(xt1,yt1,xt2,yt2,cmap1(j,:),1);
        end
    end
end

% Colorbar levels
Ncb_level=6;
delta_tick=n/(Ncb_level-1);
delta_cb=(max_field-min_field)/(Ncb_level-1);
for i=1:Ncb_level
    cb_tick(i)=1+int8((i-1)*delta_tick);
    cb_values(i)=min_field+(i-1)*delta_cb;
end
cb_tick(1)=2;
cb_tick(Ncb_level)=n;
% colorbar('YLim',[min_field max_field],'YTickLabel',cb_values);
% colorbar('ALim',[min_field max_field],'CLim',[min_field max_field],'YLim',[min_field max_field]);
% colorbar('East','YTick',cb_tick,'YTickLabel',cb_values);
% for i_bar=1:length(cb_values)

cb_values_str=strcat(num2str(cb_values,'%0.0f\n'),' dB');
c=colorbar('South','XTickLabel',cb_values_str,'FontSize',14);
c.Position=[0.1,0.05,0.83,0.05];
title_1=['RMS sound pressure level (SPL) for ', num2str(freq(nf),'%0.0f\n'),' Hz'];
% title_2 = ['average SPL in room = ',num2str(L_P_average_dB,'%0.0f'),' dB'];
if (model_Z_by_alpha)
    title_2 = ['absorption coefficient wall = ',num2str(alpha(nf,1),'%0.2f')];
else
    title_2 = ['specific impedance wall Z/Z0 = ',num2str(1/beta(nf,1)/Z0,'%0.0f')];
end
title_plot = title({title_1;title_2},'fontsize',20);
