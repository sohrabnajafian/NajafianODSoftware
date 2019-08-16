function functionOtherImageAnalysis(handles)


InputImage = handles.InputImage;
rowPoint  = handles.SelectedPoint(1);
colPoint  = handles.SelectedPoint(2);

if size(InputImage,3) == 3 
    input_gray = rgb2gray(InputImage);
else
    input_gray = InputImage;
end 
input_binary = im2bw(input_gray,graythresh(input_gray));


if input_binary(rowPoint, colPoint) == 1
    [L_ipsi,N_region_ipsi] = bwlabel(input_binary);
elseif input_binary(rowPoint, colPoint) == 0
    [L_ipsi,N_region_ipsi] = bwlabel(~input_binary);
end
%figure,imagesc(L_ipsi)

i_ipsi = L_ipsi(rowPoint, colPoint);
selected_region = L_ipsi  == i_ipsi ;
debug_thickness = 0;
pixel2um = 1;

step_hist  = 30;
initial_edge = -(step_hist/2);
last_edge = 360 + initial_edge;
edge_range = (initial_edge:step_hist:last_edge);
bin_plot= (0:step_hist:360) *(pi/180);


if ( sum(sum(selected_region))>10 ) % to remove noises
    
    [n_point,angle_line,thick,point_coordiantes] = horton_stripes_info2(selected_region,debug_thickness,pixel2um);
    %[n_point,angle_line,thick] = horton_stripes_info2(selected_region,debug_thickness,pixel2um);
    
    npoint_ipsi(i_ipsi) = n_point;
    thickness_ipsi(i_ipsi) = mean(thick);


    cla(handles.axesImage,'reset')
    axes(handles.axesImage)
    imshow(InputImage)
    title(sprintf(' Medial Axis Length = %.0f pixels',npoint_ipsi(i_ipsi)),'fontsize',18)
    [selected_region_row,selected_region_col] = find(edge(selected_region));
    hold on, plot(selected_region_col,selected_region_row,'.r')
    
    nbins_thickness_ipsi = 6; 
    cla(handles.axesResult)
    set(handles.axesResult,'visible','on')
    axes(handles.axesResult)
    thick(thick<0) = 0;
    
    if input_binary(rowPoint, colPoint) == 1
        [N_thick_ipsi,plot2] = hist_mid_line(thick,nbins_thickness_ipsi,'w');
    elseif input_binary(rowPoint, colPoint) == 0
        [N_thick_ipsi,plot2] = hist_mid_line(thick,nbins_thickness_ipsi,'k');
    end

    
    axis([.9*min(thick) 1.1*max(thick) -.05 1])
    xlabel('Width ','fontsize',14)
    ylabel('Frequency','fontsize',14)
    title(sprintf('                                               Width  %.2f pixels ',thickness_ipsi(i_ipsi)),'fontsize',16)
    ax_width = gca;
    get( ax_width );
    set( ax_width, 'Color', [0.7,0.7,0.7] )
    ax_width.TickDir = 'out';
    ax_width.TickLength = [0.02 0.02];
    box(ax_width,'off')

    axes(handles.axesOrientation)
    set(handles.axesOrientation,'visible','on')
    angle_line_0_180 = angle_line .* (angle_line>0) + (angle_line+180) .* (angle_line<0);
    angle_line_full_range_I = cat(2,angle_line_0_180,angle_line_0_180+180);
    [frequency_ipsi,bin_edge_ipsi] = histcounts(angle_line_full_range_I * (pi/180) ,edge_range * (pi/180),'Normalization','Probability'); %edge is the bar start and end points
    frequency_ipsi(1) = frequency_ipsi(floor(end/2)+1);% 0 and 180 should be equal
    frequency_ipsi(end+1) = frequency_ipsi(1); %    for smooth plotting
    
    if input_binary(rowPoint, colPoint) == 1
         polar_histogram_predefined_edge_makevideo(bin_plot,frequency_ipsi,max(frequency_ipsi),'-',[255 255 255]/255);
    elseif input_binary(rowPoint, colPoint) == 0
         polar_histogram_predefined_edge_makevideo(bin_plot,frequency_ipsi,max(frequency_ipsi),'-',[0 0 0]/255);
    end
    
    polar_histogram_predefined_edge_makevideo(bin_plot,frequency_ipsi,max(frequency_ipsi),'-',[0 0 0]/255);
    xlabel('      Angle','fontsize',14)
    

end





