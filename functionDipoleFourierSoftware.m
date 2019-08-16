function functionDipoleFourierSoftware(condition,handles)


%   Update : 021519
%   Every Patch should contain three cycles of stripes (31 by 31 pixels)
%   for comparison with database 
%   To measure the resize value for each map I used the excel sheet that
%   has the value of mean width for all the animals 

%load('database_modified040418.mat') % some convolution masks are omitted because of OD map generation problem

load('DatabaseLengthNorm.mat')
load('DatabaseNregion.mat')
load('DatabaseWidth.mat')
load('DatabaseOrientation.mat')
load('DatabaseResult.mat')
load('DatabaseTitle.mat')

%Control parameter
save_patches = 0;


%*******************************%
%condition = 1;
%*******************************%




if condition(1) == 1
    
    horton_human = imread('horton_retinotopy.jpg');%figure,imshow(horton_human)
    input_image = imbinarize(horton_human);
    resize_value = 0.364; %in comparison to 'horton_binary_gray2.jpg'
    %3 cycles per patch
    %input_image = imresize(input_image,resize_value);
    
    initial_row_range = 1;
    initial_col_range = 1;
        
    patch_width = 31;
    patch_step = 31;
    
    gray_region = horton_human >= 179 & horton_human <= 180;
    [Lhorton] = bwlabel(~gray_region);
    v1_region = Lhorton == 22; %main area
    %removing noises
    v1_region_resized = imresize(v1_region,resize_value);
    surrounding = ~ v1_region_resized;
    surrounding = imdilate(surrounding,ones(6));
    v1_region_resized = ~ surrounding;
    
    input_OD = imresize(input_image,resize_value) .* v1_region_resized;
    % figure,imshow(input_OD)
    
    horton_human_resized = imresize(horton_human,resize_value);
    %imshow(horton_human_resized)
    cla(handles.axesImage,'reset')
    axes(handles.axesImage)
    imshow(uint8(imresize(~input_image,resize_value) .* v1_region_resized *255+ double(horton_human_resized) .* (~v1_region_resized)))
    axis off
    
    %{
    [edge_row,edge_col] = find(edge(v1_region_resized));
    
    for u1 = 1:size(edge_row,1)
        hold on
        plot(edge_col(u1),edge_row(u1),'m','LineWidth',2)
    end
    
    hold on
    plot(edge_col,edge_row,'.k','MarkerSize',5)
    %}
    
    pixel2um =  209;

    
elseif condition(2) == 1 
    macaque_OD1 = imread('macaque retinotopy.jpg');
    macaque_OD1 = rgb2gray(macaque_OD1);
    
    resize_value = .463 / 1.688; %3 cycles per patch 
    
    input_image = imbinarize(macaque_OD1);
    input_OD = imresize(input_image,resize_value); %image almost the same size of horton image
    %figure,imshow(input_OD)
    
    initial_row_range = 1;%3;
    initial_col_range = 1;
    
    patch_width = 31;
    patch_step = 31;
    
    gray_region = macaque_OD1 >= 178 & macaque_OD1 <= 180; %gray region == 179
    [L] = bwlabel(~gray_region);
    %figure,imagesc(L)
    v1_region = L == 6;
    
    v1_region_resized = imresize(v1_region,resize_value);
    %figure,imshow(v1_region_resized)
    %pixel2um = 94.340 * resize_value;
    
    
    cla(handles.axesImage,'reset')
    axes(handles.axesImage)
    input_image_resized = imresize(macaque_OD1,resize_value);
    imshow(uint8(imresize(~input_image,resize_value) .* v1_region_resized *255+ double(input_image_resized) .* (~v1_region_resized)))

    %imshow(input_image_resized)
    %imshow(ones(size(v1_region_resized,1),size(v1_region_resized,2)))
    axis off
    
    %{
    [edge_row,edge_col] = find(edge(v1_region_resized));
    
    for u1 = 1:size(edge_row,1)
        hold on
        plot(edge_col(u1),edge_row(u1),'m','LineWidth',2)
    end
    
    hold on
    plot(edge_col,edge_row,'.k','MarkerSize',5)
    %}

    pixel2um = 97.5;

    
elseif condition(3) == 1 
    cat1 = imread('cat retinotopy2 with MC.jpg');
    cat1 = rgb2gray(cat1);
    input_image = imbinarize(cat1);
    
    resize_value = 0.43;%for dipole generation = .43;
    input_image_resized = imresize(input_image,resize_value); % each white/black stripe width = 28, horton = 12
    %figure,imshow(input_image_resized)
    % threshold should be
    % around around resize 12/28 =.43
    % the image could be resized
    % or size of patch
    % could be increased by
    % 1/.43
    
    %patch_width = 75;
    patch_width = 31;% 80 has the least minimum error in comparison of all the best match patches
    patch_step = 31; % 5 histogram analysis
    
    
    input_OD = ~input_image_resized;
    initial_row_range = 1;
    initial_col_range = 1;
    
    gray_region = cat1 >= 178 & cat1 <= 180; %gray region == 179
    [L] = bwlabel(~gray_region);
    %figure,imagesc(L)
    v1_region = L == 2;
    
    v1_region_resized = imresize(v1_region,resize_value);
    %figure,imshow(v1_region_resized)
    
    
    cla(handles.axesImage,'reset')
    axes(handles.axesImage)
    input_image_resized = imresize(input_image,resize_value);
    %imshow(input_image_resized)
    
    imshow( uint8(input_image_resized .* v1_region_resized*255 + double(input_image_resized)*180 .* (~v1_region_resized)) )

    %imshow(ones(size(v1_region_resized,1),size(v1_region_resized,2)))
    axis off
    [edge_row,edge_col] = find(edge(v1_region_resized));
    
    for u1 = 1:size(edge_row,1)
        hold on
        plot(edge_col(u1),edge_row(u1),'m','LineWidth',2)
    end
    
    hold on    
    %plot(edge_col,edge_row,'.k','MarkerSize',5)
    

    pixel2um = 113;
    %*************************************************%

    
end






counter1 = 0;

row_range = initial_row_range : patch_step : size(input_OD,1); % for histogram analysis the step is reduced to 30/6 = 5 to have smooth curve
col_range = initial_col_range : patch_step : size(input_OD,2);



orientation_OD_total_fourier = zeros(size(row_range,2),size(col_range,2));
mask_output{size(row_range,2),size(col_range,2)} = 0;

attract_x_OD_total= zeros(size(row_range,2),size(col_range,2));
attract_y_OD_total = zeros(size(row_range,2),size(col_range,2));
repulse_x_OD_total = zeros(size(row_range,2),size(col_range,2));
repulse_y_OD_total = zeros(size(row_range,2),size(col_range,2));

%plot and calculate for just indeces that are one 
ind_hist_plot = zeros(size(row_range,2),size(col_range,2));

%******************************************%
%*****Fourier Orientation******************%

row_input = patch_width;
col_input= patch_width;
mask_center = zeros(row_input,col_input); 
mask_center(round(row_input/2),round(col_input/2)) = 1; 
 
[x1,y1] = meshgrid(-floor(col_input/2) : floor(col_input/2), -floor(row_input/2): floor(row_input/2));  
 
matrix_angle = (180/pi) * atan(-y1./x1); 
matrix_angle = (matrix_angle < 0)*180 + (matrix_angle) ; 
matrix_angle(ceil(row_input/2),1:floor(row_input/2)) = 180; 
matrix_angle(ceil(row_input/2),ceil(col_input/2)) = 0; 
matrix_angle = matrix_angle .* (y1<=0);
 
%figure,imshow(a1)
 
mask_angle_row = zeros(36,ceil(row_input/2)); % 36 orientations starting from 5 to 180 degrees with the step of 5 degrees 
mask_angle_col = zeros(36,ceil(col_input/2));
mask_angle_row(:,1) =  ceil(row_input/2); 
mask_angle_col(:,1) =  ceil(col_input/2); 

for radius = 2:ceil(col_input/2)
    
    mask_dilation = imdilate(mask_center,ones(2*radius-1)) - imdilate(mask_center,ones(2*radius-3));
    mask_dilation = mask_dilation .* (y1<=0);
    
    mask_radius = matrix_angle .* mask_dilation;
    
    for orient_count = 1:36 % 5 to 180 degrees with the step of 5
        angle = ones(row_input,col_input)*1000;
        angle(mask_radius~=0) = orient_count * 5;
        
        [closest_angle,ind_min] = min(abs(angle(:) - mask_radius(:)));
        mask_angle_row(orient_count,radius) = mod(ind_min,row_input);
        mask_angle_col(orient_count,radius) = floor(ind_min/row_input) + 1;
    end
    
end


%{
figure,imshow(zeros(row_input,col_input))
hold on 
plot(mask_angle_col(1,:),mask_angle_row(1,:),'.r')
%}

%figure('units','normalized','outerposition',[0 0 1 1]);

for i = 1:size(row_range,2)-1 % 6 for histogram analysis (step size 5), 3 for step size 10
    for j = 1:size(col_range,2)- 1 %1 for size patch = 30 , center would go out of the image
        
        
        patch_rows = [row_range(i),row_range(i)+ (patch_width-1) ];
        patch_cols = [col_range(j),col_range(j)+ (patch_width-1) ];
        
        row_center_OD_total(i,j) = patch_rows(1) + round(patch_width/2); %patch size = 30
        col_center_OD_total(i,j) = patch_cols(1) + round(patch_width/2); %patch size = 30
        
        %if ( v1_region( row_center_OD_total(i,j),col_center_OD_total(i,j)) == 1 ) %v1_region is output of onoff_combine.m
        thresh_dipole_generation = .5 * sum(sum( ones(patch_width) ));%10 percent of patch size
        if ( sum(sum( v1_region_resized( patch_rows(1):patch_rows(2),patch_cols(1):patch_cols(2)))) > thresh_dipole_generation )
            
            
            input_patch = input_OD(patch_rows(1):patch_rows(2),patch_cols(1):patch_cols(2));
            
            patch_row_size = size(input_patch,1);
            patch_col_size = size(input_patch,2);
            
            [L_patch,L_num] = bwlabel(input_patch);
            
            num_stripes = 0;
            
            patch_npoint = zeros(1,L_num);
            patch_orientation = zeros(1,L_num);
            patch_thickness = zeros(1,L_num);
            
            
            J2 = zeros(patch_width,patch_width);
            
            if L_num > 0
                
                angle_line_vec = [];
                for i_L = 1:L_num
                    selected_region = L_patch == i_L;
                    %figure, imshow(selected_region)
                    
                    if ( sum(sum(selected_region))>3 ) % to remove noises
                        num_stripes = num_stripes + 1;
                        %[n_point,angle_line,thick] = horton_stripes_info4(selected_region);
                        [n_point,angle_line,thick] = horton_stripes_info2(selected_region);
                        
                        patch_npoint(i_L) = n_point;
                        patch_orientation(i_L) = mean(angle_line);
                        patch_thickness(i_L) = mean(thick);
                        
                        
                        angle_line_vec = cat(2,angle_line_vec,angle_line);
                    end
                    %title(sprintf('length = %.1f \n Orientation = %.1f \n Thickness = %.1f',patch_stripness(i_L),patch_orientation(i_L),patch_thickness(i_L)))
                end
                
                
                if num_stripes > 0
                    
                    n_point_patch = sum(patch_npoint);
                    length_mean_patch = n_point_patch / num_stripes;
                    patch_mean_orientation = sum( patch_orientation .* patch_npoint )/n_point_patch;
                    patch_mean_thickness = sum(patch_thickness .* patch_npoint )/n_point_patch;
                    
                    length_normalized = length_mean_patch/patch_row_size;
                    thick_normalized = patch_mean_thickness/patch_row_size;
                    
                    %measuring the difference of best match and selected
                    %patch morphologic parameters
                    diff_n_region = abs (( n_region_database_norm - ones(5,5,100)*num_stripes)/num_stripes );
                    diff_thickness = abs ( (thickness_database_norm - ones(5,5,100)*thick_normalized)/thick_normalized );
                    diff_length = abs ( (length_database_norm- ones(5,5,100)*length_normalized)/length_normalized );
                    
                    diff_total = diff_n_region + diff_thickness + diff_length;
                    [diff_value,ind_min_diff] = sort(diff_total(:));
                    
                    
                    
                    %*************orientation fourier*********************
                    thinned_image = bwmorph(input_patch,'thin',100);
                    resize_value2 = 1;
                    input_patch_resized = imresize(thinned_image,resize_value2);
                    
                    
                    y1 = fft2(input_patch_resized);
                    y2 = abs(fftshift(y1));
                    y3 = log(abs(y2)+1);
                    y3 = imrotate(y3,90); 
                    %figure,imagesc(y3)
                    
                    radial_integral = zeros(1,36);
                    for orient_count = 1:36 % 5 to 180 degrees with the step of 5
                        
                        r_temp = 0;
                        for radius_count = 1:ceil(col_input/2)
                            
                            r_temp = r_temp + y3(mask_angle_row(orient_count,radius_count),mask_angle_col(orient_count,radius_count));
                        end
                        
                        radial_integral(orient_count) = sum(r_temp);
                    end
                    
                    [~,ind_max] = max(radial_integral);
                    orienration_fourier = ind_max * 5;
                    %figure,plot(5:5:180,radial_integral)

 
                    attract_x_OD_total(i,j) = floor((ind_min_diff(1) - 725)/250)+3;
                    attract_y_OD_total(i,j) = attract_x_OD_total(i,j) ;
                    num_matrix =  mod ( mod((ind_min_diff(1) - 725),250) , 25);
                    database_col = round(num_matrix/5)+1;
                    database_row = mod(num_matrix,5);
                    repulse_x_OD_total(i,j) = database_row .* database_col .* attract_x_OD_total(i,j);
                    repulse_y_OD_total(i,j) = database_row  .* attract_x_OD_total(i,j);
                    
                    orientation_OD_total_fourier(i,j) = orienration_fourier;
                    
                    
                    if save_patches == 1
                        
                        figure(1)
                        clf
                        %figure('units','normalized','outerposition',[0 0 1 1]);
                        subplot(2,4,[1,2])
                        imshow( imresize(input_image,resize_value) )
                        hold on
                        rectangle('Position',[patch_cols(1),patch_rows(1),patch_row_size,patch_col_size],'EdgeColor','r','LineWidth',2)
                        title(sprintf('Patch Size = %.0f by %.0f ',patch_row_size,patch_col_size),'FontSize',14)
               
                        
                        subplot(243)
                        imagesc(y3)
                        
                        subplot(244)
                        plot(5:5:180,radial_integral,'LineWidth',2)
                        title(sprintf('Orientation Fourier\n %.0f degree ',orienration_fourier))
                        
                        
                        subplot(245)
                        imshow(input_patch)
                        title(sprintf('Horton Image \n N region = %.0f \n Orientation = %.1f \n Thick = %.4f \n Length = %.4f',num_stripes,patch_mean_orientation,thick_normalized,length_normalized))
                        
                        subplot(246)
                        nbins_orientation_contra_polar = 7;
                        orientation_contra_radian = angle_line_vec * (pi/180);
                        [frequency_contra,bin_edge_contra] = histcounts(orientation_contra_radian ,nbins_orientation_contra_polar,'Normalization','Probability'); %edge is the bar start and end points
                        edge1_contra = bin_edge_contra(1:end-1);
                        edge2_contra = bin_edge_contra(2:end);
                        bins_contra = (edge1_contra + edge2_contra)/2;
                        frequency_plot_contra = zeros(1,size(frequency_contra,2)+2);
                        bin_plot_contra = zeros(1,size(bins_contra,2)+2);
                        frequency_plot_contra(2:end-1) = frequency_contra;
                        bin_plot_contra(2:end-1) = bins_contra;
                        frequency_plot_contra(1) = 0;
                        frequency_plot_contra(end) = 0;
                        bin_plot_contra(1) = bins_contra(1) - abs(bins_contra(1)-bins_contra(2));
                        bin_plot_contra(end) = bins_contra(end) + abs(bins_contra(1)-bins_contra(2));
                        polar_histogram_mine(bin_plot_contra,frequency_plot_contra,'-','k');
                        title(sprintf('Mean Orientation \n every pixel \n %.2f degree ',patch_mean_orientation))
                        
                        
                        g=1; 
                        subplot(247)
                        imshow(output_result_database{ind_min_diff(g)})
                        title(sprintf('Simulation \n Diff(Percent) =%.2f \n N region = %.0f  \n Thick =%.4f \n Length = %.4f',diff_value(g)*100,n_region_database_norm(ind_min_diff(g)),...
                            thickness_database_norm(ind_min_diff(g)),length_database_norm(ind_min_diff(g)) ));
                        
                        
                        %J1 = imrotate(output_mask_database{ind_min_diff(g)},90+round(patch_mean_orientation),'crop');
                        J1 = imrotate(output_mask_database{ind_min_diff(1)},90+round(orienration_fourier),'crop');
                        
                        if ~isempty(J1)
                            J2 = J1(37:66,37:66);
                            
                            subplot(248)
                            imagesc(J2)
                            colormap('default')
                            title(output_title1_database{ind_min_diff(1)})
                            
                            mask_output{i,j} = J1;
                            mask_output_cropped{i,j} = J2;
                        end
                    end
                    
                    
                    counter1 = counter1 + 1;
                    
                    
                    
                    rx = round(repulse_x_OD_total(i,j)/1); % it could be devided by 1.5 for fitting in the image
                    ry = round(repulse_y_OD_total(i,j)/1); % it could be devided by 1.5 for fitting in the image
                    
                    x = -rx:.05:rx;
                    y1 = (ry/rx)*sqrt(rx^2 - x.^2);
                    y2 = -y1;
                    
                    theta = (90 + orienration_fourier) * (pi/180);
                    x1_rotated = cos(theta)*x + sin(theta)*y1;
                    y1_rotated = -sin(theta)*x + cos(theta)*y1;
                    
                    x2_rotated = cos(theta)*x + sin(theta)*y2;
                    y2_rotated = -sin(theta)*x + cos(theta)*y2;
                    
                    if (v1_region_resized(row_center_OD_total(i,j),col_center_OD_total(i,j)))
                        

                        %Removing huge dipoles for plotting purposes 
                        if condition(2) == 1 && repulse_x_OD_total(i,j) > 35
                            continue
                        end
                        
                        %Removing huge dipoles for plotting purposes 
                        if  repulse_x_OD_total(i,j) < 40 && repulse_x_OD_total(i,j) > 0
                            
                            axis(handles.axesImage);
                            hold on,
                            plot(x1_rotated + col_center_OD_total(i,j),y1_rotated + row_center_OD_total(i,j)  ,x2_rotated + col_center_OD_total(i,j) ,y2_rotated + row_center_OD_total(i,j) ,'Color',[180,0,180]/255,'LineWidth',2)
                            
                            
                            R = round(attract_x_OD_total(i,j)/1); % it could be devided by 1.5 for fitting in the image
                            xR = -R:0.5:R;
                            yR1 = sqrt(R^2 - xR.^2 );
                            yR2 = - yR1;
                            hold on,
                            plot(xR + col_center_OD_total(i,j),yR1 + row_center_OD_total(i,j) ,xR + col_center_OD_total(i,j) ,yR2 + row_center_OD_total(i,j),'Color',[255,183,0]/255,'LineWidth',2)
                            
                            ind_hist_plot(i,j) = 1;
                            shg
                            pause(.1)
                            
                        end
                        
                        
                    end
                    
                    
                end
            end
            

        end
        
    end
end








