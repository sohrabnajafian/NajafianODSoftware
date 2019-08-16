function [n_point,angle_line,thick,point_coordinate] = horton_stripes_info2(selected_region,debug,pixel2um)
% Morphological measurement without pre sorting of the input points 
% the first column of point_coordinate shows row number of a point and
% second column shows the number of column of a point 


if nargin == 1
    debug = 0; 
end

%a1 = imread('horton_visual_cortex.jpg');
%a2 = imread('horton.jpg');

%horton_gray = a1(49:360,395:958);
%figure,imshow(horton_gray);
%title('gray image')

%horton_binary = a2(777:1069,153:717);
%imwrite(horton_binary,'horton_binary_cropped.jpg')


thin_image = bwmorph(selected_region,'thin',100);
%thin_image = bwmorph(thin_image,'spur',10);

%[endpoints_row,endpoints_col] = find(bwmorph(thin_image,'endpoints')); 
%{
figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,3,1)

subplot(1,3,2)
imshow(selected_region)
subplot(1,3,3)
imshow(thin_image)
title('thinned image')
%}

[point_row,point_col] = find(thin_image);


n_point = size(point_row,1);
angle_line = zeros(1,size(point_row,1));
thick = zeros(n_point,1);
intersect1 = zeros(n_point,2);
intersect2 = zeros(n_point,2);

point_coordinate = zeros(n_point,2);

edge_boundary = zeros(size(thin_image,1),size(thin_image,2)); 
edge_boundary(:,1) = 1; 
edge_boundary(:,end) = 1; 
edge_boundary(1,:) = 1; 
edge_boundary(end,:) = 1; 
edge_boundary = edge_boundary & selected_region; 

edge_selected =  (imdilate(selected_region,ones(3)) - selected_region) | edge_boundary;

if debug == 1 
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(121)
    imshow(selected_region)
    subplot(122)
    imshow(thin_image)
end

for i = 1 : n_point
    
    if n_point > 2
        dist = ( point_row(i)-point_row ).^2 + ( point_col(i)-point_col ).^2;
        [sorted_dist,ind] = sort(dist);
        
        if ( abs(sorted_dist(2)-sorted_dist(3)) < 2) %it means they are  neighbourhood with each other (surrounding)
            angle_line(i) =  -(180/pi) * atan((point_row(ind(2))-point_row(ind(3)))/(point_col(ind(2))-point_col(ind(3)) ));
        else
            angle_line(i) =  -(180/pi) * atan((point_row(ind(1))-point_row(ind(2)))/(point_col(ind(1))-point_col(ind(2)) ));
        end
    end
    


    [intersect,thick1] = thickness(edge_selected,[point_row(i),point_col(i)],angle_line(i));
    
    if ~isempty(thick1)
        thick(i) = thick1;
        intersect1(i,:) = intersect(1,:);
        intersect2(i,:) = intersect(2,:);
        
        point_coordinate(i,:) = [point_row(i),point_col(i)];
    end
    
    if debug == 1 
        if (mod(i,5)==0) %for presentation purposes 
            subplot(121)
            title('Selected Region','FontSize',18)
            hold on
            %plot(intersect(1,2),intersect(1,1),'mo','MarkerFaceColor','m','MarkerSize',8)
            %plot(intersect(2,2),intersect(2,1),'mo','MarkerFaceColor','m','MarkerSize',8)
            %plot(point_col(i),point_row(i),'ro','MarkerFaceColor','r','MarkerSize',8)
            line([intersect(1,2),intersect(2,2)],[intersect(1,1),intersect(2,1)],'Color','r','LineStyle',':','LineWidth',2)
            
            subplot(122)
            title(sprintf('Mid Line \n Thickness = %.2f pixels \n Thickness = %.2f um',thick1,thick1*pixel2um),'FontSize',18)
            hold on
            %plot(intersect(1,2),intersect(1,1),'mo','MarkerFaceColor','m','MarkerSize',10)
            %plot(intersect(2,2),intersect(2,1),'mo','MarkerFaceColor','m','MarkerSize',10)
            plot(point_col(i),point_row(i),'ro','MarkerFaceColor','r','MarkerSize',8)
            
            pause
        end 
    end 
    
end

if (mean(angle_line)< 0)
    angle_line(angle_line == 90) = -90;
elseif (mean(angle_line) > 0)
    angle_line(angle_line == -90) = 90;
end

