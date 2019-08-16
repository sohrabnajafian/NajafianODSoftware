function [intersections,dist] = thickness(edge,point,angle)
%edge : outer edge of input image 
%input_point : main point [row,col]
%angle : direction of line 
%intersections : the intersections of perpendicular line from the point to
%the edge 
%dist : thickness (distance between intersections)
dist = []; 
intersections = []; 

a = zeros(size(edge)); 
a(point(1),point(2)) = 1; 

if (angle > 0)
    angle = angle - 90;
else 
    angle = 90 + angle; 
end 

[x1,y1] = meshgrid(-100 : 100, -100: 100);  
y1 = - y1; 
matrix_angle = (180/pi) * atan(y1./x1); 
matrix_angle = (matrix_angle == -90)*(180) + matrix_angle; 

if angle > 80
    mask_angle1 = (matrix_angle  < angle+10) & (matrix_angle > angle-10) & y1>0 ; %vertical lines should be excluded??
    mask_angle2 = (matrix_angle  < angle+10) & (matrix_angle > angle-10) & y1<0 ; % range of 10 degrees
else
    mask_angle1 = (matrix_angle  < angle+10) & (matrix_angle > angle-10) & x1>0 ; %vertical lines should be excluded??
    mask_angle2 = (matrix_angle  < angle+10) & (matrix_angle > angle-10) & x1<0 ; % range of 10 degrees
end

%mask_angle1(51,51) = 1; 
%mask_angle2(51,51) = 1; 

a1 = filter2(mask_angle1,a);
a1 = imdilate(a1,ones(3)); 

a2 = filter2(mask_angle2,a);
a2 = imdilate(a2,ones(3));

  
cond1 = a1 & edge; 
cond2 = a2 & edge; 

[p1,p2] = find(cond1); 
[p3,p4] = find(cond2); 

%figure,imagesc(a1+2*edge)
%figure,imagesc(a2+2*edge)

dist1 = zeros(size(p1,1),1);
if (size(p1,1)>0)
    for i = 1 : size(p1,1)
        dist1(i) = (p1(i)-point(1))^2 + (p2(i)-point(2))^2;
    end
    [~,ind1] = sort(dist1);
    
    intersections(1,1)= p1(ind1(1));
    intersections(1,2)= p2(ind1(1));
    
else 
    intersections(1,1)= point(1);
    intersections(1,2)= point(2);
end


dist2 = zeros(size(p3,1),1);
if (size(p3,1)>0)
    for i = 1 : size(p3,1)
        dist2(i) = (p3(i)-point(1))^2 + (p4(i)-point(2))^2;
    end
    [~,ind2] = sort(dist2);
    
    intersections(2,1)= p3(ind2(1));
    intersections(2,2)= p4(ind2(1));
else 
    intersections(2,1)= point(1);
    intersections(2,2)= point(2);    
end

%since two points of outer edge are selected the overal value should be
%substracted by one 
dist = sqrt((intersections(1,1) - intersections(2,1))^2 + (intersections(1,2) - intersections(2,2))^2) - 1;
