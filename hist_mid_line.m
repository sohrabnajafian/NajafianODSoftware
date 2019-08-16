 function [frequency,plot1] = hist_mid_line(data,nbins,plot_color)
%Plotting Histogram of input data with connecting the mid coordinates of
%histogram bars 
%Input nbins: number of bins that the data will be segregated into
%      data: input data 
%      line_color : color of plot line
%      spline_npoint : number of points to fit the spline 
%
%Output N : frequency of input (y axis of histogram)
%       plot1: spline curve for legend referencing 

    
    [frequency,bin_edge] = histcounts(data ,nbins,'Normalization','Probability'); %edge is the bar start and end points 
    %[frequency,bin_edge] = histcounts(data ,nbins); %edge is the bar start and end points 
 
    
    %bin_edge is the center postion of each bin
    %for histogram plotting the left and right coordinate of each bin is
    %calculated and the mean value is used for plotting 
    edge1 = bin_edge(1:end-1); 
    edge2 = bin_edge(2:end);
    bins = (edge1 + edge2) / 2;
    
    frequency_line = frequency; %for plotting purposes
    bin_edge_line = bins; %for plotting purposes
        
    outlier_ind = frequency == 0; 
    frequency(outlier_ind) = []; 
    bins(outlier_ind) = []; 


    %{
    spline_bins = linspace(round(min(bins)),round(max(bins)),spline_npoint); 
    curve_spline = spline(bins,frequency,spline_bins); 
    plot1 = plot(spline_bins,curve_spline,plot_color,'LineWidth',2);
    hold on 
    %}
    
    plot1 = plot(bin_edge_line,frequency_line,plot_color,'LineWidth',3);
    hold on
    plot(bins,frequency,[plot_color 'o'],'MarkerSize',10,'LineWidth',1.5);    
    
end 
