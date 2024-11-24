clear;close all;clc;

% pick an integer coordinate that is near the left edge
ind_left_x=148;  
ind_left_y=219;   
% pick an integer coordinate that is near the right edge 
ind_right_x=152; 
ind_right_y=236;
% set the polarity of the left edge
flag_positive_or_negative_left=0;
% set the polarity of the right edge
flag_positive_or_negative_right=1;
threshold=20;
tolerance=0.08;
% The orientation of the edge will give rise to non-trival response on 
% B1 band. See Fig. 2(b) on ACCURATE EDGE LOCATION IDENTIFICATION BASED ON
% LOCATION-DIRECTED IMAGE MODELING for details. 
level='B';
index=1;
workingDir='data/';
num_frames=35;

x_seg1=130; % The horizontal scanning line's beginning cooridnate
x_seg2=180; % The horizontal scanning line's ending coordinate
x_slice=155; % The scanning line with in the beinning and middle to view the edge location variation


width=zeros(num_frames,1);
Y_RIGHT=zeros(num_frames,1);
for i=1:num_frames
    fprintf('Processing frame index: %d\n', i);
    filename = [sprintf('test%03d',i) '.png'];
    I=imread(fullfile(workingDir,filename)); % Reading frame.
    [x_left,y_left,x_right,y_right,dis1,dis2]=location_return(I,ind_left_x,ind_left_y,ind_right_x,ind_right_y,level,index,flag_positive_or_negative_left,flag_positive_or_negative_right,threshold,tolerance);
    ind1=find(x_right==x_seg1);
    ind2=find(x_right==x_seg2);
    Y_RIGHT(i)=y_right(x_right==x_slice);
    dis=dis2(ind1:1:ind2);
    width(i)=mean(dis);
end

figure;plot(width);title("The average width of the vessel across temporal frames");
% This reproduces Figure 4.3 on ACCURATE EDGE LOCATION IDENTIFICATION BASED ON
% LOCATION-DIRECTED IMAGE MODELING. 

figure;plot(Y_RIGHT + 0.5);title("The edge location aross temporal frames");
% This reproduces Fig. 7 on Determining Accurate Locations of Edges in 
% Natural Images: A Phase-Based, Nonparametric Framework. 




