function [x_left,y_left,x_right,y_right,dis1,dis2]=location_return(I,ind_left_x,ind_left_y,ind_right_x,ind_right_y,level,index,flag_positive_or_negative_left,flag_positive_or_negative_right,threshold,tolerance)

if size(I,3)==3
I=double(I(:,:,2));
end
multi_processing_nocut

patch_ini=zeros(size(I,1),size(I,2),6);
for i=1:6
  patch_ini(:,:,i)=eval([level,num2str(i)]);  
end

[x_left,y_left]=scan_weight_update(ind_left_x,ind_left_y,patch_ini(:,:,index),index,flag_positive_or_negative_left,threshold,tolerance);
[x_right,y_right]=scan_weight_update(ind_right_x,ind_right_y,patch_ini(:,:,index),index,flag_positive_or_negative_right,threshold,tolerance);

I_dis=zeros(size(I,1),size(I,2),3);
I_dis(:,:,1)=I;
I_dis(:,:,2)=I;
I_dis(:,:,3)=I;

figure(1);imshow(uint8(I_dis));title("The location of both sides of edges.");
hold on;
for i=1:length(x_left)-1 
    plot([y_left(i) + 0.5, y_left(i+1) + 0.5], [x_left(i) + 0.5, x_left(i+1)+0.5], 'r-')
end
for i=1:length(x_right)-1 
    plot([y_right(i) + 0.5 , y_right(i+1) + 0.5], [x_right(i) + 0.5, x_right(i+1) + 0.5], 'r-')
end
hold off;

dis1=seq_compare(x_left,y_left,x_right,y_right);
dis2=seq_compare(x_right,y_right,x_left,y_left);

end