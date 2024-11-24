function [x_final_final,y_final_final]=scan_weight_update(ind1,ind2,patch,band_flag,flag_positive_or_negative,threshold, tolerance)
% this code is written to use a flexible way to determine the stroke
% location starting with the seed point. 

global tol;
tol=tolerance;  
global thres
thres=threshold;


[M,N]=size(patch);

ind_x0=zeros(M*10,1);
ind_y0=zeros(N*10,1);

[x_seed,y_seed]=find_right_point(ind1,ind2,patch,flag_positive_or_negative);

if (x_seed==-1)&&(y_seed==-1)
    error('The Seed Point Has A SMALL MAGNITUDE!');
end

ind_x0(1)=x_seed;
ind_y0(1)=y_seed;

ind=1;
for i=1:max(M,N)  
    x=ind_x0(ind);
    y=ind_y0(ind);

    if abs(patch(round(x),round(y)))<thres
       break;
    end
  [x_next_update,y_next_update,break_flag_final]=next_one(x,y,patch,band_flag,1,flag_positive_or_negative);
  if break_flag_final
      break;
  end
   ind=ind+1;
   ind_x0(ind)=x_next_update;
   ind_y0(ind)=y_next_update;
   
end
 

ind_x0=ind_x0(1:ind);
ind_y0=ind_y0(1:ind);


ind_x1=zeros(M*10,1);
ind_y1=zeros(N*10,1);

ind_x1(1)=x_seed;
ind_y1(1)=y_seed;

ind=1;

for i=1:max(M,N)   
    
    x=ind_x1(ind);
    y=ind_y1(ind);
    if abs(patch(round(x),round(y)))<thres
       break;
    end
  [x_next_update,y_next_update,break_flag_final]=next_one(x,y,patch,band_flag,0,flag_positive_or_negative);
  if break_flag_final
      break;
  end
   ind=ind+1;
   ind_x1(ind)=x_next_update;
   ind_y1(ind)=y_next_update;
end
    
ind_x1=ind_x1(2:ind);
ind_y1=ind_y1(2:ind);

%%%% finished the dots depiction above. %%%%

%%%% we need to transfrom the depiction of the dots to a single lines and gurantee
%% the that y-direction is continuous. (kind of interpolation) 
%% Since the x direction is continous already.
% we consider that x_final is increasing and y_final is decreasing by
% default.

ind_s1=find((ind_x1>=1)&(ind_x1<=M)&(ind_y1>=1)&(ind_y1<=N));
ind_x1=ind_x1(ind_s1);
ind_y1=ind_y1(ind_s1);

ind_s0=find((ind_x0>=1)&(ind_x0<=M)&(ind_y0>=1)&(ind_y0<=N));
ind_x0=ind_x0(ind_s0);
ind_y0=ind_y0(ind_s0);

x_final=[ind_x1(end:-1:1);ind_x0(1:end-1)]; % in ascending order
y_final=[ind_y1(end:-1:1);ind_y0(1:end-1)]; % in ascending order
 

x_final_final=zeros(size(x_final));
y_final_final=zeros(size(y_final));
inddd=1;
for i=1:length(x_final_final)
    if abs(patch(round(x_final(i)),round(y_final(i))))>thres
        x_final_final(inddd)=x_final(i);
        y_final_final(inddd)=y_final(i);
        inddd=inddd+1;
    end
end

x_final_final=x_final_final(x_final_final~=0);
y_final_final=y_final_final(x_final_final~=0);

end


function [orien,break_flag]=detect_orien(x,y,patch,band_flag,flag_forward_backward)

    x1=floor(x);
    x2=ceil(x);
    y1=floor(y);
    y2=ceil(y);


if (abs(x-round(x))<1e-5)&&(abs(x-round(x))<abs(y-round(y)))
    if flag_forward_backward
lo=check_availability(x1,y1,patch,band_flag);
[orien,~,~]=single_block_parameter(patch(x1:x1+1,y1:y1+1),1);
    else
lo=check_availability(x1-1,y1,patch,band_flag);
[orien,~,~]=single_block_parameter(patch(x1-1:x1,y1:y1+1),1);
    end
else
    if flag_forward_backward
lo=check_availability(x1,y1,patch,band_flag);
[orien,~,~]=single_block_parameter(patch(x1:x1+1,y1:y1+1),1);
    else
lo=check_availability(x1,y1-1,patch,band_flag);
[orien,~,~]=single_block_parameter(patch(x1:x1+1,y1-1:y1),1);
    end
end
   if lo    
    break_flag=0;
   else
    break_flag=1;  orien=0;
   end
   break_flag=logical(break_flag);
end

function lo=check_availability(x,y,patch,band_flag)
global tol;
[M,N]=size(patch);

    if (y>N-1)||(x>M-1)||(x<1)||(y<1)
        lo=0;
        lo=logical(lo);
        return;
    end

    tempp=patch(x:x+1,y:y+1);
    if ~isempty(tempp(abs(tempp)==0))
        lo=0;
        lo=logical(lo);
        return;
    end

    [orien,~,delta]=single_block_parameter(patch(x:x+1,y:y+1),1);
    if isreasonable(orien,band_flag)&&(abs(delta)<tol) 
        %% it seems to me that if we apply the code to B level band, we had better use the tolerance value bigger than 0.02.
        lo=1;
    else 
        lo=0;
        lo=logical(lo);
        return;
    end
    lo=logical(lo);
end

function [x_new_final,y_new_final]=find_right_point(x,y,patch,flag_positive_or_negative)

global thres
if abs(patch(x,y))<thres
    x_new_final=-1;
    y_new_final=-1;
    return;
end


plane=unwrap2Dphase(patch(x-1:x+1,y-1:y+1));
plane=plane';
vector=plane([2;4;5;6;8]);

if flag_positive_or_negative==1
    vector_minus_center=vector-(pi/2);
else
    vector_minus_center=vector-(-pi/2);  % disucuss the possibility that vector=pi/2 or -pi/2) 
                                        % I this the possibility of being exactly pi/2 or -pi/2 is minuscule.
end
if vector_minus_center(3)<-pi
    vector_minus_center=vector_minus_center+2*pi;
    vector=vector+2*pi;

elseif vector_minus_center(3)>pi
    vector_minus_center=vector_minus_center-2*pi;
    vector=vector-2*pi;
end

sign_vector_diff=abs(sign(vector_minus_center)-sign(vector_minus_center(3)));
        %disp([x,y]);

if isempty(find(sign_vector_diff~=0))
    
    dd=find(min(abs(vector_minus_center))==abs(vector_minus_center));dd=dd(1);
    
    if dd~=3
    while isempty(find(sign_vector_diff~=0))
        if dd==1
           x_update=x-1;
           y_update=y;
        elseif dd==2
           x_update=x;
           y_update=y-1; 
        elseif dd==4
           x_update=x;
           y_update=y+1;
        elseif dd==5
           x_update=x+1;
           y_update=y;
        end     

    plane=unwrap2Dphase(patch(x_update-1:x_update+1,y_update-1:y_update+1));
    plane=plane';
    vector=plane([2;4;5;6;8]);

if flag_positive_or_negative==1
    vector_minus_center=vector-(pi/2);
else
    vector_minus_center=vector-(-pi/2);  % disucuss the possibility that vector=pi/2 or -pi/2) 
                                        % I this the possibility of being exactly pi/2 or -pi/2 is minuscule.
end
if vector_minus_center(3)<-pi
    vector_minus_center=vector_minus_center+2*pi;
    vector=vector+2*pi;

elseif vector_minus_center(3)>pi
    vector_minus_center=vector_minus_center-2*pi;
    vector=vector-2*pi;
end

sign_vector_diff=abs(sign(vector_minus_center)-sign(vector_minus_center(3)));
       x=x_update;
       y=y_update;  
    end
    else        
        x_new_final=-1;
        y_new_final=-1;
        return;
    end  
end    
        
    mm=length(find(sign_vector_diff~=0));
    x_new=zeros(length(mm),1);
    y_new=zeros(length(mm),1);
   ind=find(sign_vector_diff~=0);
    for i=1:mm
        if ind(i)==1
           x_new(i)=x-abs(vector_minus_center(3))/(abs(vector_minus_center(3))+abs(vector_minus_center(1)));
           y_new(i)=y;
        elseif ind(i)==2
           x_new(i)=x;
           y_new(i)=y-abs(vector_minus_center(3))/(abs(vector_minus_center(3))+abs(vector_minus_center(2))); 
        elseif ind(i)==4
           x_new(i)=x;
           y_new(i)=y+abs(vector_minus_center(3))/(abs(vector_minus_center(3))+abs(vector_minus_center(4)));
        elseif ind(i)==5
           x_new(i)=x+abs(vector_minus_center(3))/(abs(vector_minus_center(3))+abs(vector_minus_center(5)));
           y_new(i)=y;
        end       
    end   
    dis=(x_new-x).^2+(y_new-y).^2;
    indd=find(min(dis)==dis);indd=indd(1);
    x_new_final=x_new(indd);
    y_new_final=y_new(indd);
end
 
function [x_next_update,y_next_update]=location_middle(x_next,y_next,patch,flag_positive_or_negative)
c1=floor([x_next,y_next]); 
c2=ceil([x_next,y_next]); 
phase1=angle(patch(c1(1),c1(2)));
phase2=angle(patch(c2(1),c2(2)));
phase_vector=unwrap([phase1,phase2]);

if flag_positive_or_negative
    phase=pi/2;
else
    phase=-pi/2;
end

if (phase_vector(1)-phase<-pi)&&(phase_vector(2)-phase<-pi)
    phase_vector=phase_vector+2*pi;

elseif (phase_vector(1)-phase>pi)&&(phase_vector(2)-phase>pi)
    phase_vector=phase_vector+2*pi;
end

phase1=phase_vector(1);
phase2=phase_vector(2);

    if (phase1-phase)*(phase2-phase)<0
        offset=abs((phase1-phase))/(abs(phase1-phase)+abs(phase2-phase));
        %offset1=abs((phase1-phase))/(abs(phase1-phase)+abs(phase2-phase));
        if isequal(c1(1),c2(1))
        x_next_update=c1(1);
        y_next_update=c1(2)+offset;
        elseif isequal(c1(2),c2(2))
        x_next_update=c1(1)+offset;
        y_next_update=c2(2);    
        end        
    else
        if abs(phase1-phase)>abs(phase2-phase)
            x=c2(1);y=c2(2);
        else
            x=c1(1);y=c1(2);
        end
        [x_next_update,y_next_update]=find_right_point(x,y,patch,flag_positive_or_negative);
    end

end


function [x_next,y_next]=next_coordinate_return(x,y,orien,flag_forward_backward)
%% this code automatically computes the next point on the grid. 
%% the unit of orien is degree, not angle. 
%% the index of x1,x2,x3 is depending on the counter-clockwise direction

     coordinates=zeros(3,2);
     
if flag_forward_backward
    
if  (abs(x-round(x))<1e-5)&&(abs(x-round(x))<abs(y-round(y)))
     x1=x+(floor(y)-y)/tan(orien*pi/180);y1=floor(y);
     dis1=sqrt(abs((x-x1)^2+(y-y1)^2)); coordinates(1,1)=x1;coordinates(1,2)=y1;
     x2=x+1;y2=tan(orien*pi/180)+y;
     dis2=sqrt(abs((x-x2)^2+(y-y2)^2)); coordinates(2,1)=x2;coordinates(2,2)=y2;
     x3=x+(ceil(y)-y)/tan(orien*pi/180);y3=ceil(y);
     dis3=sqrt(abs((x-x3)^2+(y-y3)^2)); coordinates(3,1)=x3;coordinates(3,2)=y3;
     dis=[dis1,dis2,dis3];

     coordinates_new=coordinates(coordinates(:,1)>x,:);
     dis_new=dis(coordinates(:,1)>x);
%      if min(dis_new)<1e-3
%          ind=find(dis_new==min(dis_new));
%          dis_new(ind)=[];
%          coordinates_new(ind,:)=[];
%      end
     n=find(dis_new==min(dis_new));n=n(1);
     %n
     x_next=coordinates_new(n,1);
     y_next=coordinates_new(n,2);
else
     x1=ceil(x);y1=tan(orien*pi/180)*(ceil(x)-x)+y;
     dis1=sqrt(abs((x-x1)^2+(y-y1)^2)); coordinates(1,1)=x1;coordinates(1,2)=y1;
     x2=1/tan(orien*pi/180)+x;y2=y+1;
     dis2=sqrt(abs((x-x2)^2+(y-y2)^2)); coordinates(2,1)=x2;coordinates(2,2)=y2;
     x3=-1/tan(orien*pi/180)+x;y3=y-1;
     dis3=sqrt(abs((x-x3)^2+(y-y3)^2)); coordinates(3,1)=x3;coordinates(3,2)=y3;
     dis=[dis1,dis2,dis3];
     coordinates_new=coordinates(coordinates(:,1)>x,:);
     dis_new=dis(coordinates(:,1)>x);
 
%      if min(dis_new)<1e-3
%          ind=find(dis_new==min(dis_new));
%          dis_new(ind)=[];
%          coordinates_new(ind,:)=[];
%      end

     n=find(dis_new==min(dis_new));n=n(1);
     
     x_next=coordinates_new(n,1);
     y_next=coordinates_new(n,2);
    
end

else
         
if  (abs(x-round(x))<1e-5)&&(abs(x-round(x))<abs(y-round(y)))
     x1=x+(ceil(y)-y)/tan(orien*pi/180);y1=ceil(y);
     dis1=sqrt(abs((x-x1)^2+(y-y1)^2)); coordinates(1,1)=x1;coordinates(1,2)=y1;
     x2=x-1;y2=-tan(orien*pi/180)+y;
     dis2=sqrt(abs((x-x2)^2+(y-y2)^2)); coordinates(2,1)=x2;coordinates(2,2)=y2;
     x3=x+(floor(y)-y)/tan(orien*pi/180);y3=floor(y);
     dis3=sqrt(abs((x-x3)^2+(y-y3)^2)); coordinates(3,1)=x3;coordinates(3,2)=y3;
     dis=[dis1,dis2,dis3];
     
     coordinates_new=coordinates(coordinates(:,1)<x,:);
     dis_new=dis(coordinates(:,1)<x);
     
%      if min(dis_new)<1e-isrea
%          ind=find(dis_new==min(dis_new));
%          dis_new(ind)=[];
%          coordinates_new(ind,:)=[];
%      end
     n=find(dis_new==min(dis_new));n=n(1);

     x_next=coordinates_new(n,1);
     y_next=coordinates_new(n,2);
      
else
     x1=floor(x);y1=tan(orien*pi/180)*(floor(x)-x)+y;
     dis1=sqrt(abs((x-x1)^2+(y-y1)^2)); coordinates(1,1)=x1;coordinates(1,2)=y1;
     x2=-1/tan(orien*pi/180)+x;y2=y-1;
     dis2=sqrt(abs((x-x2)^2+(y-y2)^2)); coordinates(2,1)=x2;coordinates(2,2)=y2;
     x3=1/tan(orien*pi/180)+x;y3=y+1;
     dis3=sqrt(abs((x-x3)^2+(y-y3)^2)); coordinates(3,1)=x3;coordinates(3,2)=y3;
     dis=[dis1,dis2,dis3];
     
     coordinates_new=coordinates(coordinates(:,1)<x,:);
     dis_new=dis(coordinates(:,1)<x);
 
     
%      if min(dis_new)<1e-3
%          ind=find(dis_new==min(dis_new));
%          dis_new(ind)=[];
%          coordinates_new(ind,:)=[];
%      end
     n=find(dis_new==min(dis_new));n=n(1);
     
     x_next=coordinates_new(n,1);
     y_next=coordinates_new(n,2);
    
end

end

end

function flag=isreasonable(ang,band_flag)
% the angle is the angle between the orientation and the row increasing
% direction. 
% when angle lies between 0 to 90, the edge points to the diagonal direction.
if numel(ang)==1
flag=single_return(ang,band_flag);
else 
    flag=zeros(size(ang));
    for i=1:length(flag)
        flag(i)=single_return(ang(i),band_flag);
    end
end
end

function flag=single_return(ang,band_flag)
if band_flag==4  %the constraints here are pretty loose
    %range=[-90 atan(1/3)*180/pi-90]; % to be discussed
    range=[-90 0];
elseif band_flag==5
    range=[-90 0];
elseif band_flag==6
    range=[-90 10];% allow some tolerance.
elseif band_flag==1
    range=[-10 90];
elseif band_flag==2
    range=[0 90]; 
elseif band_flag==3
    range=[0 90];
end

if (ang-range(1)>=0)&&(ang-range(2)<=0)
    flag=1;
else
    flag=0;
end
end
function [orien,w_hat,delta]=single_block_parameter(Band,flag)
% the value of orien is of degree. 是和x方向的夹角（x方向指的是矩阵的row ascending方向）.

[m,n]=size(Band);
orien=zeros(m-1,n-1);
w_hat=zeros(m-1,n-1);
delta=zeros(m-1,n-1);
for i=1:m-1
    for j=1:n-1
    PP=unwrap2Dphase(Band(i:i+1,j:j+1));
    [orien(i,j),w_hat(i,j),delta(i,j)]=block_parameter_return(PP(1),PP(2),PP(3),PP(4),Band(i:i+1,j:j+1),1,flag);
    end
end
end
function [orien,w_hat,delta]=block_parameter_return(ang1,ang2,ang3,ang4,patch,dis,flag)
% The function is written for the purpose of returning the desired
% parameters from each critically sampled(such as downsampled by 8 in B)
% block. 

% The parameters includes the w_hat, orientation and the residual of the
% flatness(ang1+ang3-ang2-ang4).

% the value of orien is of degree. 是和x方向的夹角（x方向指的是矩阵的row ascending方向）
[P,Q]=size(patch);
N=length(patch(:));
x=repmat((1:P)',Q,1);
y=zeros(N,1);
for i=1:Q
    y((i-1)*P+1:i*P)=i;
end
x=dis*x;
y=dis*y;

X=[x y ones(N,1)];
if flag==1
W=diag(abs(patch(:))).^2; % use the magnitudes as the weight to compute the parameters
else 
W=eye(N);
end
Y=[ang1 ang2 ang3 ang4];% supposed to be replaced.
Y=Y(:);
beta=(X'*W*X)\(X'*W*Y);
Z=zeros(P,Q);
for i=1:P
    for j=1:Q
        Z(i,j)=beta(1)*i+beta(2)*j+beta(3);
    end
end
delta=ang1+ang4-ang2-ang3;
x1=[ beta(1) beta(2) -1];
x2=[-beta(2) beta(1) 0];
x3=null([x1;x2]);
w_hat=x3(3)/sqrt(abs(x3(1)^2+x3(2)^2));
orien=atan(x2(2)/x2(1))*180/pi;
end
function [x_next_update,y_next_update,break_flag_final]=next_one(x,y,patch,band_flag,flag_forward_backward,flag_positive_or_negative)
    
    [M,N]=size(patch);
    if (x>=M-1) ||(y>=N-1)||(x<=1)||(y<=1)
       break_flag_final=1;x_next_update=-1;y_next_update=-1;
       return;
    end
    
    [orien,break_flag]=detect_orien(x,y,patch,band_flag,flag_forward_backward);
    if break_flag
       break_flag_final=1;x_next_update=-1;y_next_update=-1;
       return;
    end
    
   [x_next,y_next]=next_coordinate_return(x,y,orien,flag_forward_backward); % find the next point in a backward way
   [x_next_update,y_next_update]=location_middle(x_next,y_next,patch,flag_positive_or_negative);

   if (x_next_update==x)&&(y_next_update==y)
       if abs(x-round(x))<1e-5
          if flag_forward_backward
              if round(y)>=y
                  y_input=ceil(y);
              else
                  y_input=floor(y)-1;
              end
          [orien,break_flag]=detect_orien(x,y_input,patch,band_flag,flag_forward_backward);
          else
              if round(y)>=y
                  y_input=ceil(y);
              else
                  y_input=floor(y)-1;
              end
          [orien,break_flag]=detect_orien(x-1,y_input,patch,band_flag,flag_forward_backward);
          end
       else
          if flag_forward_backward
              if round(x)>=x
                  x_input=ceil(x);
              else
                  x_input=floor(x)-1;
              end
          [orien,break_flag]=detect_orien(x_input,y,patch,band_flag,flag_forward_backward);
          else
              if round(x)>=x
                  x_input=ceil(x);
              else
                  x_input=floor(x)-1;
              end
          [orien,break_flag]=detect_orien(x_input,y-1,patch,band_flag,flag_forward_backward);
          end
       end
          if break_flag || (x>M-1) ||(y>N-1)||(x<1)||(y<1)
            break_flag_final=1;x_next_update=-1;y_next_update=-1;
            return;
          end
          [x_next,y_next]=next_coordinate_return(x,y,orien,flag_forward_backward); % find the next point in a backward way
          [x_next_update,y_next_update]=location_middle(x_next,y_next,patch,flag_positive_or_negative);
   end
   
   if (x_next_update==-1)&&(y_next_update==-1)
    break_flag_final=1;
   else
    break_flag_final=0;
   end

end





