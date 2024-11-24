function phase_matrix=unwrap2Dphase(patch,x_finally,y_finally)
phase_unwrap=zeros(size(patch));
phase_unwrap_transpose=zeros(size(patch)).';
for j=1:size(patch,1)
    phase_unwrap(j,:)=unwrap(angle(patch(j,:)));
end

for i=1:size(patch,2)
    phase_unwrap_transpose(i,:)=unwrap(angle(patch(:,i)));
end

phase_matrix=phase_unwrap;

[ind,~]=find(abs(patch)==max(abs(patch(:)))); % correction here. 

offset=phase_unwrap_transpose(:,ind)-phase_unwrap(ind,:).';

    

for j=1:size(patch,2)          
phase_matrix(:,j)=phase_unwrap_transpose(j,:)-offset(j);
end

if nargin==1
return;

elseif nargin==3
phase_unwrap=zeros(size(patch));
for j=1:size(patch,1)
    phase_unwrap(j,:)=unwrap(angle(patch(j,:)));
end

phase_matrix=phase_unwrap;

P=phase_unwrap(sub2ind(size(patch),x_finally,y_finally));
L=unwrap(angle(patch(sub2ind(size(patch),x_finally,y_finally))));
offset=P-L;

for j=1:length(x_finally)          
phase_matrix(x_finally(j),:)=phase_unwrap(x_finally(j),:)-offset(j);
end

% 
% phase_matrix2=phase_unwrap_transpose;
% 
% for j=1:length(y_finally)          
% phase_matrix2(:,y_finally(j))=phase_unwrap(:,y_finally(j))-offset(j);
% end
%       
end
end