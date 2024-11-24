function I_filter=cconv2(hh1,hh2,I)
N=length(hh1);
M1=size(I,1); % ????
M2=size(I,2); % 列数
L1=size(I,1)+N-1;
L2=size(I,2)+N-1;
I_filter = conv2(hh1,hh2,I,'full');
I_filter(1:N-1,:)=I_filter(1:N-1,:)+I_filter((M1+1):L1,:);
I_filter(:,1:N-1)=I_filter(:,1:N-1)+I_filter(:,(M2+1):L2);
I_filter=I_filter(1:M1,1:M2);
I_filter=circshift(I_filter,[-N/2 -N/2]);
end
