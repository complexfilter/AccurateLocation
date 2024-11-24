function dis=seq_compare(x_base,y_base,x1,y1)
MM=zeros(length(x_base),length(x1));
for i=1:length(x_base)
    MM(i,:)=(sqrt((x_base(i)-x1).^2+(y_base(i)-y1).^2))';
end
dis=min(MM,[],2);
end