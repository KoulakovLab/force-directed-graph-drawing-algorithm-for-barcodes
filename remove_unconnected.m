function M = remove_unconnected( M0 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

ind = find(sum(M0,2));
M = M0(ind,:);
ind = find(sum(M,1));
M = M(:,ind);

end

