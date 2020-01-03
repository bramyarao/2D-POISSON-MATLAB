function [P] = required_nodes(x,y,NS,ss)

% Required nodes for constructing the moment matrix is P
% Selecting the nodes around xI, till 4 times the ss
integer2 =1;
P = [];
for interger1 = 1:size(NS,1)
    x_m = NS(interger1,1);
    y_m = NS(interger1,2);
    if ((sqrt((x-x_m)^2 + (y-y_m)^2))<=ss)
        P(integer2,1) = x_m;
        P(integer2,2) = y_m;
        integer2 = integer2 + 1;
    else
        continue;
    end
end
