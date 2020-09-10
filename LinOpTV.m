function [Dx,Dy] = LinOpTV(u)
%Input is just the vector so we know how to calculate size.
%First order TV regularizer. Dealing with square matrices.
L = round(sqrt(length(u)));
%Horizontal operator.
Bin = zeros(L^2,2);
d1 = [-ones(L*(L-1),1);zeros(L,1)];
d2 = ones(L^2,1);%Don't include zeros.
Bin(:,1) = d1;
Bin(:,2) = d2;
d = [0 L];
Dx = spdiags(Bin,d,L^2,L^2);

%Vertical Operator
Bin = zeros(L^2,2);
d3 = -ones(L^2,1);
d3(L:L:L^2) = 0;
d4 = ones(L^2,1);
d4(1:L:L^2) = 0;
Bin(:,1) = d3;
Bin(:,2) = d4;
d = [0 1];
Dy = spdiags(Bin,d,L^2,L^2);

end

