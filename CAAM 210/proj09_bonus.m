function [x] = lu_solver(A,b)
[L,U,P] = lu(A);
c = P*b;
sub = 0;
for i = 1:length(c)
    y(i) = c(i) - sub;
    if i ~= length(c)
        sub = sub + L(i+1,i)*c(i);
    end
end
y = y';

x = zeros(length(c),1);
ssub = 0;
for j = 0:(length(c)-1)
    x(end-j) = (y(end-j) -ssub)/U(end-j,end-j);
    if j ~= length(c)-1
        ssub = ssub + U(end-j-1,end-j)*x(end-j);
    end
end


end