function X = constructX(x, dim)
if dim == 2
    X = zeros(length(x),3); % dim*I by 3 matrix
    for i = 1:length(x)/2 % i to I
        x1 = x(2*i-1);
        x2 = x(2*i);
        X(2*i-1,1) = x1;
        X(2*i,2) = x2;
        X(2*i-1,3) = 0.5*x2;
        X(2*i,3)   = 0.5*x1;
    end
end
if dim == 3
    X = zeros(length(x),4); % dim*I by 4 matrix
    for i = 1:length(x)/3 % i to I
        x1 = x(3*i-2);
        x2 = x(3*i-1);
        x3 = x(3*i);
        X(3*i-2,1) = x1;
        X(3*i-1,2) = x2;
        X(3*i-2,3) = 0.5*x2;
        X(3*i-1,3) = 0.5*x1;
        X(3*i,4) = x3;
    end
end
end