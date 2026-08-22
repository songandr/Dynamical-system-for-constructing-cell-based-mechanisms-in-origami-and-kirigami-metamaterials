function X = constructX(x)
    X = zeros(length(x),3); % 2I by 3 matrix
    for i = 1:length(x)/2
        x1 = x(2*i-1);
        x2 = x(2*i);
        X(2*i-1,1) = x1;
        X(2*i,2) = x2;
        X(2*i-1,3) = 0.5*x2;
        X(2*i,3)   = 0.5*x1;
    end
end