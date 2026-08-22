function Ystar = constructY(ystar, group)

n = length(ystar); % number of vertices * 3
e1 = [1; 0; 0];
e2 = [0; 1; 0];

if group == "translation"
    Y_tstar = zeros(n, 3);
    for i=1:n/3
        y_istar = ystar(3*i-2:3*i);
        Y_tstar(3*i-2:3*i,1) = y_istar(1)*e1;
        Y_tstar(3*i-2:3*i,2) = y_istar(2)*e2;
        Y_tstar(3*i-2:3*i,3) = 0.5*(y_istar(2)*e1 + y_istar(1)*e2);
    end
    Ystar = Y_tstar;

elseif group == "helical"
    Y_hstar = zeros(n, 3);
    for i=1:n/3
        y_istar = ystar(3*i-2:3*i);
        if i == 1
            r = norm(y_istar); % effective radius of the tube is determined by y_1
        end
        Y_hstar(3*i-2:3*i,1) = y_istar - y_istar(2)*e2;
        Y_hstar(3*i-2:3*i,2) = y_istar(2)*e2;
        Y_hstar(3*i-2:3*i,3) = y_istar(2)*cross(e2,y_istar)/r;
    end
    Ystar = Y_hstar;
end

end