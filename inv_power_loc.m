function y = inv_power_loc(n,a,b,size)
k = zeros(1,n);
if size==60
    for i = 1:n
       % k(i) = 8.51*(i+.5)^-2.33; %converted to 5' blocks
       k(i) = a*(i+.5)^-b;
    end
else %size==1
    for i = 1:n
       % k(i) = .0356*(i+0.14*39.37)^-2.503; %converted to inches
       k(i) = a*(i+6.3)^-b;
    end
end
y = k./sum(k(1,:));
