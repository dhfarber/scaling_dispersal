function k = kernel_spatially_explicit_loc(source, field, a,b,size)
%returns proportion of total new infections arising in each cell of a
%field, from a single source. calls inv_power
%size is length of each compartment of field in inches
%calculates kernel by single inch, then puts 'em in bins of 'size' inches
%z = pdf_lomax(7260);
z = inv_power_loc(2000000,a,b,size);
y = zeros(1,2000000);
for i = 1:2000000
    if mod(size,2) == 0 %test if bin size in inches is even
        a = abs(i-(source*size-size/2)); %subtract size/2 to center source in compartment
    else
         a = abs(i-(source*size)); 
    end
    y(i) = z(a+1);
end
y = y./sum(y(1,:));
k = zeros(1,field);
for i = 1:length(k)
    a = size*(i-1)+1;
    b = size*i;
    k(i) = sum(y(a:b));
end