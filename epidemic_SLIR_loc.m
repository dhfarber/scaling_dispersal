function [p,u,v,rem,dis] = epidemic_SLIR_loc(n,T,l,s,r,m,f,u_0,e,d,lambda,a,b,size)
%EPIDEMIC_SIR
%epidemic spread with removal of sporulating infections after i days of
%infectiousness
%size - size in meters of each compartment n
u = zeros(T,n,l); %latent array
v = zeros(T,n,d); %sporulating array
p = zeros(T,n); %uninfected plant array
rem = zeros(T,n); %removed (no longer sporulating) array
kernel = square_trial_loc(n,a,b,size);
%m sites per compartment
%lambda - weight from (0,1) of dispersal;
%(1-lambda) - weight of somatic growth
for i = 1:s:n
    p(1:l,i) = m;
end

%u_0 initial severity in compartment(s) f
for i=f
    p(1,i) = m - m*u_0;
    u(1,i,1) = u_0*m;
end

%growth for all time 1:T for all compartments 1:n
for t=1:T-1
    fprintf('Day %d\n',t)
    for i=1:n
        a = kernel(i,:); %row array produced by kernel_square
        k = min(1,(sum(a(1:n).*sum(v(t,:,:),3)*r))/m);
        p(t,i) = max(0, m - sum(u(t,i,:))-sum(v(t,i,:))-rem(t,i));
        g = growth_bev(sum(u(t,i,:))+sum(v(t,i,:)),e,m);

        u(t+1,i,1) = p(t,i)*(lambda*k+(1-lambda)*g);
        %growth of immature latent infections
        for j = 1:l-1
            u(t+1,i,j+1) = u(t,i,j);
        end
        if u(t,i,l)+sum(v(t,i,:))<m
            v(t+1,i,1) = u(t,i,l);
        else
            v(t+1,i,1) = m - sum(v(t,i,2:d));
        end
        for j = 1:d-1
            v(t+1,i,j+1) = v(t,i,j);
        end
        rem(t+1,i) = rem(t,i) + v(t,i,d);
    end         
end
dis = sum(v,3)+rem;
%plotting:
%format bank
%figure()
%hold on
%for i = 2:l:t
%    plot(.0254*size:.0254*60:n*size*.0254,dis(i,:));
%end        
%str=sprintf('%.2f m compartments: R_0: %.2f P_0: %.2f',size*.0254,round(r,2),round(u_0,2));
%title(str)