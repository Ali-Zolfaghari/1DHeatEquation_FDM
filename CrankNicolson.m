function [ A,B ] = CrankNicolson( T,x,dx,dt,M,Tl,Tr,alpha )

A = zeros(M+1,M+1);
B = zeros(M+1,1);

r = (alpha*dt)/(dx*dx);

% left b.c.
i = 1;
A(i,i) = 1.0;
B(i,1) = Tl;

% right b.c.
i = M+1;
A(i,i) = 1.0;
B(i,1) = Tr;

% center
for i = 2:M
    A(i,i+1) = -0.5*r;
    A(i,i-1) = -0.5*r;
    A(i,i) = 1.0+r;
    B(i,1) = T(i)+0.5*r*(T(i+1)-2.0*T(i)+T(i-1));
end

end
