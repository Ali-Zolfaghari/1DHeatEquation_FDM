function [ A,B ] = DufortFrankel( T1,T2,x,dx,dt,M,Tl,Tr,alpha )

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
    A(i,i+1) = 0.0;
    A(i,i-1) = 0.0;
    A(i,i) = 1.0+2.0*r;
    B(i,1) = T1(i)+2.0*r*(T2(i+1)-T1(i)+T2(i-1));
end

end