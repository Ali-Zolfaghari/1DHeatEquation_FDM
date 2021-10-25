function [ U ] = SOR( A,B,Ws,MAXERROR )

n = length(B);

U = zeros(1,n);
Us = zeros(1,n);

Iteration = 1;
Residual = 1.0;
while (Residual > MAXERROR)

    for i = 1:n
        sum = 0.0;
        for j = 1:i-1
            sum = sum+A(i,j)*U(j);
        end
        for j = i+1:n
            sum = sum+A(i,j)*Us(j);
        end
        Us(i) = (B(i)-sum)/A(i,i);
    end
    Iteration = Iteration+1;
    Residual = norm(Us-U);
	U = Ws*U+(1.0-Ws)*Us;
	
end

end



