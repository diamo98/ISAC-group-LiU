function [H_RR] = gen_H_RR(N,Repeater_pos)
H_RR = zeros(N,N); c=3e8; fc = 15e9; lambda = c/fc;
onelink_beta = @(l1) (lambda^2./(4*pi.*l1).^2);
for i=1:(N-1) % i is row
    for j=(i+1):N % j is column
        dis_nn = norm(Repeater_pos(:,i)-Repeater_pos(:,j));
        beta_nn = onelink_beta(dis_nn);
        H_RR(i,j) = sqrt(beta_nn).*exp(-1i*2*pi*fc*dis_nn/c); 
        %exp(-1i*2*pi*fc*dis_nn/c)
        H_RR(j,i) = H_RR(i,j);
    end
end

end