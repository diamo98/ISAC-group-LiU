a = 0;
for i=1:1000
    a = a + abs(precoderSSB_GoB(:,1)'*(1/sqrt(2).*( randn(M,M) + 1i*randn(M,M)))*precoderSSB_GoB(:,100))^2;
end

a/1000

%%
b=0;
for i=1:1000
    b= b+ abs(precoderSSB_GoB(:,1)'*(1/sqrt(2).*( randn(M,1) + 1i*randn(M,1))))^2;

end
b/1000

%%
c =0;


for i=1:1000
    h = 1/sqrt(2).*( randn(M,1) + 1i*randn(M,1));
    cha = precoderSSB_GoB(:,1)*precoderSSB_GoB(:,1)';
    c= c+ abs(h'*cha*h);
    %c = c + abs((h)'*precoderSSB_GoB(:,1))^2;
end
c/1000

%%
c =0;


for i=1:1000
    h = 1/sqrt(2).*( randn(M,1) + 1i*randn(M,1)); 
    h2 = precoderSSB_GoB(:,1)*2; 
    cha = precoderSSB_GoB(:,1) * precoderSSB_GoB(:,1)'; 
    %c= c + abs(h'*cha*h); 
    c = c + abs((h)'*h2)^2; 

end

c/1000
%%
c =0;


for i=1:1000
    h = 1/sqrt(2).*( randn(M,1) + 1i*randn(M,1));
    h2 = precoderSSB_GoB(:,1)*2;
    cha = precoderSSB_GoB(:,1) * precoderSSB_GoB(:,1)';
    c= c + abs(h'*cha*h);
    %c = c + abs((h)'*h2)^2;
    
end

c/1000


%%
c =0;
n_var = 1;

for i=1:1000
    h = sqrt(n_var/2).*( randn(M,1) + 1i*randn(M,1));
    %h2 = precoderSSB_GoB(:,1)*2;
    %cha = precoderSSB_GoB(:,1) * precoderSSB_GoB(:,1)';
    c= c + norm(h)^2;
    %c = c + abs((h)'*h2)^2;
    
end

c/1000


%%
i=1;
% received signal y at APr [Mxq]
MC=5000;
SSB_symbols =  qammod(randi([0,mod_order-1],MC,1),mod_order,UnitAveragePower=true);
Sens_symbols =  qammod(randi([0,mod_order-1],MC,1),mod_order,UnitAveragePower=true);
%y_received_nu = zeros(MC,1)
y_received_nu_pwr = 0; y_received_denom_pwr = 0;
for mc=1:MC
        rcs1 = sqrt(rcs_var/2).*(randn()+1i*randn());
        rcs2 = sqrt(rcs_var/2).*(randn()+1i*randn());
        rcs3 = sqrt(rcs_var/2).*(randn()+1i*randn()); % H_all_fullloss H_all
        y_received_nu = Sens_symbols(mc,1).*( rcs1.*H_all_fullloss(:,1:M,i)*precoder_vec_por(:,1,i) + ...
            rcs1.*H_all_fullloss(:,M+1:2*M,i)*precoder_vec_por(:,2,i) + ...
            rcs1.*H_all_fullloss(:,2*M+1:3*M,i)*precoder_vec_por(:,3,i));
        y_received_denom =  SSB_symbols(mc,1).*sqrt(beta_jt_all/2).*sqrt(rho_1(1)).*(randn(M,M)+1i*randn(M,M) + randn(M,M)+1i*randn(M,M) + randn(M,M)+1i*randn(M,M))*conj(precoderSSB_GoB(:,i)) + sqrt(noise_var/2).*(randn(M,1)+1i*randn(M,1)) ;%.*((sqrt(beta_jt_all/2)).*(randn(M,M)+1i*randn(M,M) + randn(M,M)+1i*randn(M,M)+ ...
           % randn(M,M)+1i*randn(M,M) )*conj(precoderSSB_GoB(:,i))  ) + sqrt(noise_var/2).*(randn(M,1)+1i*randn(M,1));
        %SCIR_sim(mc,i) = norm(V(:,i)'*y_received_nu(:,i,mc))^2 / norm(V(:,i)'*y_received_denom(:,i,mc))^2;
        y_received_nu_pwr = norm(V(:,i)'*y_received_nu)^2 + y_received_nu_pwr;
        y_received_denom_pwr = y_received_denom_pwr + norm(V(:,i)'*y_received_denom)^2;
end

% mc
%y_received_nu_pwr/MC
y_received_denom_pwr_sim = y_received_denom_pwr/MC
% expected
%precoder_temp = reshape(precoder_vec(:,:,i),[J*M,1]);
%rcs_var*abs( V(:,i)'*                                                                                                                                     H_all_fullloss(:,:,i)*precoder_temp )^2
y_received_denom_pwr_expt = beta_jt_all*3*rho_1(1) + noise_var

% rcs at numerator should be similar
% denom >> 

%%

MC=1000;pwr_acc=0;
for i=1:MC
    n_temp = sqrt(noise_var/2).*(randn(M,1)+1i*randn(M,1));
    pwr_acc = pwr_acc + norm(V(:,1)'*n_temp)^2;
end

pwr_acc/MC

noise_var


%%
MC=1000;
temp=0; V1 = V(:,1);
for i=1:MC
    
end
temp/MC
