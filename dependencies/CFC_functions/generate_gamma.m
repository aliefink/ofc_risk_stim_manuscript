
for a = 1:200
    for p =1:50
        
        gam_shape = rand;
        gam_scale = rand;
        surr_matrix(a,p,:) = gamrnd(gam_shape,gam_scale,[1 200]);
    
    end 
end