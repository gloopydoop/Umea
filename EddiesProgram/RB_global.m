function r_ij = RB_global(x,y,rho,omega,NEL)

r_ij = 0;
switched = false;
for j = 1:size(rho,1)
    for i = 1:size(rho,2)
        if omega(j,i) == 1 && rho(j,i) == 1 && x~=i && y~=j
            r_local = RB_local(x,y,i,j,rho,NEL);
            if r_ij < r_local
                r_ij = r_local;
                switched = true;
            end
        end
    end
end

if switched == false
    r_ij = NaN;
end