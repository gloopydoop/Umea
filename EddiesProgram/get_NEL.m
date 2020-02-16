function [NEL,NEL_x,NEL_y] = get_NEL(rho,omega)
NEL = 99999999999999;

for x = 1:size(rho,2)
    for y = 1:size(rho,1)
        if omega(y,x) == 1 && rho(y,x) == 1
            r_ij = RB_global(x,y,rho,omega,NEL);
            if r_ij < NEL
                NEL = r_ij;
                NEL_x = x;
                NEL_y = y;
                disp(['yo fam, NEL = ',num2str(NEL)])
            end
        end
    end
end