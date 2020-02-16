function r_local = RB_local(x,y,i,j,rho,NEL)

% find the smallest ball with (x,y) in it
r_first = ceil(norm([(x-i), (y-j)],2));

if r_first < NEL || NEL == 99999999999999
    r = r_first - 1;
    finished = false;

    while finished == false && r < NEL
        % search area
        r = r+1;
        N = min(size(rho,1),j+r);
        S = max(1,j-r);
        E = min(size(rho,2),i+r);
        W = max(1,i-r);

        for m = W:E
            for n = S:N
                if norm([(n-j), (m-i)],2)<=r && rho(n,m) == 0
                    finished = true;
                    break
                end
            end
        end

        % also break if we're at the wall of rho
        % NOTE: regarding the intersection with omega. I believe this works
        % currently for a rectangular domain. 
        % If something is cookie cut from the center i THINK you just need to
        % ensure that those cells are marked as design material in rho matrix
        if (N == size(rho,1) || S == 1 || E == size(rho,2) || W == 1)
            if finished == false
                r = r + 1;
            end
            finished = true;
        end
    end

    if r == r_first
        r_local = NaN;
    else
        r_local = r-1;
    end
     
else
    r_local = NaN;
end
            
        
