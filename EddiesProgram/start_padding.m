function x_padded = start_padding(x,param)
% HARRY BEING SNEAKY AGAIN!
pad = param.rFactor + 1;
line_BCs = reshape(1-param.nullel,param.nely,param.nelx);
if param.boundary_treatment == 1
    background = zeros(param.sizey,param.sizex);
    x(find(line_BCs == 0)) = 0;
    background(pad+1:end-pad,pad+1:end-pad) = reshape(x,param.nely,param.nelx);
    x = reshape(background,(param.nelx +2*pad)*(param.nely +2*pad),1);
    % Corners
    %background(1:pad,1:pad) = x(1,1)*ones(pad,pad);                 %NW
    %background(1:pad,end-pad:end) = x(1,end)*ones(pad,pad);         %NE
    %background(end-pad:end,1:pad) = x(end,1)*ones(pad,pad);         %SW
    %background(end-pad:end,end-pad:pad) = x(end,end)*ones(pad,pad); %SW
    % Sides
    
elseif param.boundary_treatment == 2
    %filterParam.line_BCs = ones(param.nelx*param.nely,1);
    % We set the interior be zero
    x(find(line_BCs == 0)) = 0;
end

x_padded = x;