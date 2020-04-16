function x_filt = filter_function(x,beta,R,param)
param.r = R;
filterParam =  initFilter_2D(param,beta);

n = 1;
s{1,n} = filterParam.cascade{1}.G{1}(...
        filterParam.cascade{1}.f{1}(x))./...
        filterParam.cascade{1}.Ni{1}; 
    s{2,n} = filterParam.cascade{1}.G{2}(filterParam.cascade{1}.f{2}(...
                    filterParam.cascade{1}.g{1}(s{1,n})) )./ ...
                    filterParam.cascade{1}.Ni{2};
    s{2,n} = filterParam.cascade{1}.g{2}(s{2,n});
    
    x_filt = s{2,1};