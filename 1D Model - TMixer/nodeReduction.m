function [isRealizable,rN] = nodeReduction(zeta,n)
    isRealizable = true;
    rN=1;
    
    if zeta(1)<=0
        isRealizable=false;
    end
    
    for i=3:2:2*n-1
        if zeta(i) <=0 || zeta(i-1)<=0
            isRealizable=false;
            break
        end
        rN=rN+1;
    end
end

