clear
clc
rng('shuffle')
inA = 0;
inB = 0;
for i = 1:10000
    rando = rand(2,1);
    ra = sqrt(1-rand(1,1)^2);
    rb = 1/(1+rand(1,1)^2);
    if rando(2,1)<=ra
        inA = inA + 1;
    end
    if rando(2,1)<=rb
        inB = inB + 1;
    end
    
    if i == 10
        pi10a = 4 * inA / 10
        err10a = abs(3.14159-pi10a)/(3.14159);
        pi10b = 4 * inB / 10
        err10b = abs(3.14159-pi10b)/(3.14159);
    end
    
    if i == 100
        pi100a = 4 * inA / 100
        err100a = abs(3.14159-pi100a)/(3.14159)
        pi100b = 4 * inB / 100
        err100b = abs(3.14159-pi100b)/(3.14159)
    end
    
    if i == 1E3
        pi1000a = 4 * inA / 1E3
        err1000a = abs(3.14159-pi1000a)/(3.14159)
        pi1000b = 4 * inB / 1E3
        err1000b = abs(3.14159-pi1000b)/(3.14159)
    end
    
    if i == 1E4
        pi10000a = 4 * inA / 1E4
        err10000a = abs(3.14159-pi10000a)/(3.14159)
        pi10000b = 4 * inB / 1E4
        err10000b = abs(3.14159-pi10000b)/(3.14159)
    end
end