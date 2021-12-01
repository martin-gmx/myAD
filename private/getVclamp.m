function [Vclamp00, addToName] = getVclamp(chooseVclamp)

Vhold = -80;

switch chooseVclamp
    case 0
        addToName = [];
        Vclamp00 = [100 Vhold; 1 30; 799 30; 1 -120; 599 -120; 1 0; 19 0; 1 -100; 399 -100; 1 Vhold; 300 Vhold; 1 -20; 300 -20; 1 Vhold; 300 Vhold];
    case 1
        addToName = '_longClamp';
        Vclamp00 = [];
        % activation & inactivation
        for i = 1:5;  Vclamp00 = [Vclamp00; 500 Vhold; 1 Vhold+i*20; 999 Vhold+i*20; 1 Vhold; 99 Vhold]; end
        % deactivation
        for i = 1:10; Vclamp00 = [Vclamp00; 1 60; 999 60; 1 -120+i*20; 999 -120+i*20]; end
        % reactivation
        for i = 1:9;  Vclamp00 = [Vclamp00; 1 60; 999 60; 1 -120+i*20; 9 -120+i*20; 1 40; 99 40; 1 -60; 199 -60]; end
    case 2
        addToName = '_midClamp';
        Vclamp00 = [];
        % activation & inactivation
        for i = 1:2;  Vclamp00 = [Vclamp00; 500 Vhold; 1 Vhold+i*20; 999 Vhold+i*20; 1 Vhold; 99 Vhold]; end
        for i = 5;  Vclamp00 = [Vclamp00; 500 Vhold; 1 Vhold+i*20; 999 Vhold+i*20; 1 Vhold; 99 Vhold]; end
        % deactivation
        for i = 1:3; Vclamp00 = [Vclamp00; 1 60; 999 60; 1 -130+i*30; 999 -130+i*30]; end
        % reactivation
        for i = 1:2;  Vclamp00 = [Vclamp00; 1 60; 999 60; 1 -120+i*20; 9 -120+i*20; 1 40; 99 40; 1 -60; 199 -60]; end
    case 3
        addToName = '_INaClamp';
        Vclamp00 = [];
        % activation & inactivation
        for i = 1:2;  Vclamp00 = [Vclamp00; 500 Vhold; 1 Vhold+i*20; 999 Vhold+i*20; 1 Vhold; 99 Vhold]; end
        for i = 5;  Vclamp00 = [Vclamp00; 500 Vhold; 1 Vhold+i*20; 999 Vhold+i*20; 1 Vhold; 99 Vhold]; end
        % deactivation
        for i = 1:3; Vclamp00 = [Vclamp00; 1 60; 999 60; 1 -130+i*30; 999 -130+i*30]; end
        % reactivation
        for i = 1:2;  Vclamp00 = [Vclamp00; 1 60; 999 60; 1 -120+i*20; 9 -120+i*20; 1 40; 99 40; 1 -60; 199 -60]; end
        % different holding potential
        for i = 1:3; Vclamp00 = [Vclamp00; 1 -190+i*30; 999 -190+i*30; 1 -20; 49 -20]; end
        % different time for reactivation (plateau time)
        for i = 1:3; Vclamp00 = [Vclamp00; 1 -120; 499 -120; 1 -20; 5+(i-1)*50 -20; 1 -120; 9 -120; 1 -20; 9 -20]; end
        % different time for reactivation (react time)
        for i = 1:3; Vclamp00 = [Vclamp00; 1 -20; 499 -20; 1 -120; 5+(i-1)*50 -120; 1 -20; 9 -20]; end
end
