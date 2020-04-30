function [yk] = MakeGoodY()
input = load('orbitdeterm_finalproj_KFdata.mat');
ybad = input.ydata;
yk = zeros(36,length(ybad));

for i = 2:length(ybad)
    obs = ybad{i};
    NumObs = size(obs,2);
    for j = 1:NumObs
        SID = obs(4,j);
        rho_i = obs(1,j);
        rho_dot_i = obs(2,j);
        phi_i = obs(3,j);
        
        y_i = [rho_i;rho_dot_i;phi_i];
        
        yk(3*SID-2:3*SID,i) = y_i;
    end
end
end