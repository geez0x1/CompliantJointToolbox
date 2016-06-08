%%__________________________________________________________________
%% Get state-space representation of linear dynamics
function sys = getStateSpace(obj)
    [A, B, C, ~, ~, ~] = obj.getDynamicsMatrices();
    D = 0;
    sys = ss(A, B, C, D);
end