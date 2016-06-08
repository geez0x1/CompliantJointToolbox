function [ dimSelect ] = getDimSelectorMatrix( jointObj )
% GETDIMSELECTORMATRIX Creates a mask matrix to select the proper joint angles and velocities for a joint model from the joint bus.
%
% [ dimSelect ] = getDimSelectorMatrix( jointObj )
%
% jointObj is the object defining the model to use.
%
% dimSelect is a ['linear system order'/2 x 3]  matrix. The matrix can be used to select the joint angles relevant to the joint model out of a full 
% 			vector of joint values (motor, gear, link). Same story for joint velocities.
%
% Examples: 
% For a completely rigid joint model dimSelect would be [1 0 0].
% For a rigid gearbox joint model with torsion bar dimSelect would be [1 0 0; 0 0 1].
%

    % Get jointObj dynamics matrices
    [A, ~, C] = jointObj.getDynamicsMatrices();

    % Parameters
    nEl = 3;                    % Number of position variables in the joint bus
    halfSysOrder = size(A,2)/2; % Half of actual system order

    
    %% Now compute the dimension selector matrix, which extracts the required 
    % motionvariables from the joint bus. The information is already encoded in
    % the dynamical system C matrix. We just have to find the unique rows. We
    % do that by dropping the non-unique ones to zero.

    dimSelect = C(1:nEl, 1:halfSysOrder); % transpose of dimension selector matrix to 
                                      % be determined, based on the dynamical system C matrix

    selectMask = zeros(1,nEl);        % this will be the mask that drops non-unique ones to zero
    [~, idx] = unique(dimSelect,'rows','stable'); % find indexes of unique ones
    selectMask(idx) = 1;              % mark unique ones
    selectMask = repmat(selectMask,halfSysOrder,1); % duplicate mask rowwise to obtain a mask 
                                        % matrix of appropriate system order

    % Apply th emask
    dimSelect = dimSelect.'  .* selectMask; % now transpose the dimension selector and apply the mask

end

