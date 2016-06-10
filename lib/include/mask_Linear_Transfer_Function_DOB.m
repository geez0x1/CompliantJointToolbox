% modelFname = 'getLinearDOB_fromData_results.mat';
% loadData =1 ;
% jointObj = IROS_Motor_Model_Output_Fixed;

% Make sure local values are defined
Pc      = 0;
Q_td 	= 0;
PQ_td	= 0;

% If filename is given, assume the model and filters can be found in there
if (loadData)
    if (isempty(modelFname))
        error('Error: Model/filters data file not found!');
        return;
    else
        % Load file
        load(modelFname);
    end
else
	% Get linear DOB based on a model
    [Pc, Q_td, PQ_td] = getLinearDOB(jointObj, 2*pi*f_c, measIdx, doPlot);
end

% Discretize using Tustin transform
Q_td_tust	= c2d(Q_td, jointObj.Ts/2, 'tustin');
PQ_td_tust	= c2d(PQ_td, jointObj.Ts/2, 'tustin');
% Q_td_tust	= c2d(Q_td, jointObj.Ts, 'matched');
% PQ_td_tust	= c2d(PQ_td, jointObj.Ts, 'matched');