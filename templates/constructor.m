%%__________________________________________________________________
%% Constructor
function this = %s(verbose, debug)
    %% Verbosity and debugging
    params.verbose		= 0;
    params.debug		= 0;
    if exist('verbose', 'var') params.verbose = verbose; end
    if exist('debug', 'var') params.debug = debug; end

    %% Model parameters
%s

    %% Sourced params and models
    params.('name')                 = '%s'; %% Joint name
    params.('paramName')            = '%s'; %% Parameter name
    params.('modelName')            = '%s'; %% Model name
    params.('nonlinearModelName')   = %s; %% Nonlinear model name(s)
    
    %% Build joint
    this = this@genericJoint(params);
end