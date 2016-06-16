function [ ] = deleteStateflowChart( )
% DELETESTATEFLOWCHART
% Deletes an orphaned StateflowChart from Simulink models

    thisBlock = gcb;
    funName = 'nonLinDynamicsFcn';
    fcnBlock = [thisBlock, '/', funName];

%     r = slroot;
%     scriptBlock = r.find('-isa','Stateflow.EMChart','path',fcnBlock);

%     if ~isempty(scriptBlock)
%         % For some reason Simulink does not recognize the EMChart
%         % as a proper Simulink object. It deletes it properly anyways,
%         % but issues a warning about this. Since delete_block is currently
%         % the only function known to properly delete the EMChart, we
%         % deactivate the warning here and reenable it afterwards as a
%         % workaround...
%         warning off
%         delete_block(scriptBlock);
%         warning on
%     end

    delete_line(thisBlock,[funName '/1'], ['Sum/1']);
    delete_line(thisBlock,[funName '/1'], ['nonLinearDynamicsScope/1']);
    delete_line(thisBlock,'Integrator/1',[funName '/1'])
    delete_block(fcnBlock)
    
end