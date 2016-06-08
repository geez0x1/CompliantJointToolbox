% clear all;
close all;

% Create jointBuilder
jb = jointBuilder;
jb.overwrite = 1;

%% Build various types of joint models

allParams = {'CE_Big';
    'WMBig10k_ds';
    'WMBig10k';
    'WMBig2300_ds'};
nPar = numel(allParams);

disp('Creating various joint types..');

for iPar = 1:nPar
    % No friction
    jb.buildJoint(allParams{iPar}, 'continuous_full_no_friction');
    jb.buildJoint(allParams{iPar}, 'continuous_output_fixed_no_friction');
    jb.buildJoint(allParams{iPar}, 'continuous_rigid_gearbox_no_friction');
    jb.buildJoint(allParams{iPar}, 'continuous_output_fixed_rigid_gearbox_no_friction');
    jb.buildJoint(allParams{iPar}, 'continuous_rigid_no_friction');
    
    % No nonlinearities
    jb.buildJoint(allParams{iPar}, 'continuous_full');
    jb.buildJoint(allParams{iPar}, 'continuous_output_fixed');
    jb.buildJoint(allParams{iPar}, 'continuous_rigid_gearbox');
    jb.buildJoint(allParams{iPar}, 'continuous_output_fixed_rigid_gearbox');
    
    % Coulomb friction - symmetric
    jb.buildJoint(allParams{iPar}, 'continuous_full', 'coulomb');
    jb.buildJoint(allParams{iPar}, 'continuous_output_fixed', 'coulomb');
    jb.buildJoint(allParams{iPar}, 'continuous_rigid_gearbox', 'coulomb');
    jb.buildJoint(allParams{iPar}, 'continuous_output_fixed_rigid_gearbox', 'coulomb');
    
    % Coulomb friction - asymmetric
    jb.buildJoint(allParams{iPar}, 'continuous_full', 'coulomb_asym');
    jb.buildJoint(allParams{iPar}, 'continuous_output_fixed', 'coulomb_asym');
    jb.buildJoint(allParams{iPar}, 'continuous_rigid_gearbox', 'coulomb_asym');
    jb.buildJoint(allParams{iPar}, 'continuous_output_fixed_rigid_gearbox', 'coulomb_asym');
    
    % Position-dependent viscous friction
    jb.buildJoint(allParams{iPar}, 'continuous_full', 'viscous_pos_dep');
    jb.buildJoint(allParams{iPar}, 'continuous_output_fixed', 'viscous_pos_dep');
    jb.buildJoint(allParams{iPar}, 'continuous_rigid_gearbox', 'viscous_pos_dep');
    jb.buildJoint(allParams{iPar}, 'continuous_output_fixed_rigid_gearbox', 'viscous_pos_dep');
    
    % Multiple nonlinear terms
    jb.buildJoint(allParams{iPar}, 'continuous_full', {'coulomb', 'viscous_asym'});
    jb.buildJoint(allParams{iPar}, 'continuous_output_fixed', {'coulomb_asym', 'viscous_asym'});
    jb.buildJoint(allParams{iPar}, 'continuous_rigid_gearbox', {'coulomb', 'viscous_asym'});
    jb.buildJoint(allParams{iPar}, 'continuous_output_fixed_rigid_gearbox', {'coulomb_asym', 'viscous_asym'});
end
