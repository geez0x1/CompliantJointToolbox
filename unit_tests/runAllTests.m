clc;
setCJTPaths;

%% Run Tests
results_JointBuilder        = runtests('jointBuilderTest');
results_GenericJoint        = runtests('genericJointTest');
results_DataSheetGenerator  = runtests('dataSheetGeneratorTest');
results_bodeTest            = runtests('bodeTest');


%% Display Results
display('#####################################')
display('TEST SUMMARY')
display('#####################################')

results_JointBuilder
results_GenericJoint
results_DataSheetGenerator
results_bodeTest