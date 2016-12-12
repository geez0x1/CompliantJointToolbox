%% Run Tests
results_JointBuilder       = runtests('jointBuilderTest');
results_GenericJoint       = runtests('genericJointTest');
results_DataSheetGenerator = runtests('dataSheetGeneratorTest');


%% Display Results
clc;

results_JointBuilder
results_GenericJoint
results_DataSheetGenerator
