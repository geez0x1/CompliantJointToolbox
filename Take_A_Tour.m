%% THE _COMPLIANT JOINT TOOLBOX_
%
% Note: This is the Code Ocean version of this example.
%
% The _Compliant Joint Toolbox_ (_CJT_) for MATLAB is being developed to ease the modelling and design of (compliant) robotic 
% joints. The goal is to provide rapid iteration of different models and control architectures by providing a number of 
% pre-built components with consistent interfaces, together with a set of tools to build new ones.
%
% This tour makes use of the MATLAB built-in publishing function. A separate run-script |Run.m| executes this script and
% produces a PDF with the console output. The authors encourage the reader to play with the contents of this script
% directly in their own MATLAB environment.
%
% Alternatively the reader can play online in the code capsule on |CodeOcean|. As a restriction to
% |CodeOcean|, the calls to |doc| and |edit| to illustrate the contents of certain files will not be displayed. However,
% feel free to open the files manually in the code section of the |CodeOcean| environment.
%
%% 

%% Getting started
% Clone the _CJT_ source files from GitHub using the URL:
%
% |https://github.com/geez0x1/CompliantJointToolbox|
%
% If you are not familiar with Git, these links might be a first resource to dive into its world:
%
% |http://rogerdudler.github.io/git-guide/|
%
% or
%
% |https://git-scm.com/book/en/v1/Getting-Started|
%
% Once you have obtained the toolbox source, open MATLAB and switch to the toolbox directory. In the toolbox directory
% Run | setCJTPaths| to add the _Compliant Joint Toolbox_ to your MATLAB search path.
setCJTPaths

%% Parameter Files
% The starting point for using the Compliant Joint Toolbox are parameter files. Parameter files are nothing but 
% m-scripts defining a struct named param. The class |genericJoint| already assigns default values to all parameters. 
% To see a full list parameters and methods, open the documentation of the |genericJoint| class by typing:
%
help genericJoint

%%
% A parameter script only requires to specify deviations from these default values. An example of such
% a parameter file is given in |example_params.m|.
%%
type example_params.m

%% Instantiate a jointBuilder
%
%
%
%% 
jb = jointBuilder;

%%
% The default build directory is |toolbox_root/build|. Let's change it to |../results|. We ommit the semicolon, to
% inspect a few more properties of the jointBuilder class. For the full class documentation call |doc jointBuilder|.
%
%%
jb.buildDir = ['..', filesep,'results']

 
%% Build a joint model class
% Here we reuse the parameters stored in |example_params.m|. We specify to use a model with rigid gearbox and chose
% Coulomb friction as well as asymmetric viscous frection as nonlinear dynamic effects. The model should incorporate
% electrical dynamics and eventually the generated class should be named |Example_Joint|.
%
%%
jb.buildJoint ( 'example_params' , ... parameters from example_params.m
'output_fixed_rigid_gearbox' , ... linear dynamics with rigid gearbox
{ 'coulomb' ,... {nonlinear
'viscous_asym'}, ... dynamics}
'electric_dyn',... electro-dynamics
'Example_Joint') % custom class name

%% Instantiate joint models
% Now that the |Example_Joint| class has been created, it is time to create a first instance of this class. 
%%

%Therefore, add the build directory to search path:
addpath(jb.buildDir)
% Then create a joint object. Ommit the semi-colon to learn more about the new |Example_Joint| instance:
exJoint = Example_Joint

%% Transfer Functions of the Example Joint 
% The |genericJoint| class builds upon the MATLAB core capabilities for linear system analysis via transfer functions and 
% state-space systems in continuous and discrete time domain. The benefit offered by the _Compliant Joint Toolbox_ is to 
% waive the need to manually equate and insert the model parameters into the corresponding built-in MATLAB functions (tf, 
% ss, etc.). Using the generated classes _Compliant Joint Toolbox_ offers direct access to the transfer functions and 
% state-space matrices in continuous and discrete time domain through a single line of code, independent of the actually 
% selected model structure.
%%

% Get all transfer functions 
exTF = exJoint.getTF;
% Look at the torque output (row index 7) 
% w. r. to the input current ( col index 1 )
exTF(7, 1)

% We obtain the same in the discrete 
% time domain with:
exTFd = exJoint.getTFd;
exTFd(7, 1)

% State-Space System the Example Joint 
% - In continuous time domain
exSS = exJoint.getStateSpace
% - In discrete time domain
exSSd = exJoint.getStateSpaceD

%% Symbolic Analysis
% With the Symbolic Math Toolbox™ installed, the Compliant Joint Toolbox offers to inspect the dynamics also in 
% symbolic form. This eases the analytical understanding of how individual parameters affect the dynamics. 
% Implementationwise, the Compliant Joint Toolbox offers the genericJoint methods makeSym and makeNum to convert 
% instances of joint models between numeric and symbolic representations. The following example looks again into the 
% transfer function of the previous example, but this time in symbolic form.
%%

% Convert joint to symbolic form
exJoint.makeSym; 

% Get all transfer functions 
exTF = exJoint.getTF ; 

% Look at the torque output (row index 7) w.r.t the input current ( col index 1) and pretty print the result.
pretty (exTF(7, 1) )

% Return to numeric form 
exJoint.makeNum;

%% Graphical Actuator Characteristics and Datasheet Generation
% The |datasheetGenerator| class implements the functionality to draw different analysis plots that permit inspection of
% the actuator's capabilities and performance. Therefore, the |datasheetGenerator| class is instantiated for a given 
% joint class. Then the drawing routines can be called. 
%%

% Instantiate a datasheetGenerator for the example joint:
dsg = dataSheetGenerator( exJoint );
% The default output directory of the dataSheetGenerator is the current directory. Let's also change this directory to 
% |../results|. Again, we ommit the semi-colon to learn more about the dataSheetGenerator properties.
dsg.outputDir = ['..', filesep, 'results']

% Draw a torque-speed diagram for the actuator in figure(1)
figure(1); clf;
dsg.drawTorqueSpeedCurve
print([dsg.outputDir,filesep,'torqueSpeedCurve'],'-dpdf','-fillpage') % save the figure in the output directory

% Draw an efficiency curve for the actuator in figure(2)
figure(2); clf;
dsg.drawEfficiencyCurve
print([dsg.outputDir,filesep,'efficiencyCurve'],'-dpdf','-fillpage') % save the figure in the output directory

% Draw a diagram with the thermal operation characteristics for the actuator in figure(3)
figure(3); clf;
dsg.drawThermalCharacteristics
print([dsg.outputDir,filesep,'thermalCharacteristics'],'-dpdf','-fillpage') % save the figure in the output directory

% Draw a torque-frequency curve for the actuator with locked output in figure(4)
figure(4); clf;
dsg.drawTorqueFrequencyCurveLocked
print([dsg.outputDir,filesep,'torqueFrequencyCurveLocked'],'-dpdf','-fillpage') % save the figure in the output directory

%% Datasheet Generation
% Provided, that a _LaTeX_ installation is present on the user’s computer, the |datasheetGenerator| class can assemble 
% the above plots into a datasheet PDF file summarizing the properties of for the respective actuator. The datasheet 
% generation is triggerd with a single command and the datasheet file is stored in the |datasheetGenerator| output 
% directory.
%% 
fName = dsg.generateDataSheet;

% Look at the output:
open(fName)

%% Simulation and Control
% The _Compliant Joint Toolbox_ features a Simulink library named |cjt_library|, which is located in the |toolboxroot/lib|
% directory. The block library includes mechanical and electrical subsystem models, state-of-the-art torque controllers
% as well as state and disturbance observers.
%
% All blocks are Simulink Real-Time compatible, so that they can be used in real-time applications and 
% deployed on physical target hardware systems. The library blocks make use of the joint model classes generated by 
% the |jointBuilder|. Their principal mask parameter is user-specified joint object or joint class name. The blocks 
% adapt their behaviour according to the dynamics and parameters specified in these derived joint classes.
%
% The Compliant Joint Toolbox comes with many examples on how to use these blocks.
%%

% Run an example of a simulation using the mechanical actuator model block and controller using a disturbance observer
% inside.

dispText = ...
{'-----'
'Note: Unfortunately this does not yet run properly in Code Ocean environment...'
'The User System version |Take_A_Tour.m| however will display a functioning GUI by running the command:'
'|Ex_12_Passivity_Based_Control_Plus_Disturbance_Observer_run|'
'-----'
};
displayFormattedText(dispText)


%% More Examples
% The _Compliant Joint Toolbox_ features many more examples in the |example| directory. A small graphical user interface
% helps to get an overview over existing examples.
%%

%%
% Note: Unfortunately the GUI cannot be displayed in Code Ocean. When trying to print the GUI figure into a file, it
% throws the error: "Printing of uicontrols is not supported on this platform." The User System version |Take_A_Tour.m|
% however will display a functioning GUI by running the command:
% |cjtExamples|
%%

dispText = ...
{'-----'
'Unfortunately the GUI cannot be displayed in Code Ocean. When trying to print the GUI figure into a file, it'
'throws the error: "Printing of uicontrols is not supported on this platform." The User System version |Take_A_Tour.m|'
'however will display a functioning GUI by running the command:'
'|cjtExamples|'
'-----'
};

displayFormattedText(dispText)


%% That's it!
% We hope it will be helpful for you. If it is, tell it to others. If you run into problems and find problems, please 
% tell it to us by opening an issue on GitHub! 
%
% |https://github.com/geez0x1/CompliantJointToolbox|
%
% Enjoy experimenting with the _Compliant Joint Toolbox_. 
%%

