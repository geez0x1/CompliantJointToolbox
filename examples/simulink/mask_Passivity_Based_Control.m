% Get mask state
maskObj = Simulink.Mask.get(gcb);

% Activate / Deactivate DOB
maskObj.Parameters(5).Enabled = maskObj.Parameters(4).Value;
dob_switch = OnOff2Bool(maskObj.Parameters(4).Value);

% Gather information about 
[Amat, Bmat, Cmat, Dmat, Imat, Rmat, Kmat] = rigid_gearbox(jointObj);

D = jointObj.d_gl;  % Intrinsic spring damping
Jm = Imat(1,1);   % Effective physical inertia
invJm = 1/Jm;