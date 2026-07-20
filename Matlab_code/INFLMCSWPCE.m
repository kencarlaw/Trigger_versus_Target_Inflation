%clear
load PCE.mat

% Create a partially specified model for estimation (NaNs for parameters to estimate)
PEst = NaN(2); % Estimate transition probabilities
mcEst = dtmc(PEst); 

% Create partially specified AR(1) submodels
mdlEst1 = arima(1,0,0); % AR(1) structure with NaN values for Constant, AR, Variance
mdlEst2 = arima(1,0,0); 
MdlEst = msVAR(mcEst, [mdlEst1; mdlEst2]);

% Create a fully specified model with initial values for the EM algorithm
P0 = [0.5 0.5; 0.5 0.5]; % Initial guess for P
mc0 = dtmc(P0);
mdl01 = arima('Constant', 0.5, 'AR', 0.9, 'Variance', 5);
mdl02 = arima('Constant', 0.5, 'AR', 0.5, 'Variance', 1);
Mdl0 = msVAR(mc0, [mdl02; mdl01]);

% Estimate the model using your actual data (replace 'y' with your data variable)
[EstMdl, SS, logL] = estimate(MdlEst, Mdl0, PCEPI(373:797));

% Display the estimated transition matrix
disp("Estimated Transition Matrix:");
disp(EstMdl.Switch.P);


mctest=dtmc(EstMdl.Switch.P);

mdltest=msVAR(mctest,[mdl01;mdl02]);

[ytest,~,sp]=simulate(mdltest,901);

fs=filter(mdltest,ytest);
ss=smooth(mdltest,ytest);
%figure
%hold on
%%plot(ss(:,2),"K");
%%plot(fs(:,2),"K",'LineStyle','--');
%plot(sp(:,1),"K");
%hold off

% Display the estimated submodel parameters for State 1
disp("Estimated State 1 Parameters:");
disp(EstMdl.Submodels(1));

% Display the estimated submodel parameters for State 2
disp("Estimated State 2 Parameters:");
disp(EstMdl.Submodels(2));

summarize(EstMdl)
