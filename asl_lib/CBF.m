
function imgOut = CBF(Q,arg)

arg.null = 0;
opt.TR = 4000;
opt.lambda=0.9;
opt.alpha=0.85;
opt.LabDur=1800;% or Bolus Duration in Contrast/ASL card

opt.PLD=1800;
opt.InversionTime=opt.PLD + opt.LabDur;%(this is LabDur+PLD)
opt.T1blood=1650;
%BloodExp=0.3359110;
opt.TI1=800;% #nominal
opt.CBFthresh=400;
opt.M0thresh=10;
opt.CBFthresh = 400;

opt.scalingM0 = 1;

opt = getFiledsFromUsersOpt(opt,arg);

% #Compute CBF  by using white paper values for PCASL
imgOut = 1000/opt.scalingM0*(6000/2)*(opt.lambda/opt.alpha)*Q * 1/opt.T1blood * exp(opt.PLD/opt.T1blood)/(1-exp(-opt.LabDur/opt.T1blood));







% % # TopExp= exp(-TI2/T1blood)=exp(-1800/1650)=0.335911
% #
% #
% scalingM0=1 #for JJ ==1; for Siemens=10!    #this multiplies the M0map
% 
% TopExp=`bc -l <<< "e(${PLD}/${T1blood})"`
% echo $TopExp
% #2.97697918749497083586
% #BotExp=(1-exp(-LabDur/T1blood))
% BotExp=`bc -l <<< "1-e(-${LabDur}/${T1blood})"`
% echo $BotExp
% #.59710967847086700159