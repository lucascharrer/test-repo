function CoolingSim
%we start with a static Hamiltonian H0, and 
%the pulses are compiled on a fixed time grid. 
%showCrystalize();
%showCrystalizeMolecule();
%showNoCrystalizeMolecule();
%showCrystalizeMoleculeChange();
%show1dSqueezingMoleculeForce();
%show1dSqueezingMoleculeNoForce();
%show1dSqueezing();
%show1dSqueezUnsqueeze();
%show1dSqueezUnsqueezeImpulse();
%performanceCheck();
%tickleIon();
%tickleScan();
%tickleChirp();
%newTickleIon();
%newTickleScan();
%energyTickleScan();
%findFreq();
%testOptical();
%fullSequenceScan();
%animationTesting();
%recoolingScan();
%recooling();
%showSlapShot;
%flipRate();
%showpulseLaser();
%overdamping();
%pulseLaserSlapshot();
%showCWLaser();
%showBufferGas();
%bufferGasFlipRate();
%test();
%eraseDarkIons();
showChirp();
%chirpTickleScan();

%MUST ADD TO TICKLECLOUD TO ADD MORE THAN 18 IONS

function performanceCheck()
showCrystalize();
show2dSqueezing();
show1dSqueezing();

function show1dSqueezUnsqueezeImpulse()

ionCloud = initializeCloud(1,1,'looseZ');
pulseSet = addThermalize({},[0 4e-4],0.05);
pulseSet{end+1} = squeezePulse(ionCloud.atoms{1},'x',0,[4e-4 8e-4],10);
pulseSet{end+1} = squeezePulse(ionCloud.atoms{1},'x',pi,[10e-4 14e-4],10);
%pulseSet{end+1} = squeezePulse(ionCloud.atoms{1},'x',pi,[14e-4 18e-4],10);

finalCloud = evolveCloud(ionCloud,pulseSet,18e-4);

noiseAmp = 4e-22* sqrt(.005);
pulseSet{end+1} = absImpulse('x',[9e-4 9.001e-4],10000e-20,0);

finalCloudWithPush = evolveCloud(ionCloud,pulseSet,18e-4);
plotFields(finalCloudWithPush);
plotTrajectories(finalCloudWithPush,[5e-4 Inf]);
subplot(3,1,1);
p = plotJustZ(finalCloud,1,'x');
p.Color = 'k';
p.LineStyle = '-';
p.LineWidth = 1.0;
addText('Squeeze/unsqueeze with impulse at at 900 usec, black is without impulse');
subplot(4,1,2);
hold all;
plotEnergy(finalCloud,'k');
subplot(3,1,3);
plotXYZ(finalCloud,1)
addText('without heating at at 800 usec');

function show1dSqueezUnsqueeze()

ionCloud = initializeCloud(1,1,'looseZ');
pulseSet = addThermalize({},[0 4e-4],0.05);
pulseSet{end+1} = squeezePulse(ionCloud.atoms{1},'x',0,[4e-4 8e-4],10);
pulseSet{end+1} = squeezePulse(ionCloud.atoms{1},'x',pi,[10e-4 14e-4],10);
pulseSet{end+1} = squeezePulse(ionCloud.atoms{1},'x',pi,[14e-4 18e-4],10);

finalCloud = evolveCloud(ionCloud,pulseSet,18e-4);

noiseAmp = 4e-22* sqrt(.005);
pulseSet{end+1} = absNoisePulse('x',[8e-4 10e-4],noiseAmp * 4);

finalCloudWithPush = evolveCloud(ionCloud,pulseSet,18e-4);
plotFields(finalCloudWithPush);
plotTrajectories(finalCloudWithPush,[5e-4 Inf]);
subplot(3,1,1);
plotJustZ(finalCloud,1,'x');
addText('Squeeze/unsqueeze with heating at at 800 usec');
subplot(4,1,2);
hold all;
plotEnergy(finalCloud,'k');
subplot(3,1,3);
plotXYZ(finalCloud,1)
addText('without heating at at 800 usec');

function show1dSqueezing()

ionCloud = initializeCloud(1,1,'looseZ');
pulseSet = addThermalize({},[0 4e-4],0.05);
pulseSet{end+1} = squeezePulse(ionCloud.atoms{1},'x',0,[8e-4 1.6e-3],8);
finalCloud = evolveCloud(ionCloud,pulseSet,2.0e-3);

noiseAmp = 4e-22* sqrt(.005);
pulseSet{end+1} = absNoisePulse('x',[4e-4 5e-4],noiseAmp * 4);

finalCloudWithPush = evolveCloud(ionCloud,pulseSet,2.0e-3);
plotFields(finalCloudWithPush);
plotTrajectories(finalCloudWithPush,[0.3e-3 0.5e-3],[0.5e-3 0.8e-3],[0.8e-3 1.0e-3],[1.0e-3 1.5e-3],[1.5e-3 Inf]);
subplot(3,1,1);
plotJustZ(finalCloud,1,'x');
addText('with heating at at 500 usec');
subplot(4,1,2);
hold all;
plotEnergy(finalCloud,'k');
subplot(3,1,3);
plotXYZ(finalCloud,1)
addText('without heating at at 500 usec');

function show1dSqueezingMoleculeForce()
ionCloud = initializeCloud(2,1,'tightZ');
pulseSet = addThermalize({},[0 4e-4],0.05,6);
pulseSet{end+1} = absOpticalPulse('x',[3.15e-3 4.15e-3],1e5);
pulseSet{end+1} = squeezePulse(ionCloud.atoms{2},'x',(0.0 * pi),[4.2e-3 4.9e-3],6);
pulseSet{end+1} = squeezePulse(ionCloud.atoms{1},'x',(0.0 * pi),[6.9e-3 7.9e-3],32);
finalCloud = evolveCloud(ionCloud,pulseSet,8.0e-3,2.05e-3,1e4,5.45e-3,5.5e4); %Changing z potential only works for mixed clouds right now
%plotFields(finalCloud);
%plotTrajectories(finalCloud,[1.15e-3 1.65e-3],[1.65e-3 2.5e-3],[2.5e-3 3.0e-3],[3.0e-3 3.3e-3],[3.3e-3 Inf]);
%plotJustTrajectories(finalCloud,[7.9e-3 Inf]);
animated2dTrajectory(finalCloud,'z','x',[0.0e-3 3.0e-3],'fullSequenceForce1.gif');
animated2dTrajectory(finalCloud,'z','x',[3.0e-3 6.0e-3],'fullSequenceForce2.gif');
animated2dTrajectory(finalCloud,'z','x',[6.0e-3 8.0e-3],'fullSequenceForce3.gif');

function show1dSqueezingMoleculeNoForce()
ionCloud = initializeCloud(2,1,'tightZ');
pulseSet = addThermalize({},[0 4e-4],0.05,6);
pulseSet{end+1} = squeezePulse(ionCloud.atoms{2},'x',(0.0 * pi),[4.2e-3 4.9e-3],6);
pulseSet{end+1} = squeezePulse(ionCloud.atoms{1},'x',(0.0 * pi),[6.9e-3 7.9e-3],32);
finalCloud = evolveCloud(ionCloud,pulseSet,8.0e-3,2.05e-3,1e4,5.45e-3,5.5e4); %Changing z potential only works for mixed clouds right now
%plotFields(finalCloud);
%plotTrajectories(finalCloud,[1.15e-3 1.65e-3],[1.65e-3 2.5e-3],[2.5e-3 3.0e-3],[3.0e-3 3.3e-3],[3.3e-3 Inf]);
%plotJustTrajectories(finalCloud,[7.9e-3 Inf]);
animated2dTrajectory(finalCloud,'z','x',[0.0e-3 3.0e-3],'fullSequenceNoForce.gif');
animated2dTrajectory(finalCloud,'z','x',[3.0e-3 6.0e-3],'fullSequenceNoForce2.gif');
animated2dTrajectory(finalCloud,'z','x',[6.0e-3 8.0e-3],'fullSequenceNoForce3.gif');

function showCrystalize()
ionCloud = initializeCloud(10,10,'tightZ');
pulseSet = addThermalize({},[0 4e-4],0.05,0);
pulseSet = addFrictionPulse(pulseSet,[5e-4 15e-4]);
finalCloud = evolveCloud(ionCloud,pulseSet,1.5e-3);
%plotFields(finalCloud);
%plotTrajectories(finalCloud,[6e-4 Inf]);
%plotJustTrajectories(finalCloud,[5e-4 6e-4],[9e-4 Inf);
animated2dTrajectory(finalCloud,'z','x',[0e-3 1.5e-3],'crystal10.gif');
%plotEnergy(finalCloud);

function showCrystalizeMolecule()
ionCloud = initializeCloud(4,3,'tightZ');
pulseSet = addThermalize({},[0 4e-4],0.05,4);
finalCloud = evolveCloud(ionCloud,pulseSet,0.5e-3); 
%plotFields(finalCloud);
%plotJustTrajectories(finalCloud,[1.5e-3 1.6e-3]);
animated2dTrajectory(finalCloud,'z','x',[0.0e-3 0.5e-3],'DemoMoleculeCrystal.gif',0,1);

function showNoCrystalizeMolecule()
ionCloud = initializeCloud(2,1,'looseZ');
pulseSet = addThermalize({},[0 4e-4],0.05,4);
finalCloud = evolveCloud(ionCloud,pulseSet,1.0e-3);
%plotFields(finalCloud);
%plotJustTrajectories(finalCloud,[0.8e-3 1.0e-3]);
animated2dTrajectory(finalCloud,'z','x',[0e-3 1.0e-3],'moleculeNoCrystalv2.gif');

function showCrystalizeMoleculeChange()
ionCloud = initializeCloud(2,1,'tightZ');
pulseSet = addThermalize({},[0 4e-4],0.05,2);
pulseSet = addFrictionPulse(pulseSet,[1.75e-3 3.0e-3]);
pulseSet = addFrictionPulse(pulseSet,[3.5e-3 4.0e-3]);
finalCloud = evolveCloud(ionCloud,pulseSet,4.0e-3,1.0e-3,3e4,2.0e-3,5e4); %Changing z potential only works for mixed clouds right now
plotFields(finalCloud);
plotTrajectories(finalCloud,[1.15e-3 1.65e-3],[1.65e-3 2.1e-3],[2.8e-3 3.0e-3],[3.0e-3 3.8e-3],[3.9e-3 4.0e-3]);
animated2dTrajectory(finalCloud,'z','x',[1.3e-3 1.6e-3]);plotFluorescence(finalCloud);

function show2dSqueezing()
ionCloud = initializeCloud(2,2,'looseZ');
pulseSet = addThermalize({},[0 4e-4],0.1);
pulseSet{end+1} = squeezePulse(ionCloud.atoms{1},'x',0,[5e-4 1.6e-3],8);
finalCloud = evolveCloud(ionCloud,pulseSet,2.0e-3);
%plotFields(finalCloud);
plotTrajectories(finalCloud);

function tickleIon() % This function I was using to test some things, so it doesn't really do what it should right now.
ionCloud = initializeCloud(1,1,'customZ',100,100.1,50e3);
pulseSet = addThermalize({},[0 4e-4],0.05,250,1.5,-10e6);
% pulseSet = addFrictionPulse(pulseSet,[5.77e-4 10e-4],1,5e6);
% pulseSet{end+1} = absImpulse('z',[10e-4 14.0e-4],1e-19,0);
%pulseSet{end+1} = absImpulse('z',[3.505e-4 4.0e-4],1e-19,3e4);
%pulseSet = addFrictionPulse(pulseSet,[6.5e-4 7.5e-4]);
finalCloud = evolveCloud(ionCloud,pulseSet,100e-3);
%plotFields(finalCloud);
%plotTrajectories(finalCloud,[6e-4 8.5e-4]);
% animated2dTrajectory(finalCloud,'z','x',[0e-4 Inf],'testNewCool.gif');
figure(1);
%plotEnergy(finalCloud);
coord = pullTrajectory(finalCloud,'vx1');
% plot(finalCloud.times,coord);
%figure(2);
%fluor = idealFluorescence(finalCloud);
plot(finalCloud.times,coord);

% True intensity: ~ 0.8-1mW over a beam waist of 250 micrometers???

function tickleScan()
numSteps = 200;
finalTemps = zeros(1,numSteps);
tickleFrequencies = zeros(1,numSteps);
parfor i = 1:numSteps
    tickleFreq = (0.90 + i/500)*1e5;
    ionCloud = initializeCloud(1,1,'tightZ');
    pulseSet = addThermalize({},[0 4e-4],0.05,1.5);
    pulseSet{end+1} = absImpulse('x',[6e-4 12e-4],1e-19,tickleFreq);
    finalCloud = evolveCloud(ionCloud,pulseSet,12e-4);
    %plotFields(finalCloud.fieldSet);
    %plotTrajectories(finalCloud,[6e-4 8.5e-4]);
    finalTemps(i) = max(finalCloud.energies)/1.3e-23;
    tickleFrequencies(i) = tickleFreq;
end
plotFinalTemps(tickleFrequencies,finalTemps);

function newTickleIon()
firstCloud = initializeCloud(1,1,'tightZ');
pulseSet1 = addThermalize({},[0 4e-4],0.05,5);
ionCloud = evolveCloud(firstCloud,pulseSet1,10e-4);
pulseSet = {};
pulseSet = addFrictionPulse(pulseSet,[0 2.5e-4]);
pulseSet = addFrictionPulse(pulseSet,[7.5e-4 12.5e-4]);
pulseSet = addFrictionPulse(pulseSet,[17.5e-4 22.5e-4]);
pulseSet = addFrictionPulse(pulseSet,[27.5e-4 30e-4]);
pulseSet{end+1} = absImpulse('x',[2.5e-4 7.5e-4],1e-18,1e5);
pulseSet{end+1} = absImpulse('x',[12.5e-4 17.5e-4],1e-18,1.1169e5);
pulseSet{end+1} = absImpulse('x',[22.5e-4 27.5e-4],1e-18,1e5);
finalCloud = evolveCloud(ionCloud,pulseSet,30e-4);
%animated2dTrajectory(finalCloud,'z','x',[0e-4 Inf],'newTickle.gif');
plotEnergy(finalCloud);

function newTickleScan()
firstCloud = initializeCloud(1,1,'tightZ');
pulseSet1 = addThermalize({},[0 4e-4],0.05,5);
ionClouds(1) = evolveCloud(firstCloud,pulseSet1,20e-4);
numSteps = 200;
tickleFrequencies = zeros(1,numSteps);
finalTemps = zeros(1,numSteps);
for i = 1:numSteps
    tickleFreq = (0.5 + i/200)*1e5;
    pulseSet = {};
    pulseSet = addFrictionPulse(pulseSet,[0 0.5e-4]);
    pulseSet{end+1} = absImpulse('z',[0.5e-4 5e-4],1e-19,tickleFreq);
    tickleFrequencies(i) = tickleFreq;
    ionClouds(i+1) = evolveCloud(ionClouds(i),pulseSet,5e-4); 
    i
    finalTemps(i) = max(ionClouds(i+1).energies)/1.3e-23;
    %animated2dTrajectory(ionCloud,'z','x',[0 2e-4],'tickleScan.gif');
end
plotFinalTemps(tickleFrequencies,finalTemps);

function tickleChirp()
firstCloud = initializeCloud(1,1,'customZ',200,200.2,10e4);
pulseSet1 = addThermalize({},[0 4e-4],0.05,5);
ionCloud = evolveCloud(firstCloud,pulseSet1,20e-4);
numSteps = 50;
pulseSet = {};
phase = 0;
freqs = {};
for j = 0:numSteps
    freqs{end+1} = (11.0 - j * 2/50)*1e4;
end
for i = 0:numSteps
    tickleFreq = freqs{i+1};
    pulseSet{end+1} = absImpulse('z',[(i*5e-4 +1e-8) (i+1)*5e-4],1e-20,tickleFreq,phase);
    if i < numSteps
        oldPhase = phase;
        phase = -1 * ( ((i+1)*5e-4 * tickleFreq * 2 * pi) + oldPhase );
    end
end
pulseSet = addFrictionPulse(pulseSet,[0 (5*(numSteps+1)+0.5)*1e-4], 2e-21);
finalCloud = evolveCloud(ionCloud,pulseSet,(5*(numSteps+1)+0.5)*1e-4); 
plotEnergy(finalCloud);
% fieldSet = compileFields(pulseSet,0:1e-8:(5*(numSteps+1)+0.5)*1e-4);
% plot(fieldSet.Ez);


%TAKE FOURIER TRANSFORM OF NON-DRIVEN ION MOTION, FIND SECULAR FREQ

function energyTickleScan()
%parfor
    firstCloud = initializeCloud(3,2,'customZ',550,550.55,500e3);
    pulseSet1 = addThermalize({},[0 4e-4],0.05,5);
    ionCloud = evolveCloud(firstCloud,pulseSet1,20e-4);
    numSteps = 200;
    pulseSet = {};
    for i = 0:numSteps
        tickleFreq = (9.85 + i * 0.16/200)*1e4;
        pulseSet = addFrictionPulse(pulseSet,[i*5e-4 (5*i+0.5)*1e-4]);
        pulseSet{end+1} = absImpulse('z',[(5*i+0.5)*1e-4 (5*i+5)*1e-4],1e-20,tickleFreq);
    end
    % tickleFreq = 5e4;
    pulseSet = addFrictionPulse(pulseSet,[(numSteps+1)*5e-4 (5*(numSteps+1)+0.5)*1e-4]);
    % pulseSet{end+1} = absImpulse('x',[(5*(numSteps+1)+0.5)*1e-4 (5*(numSteps+1)+5)*1e-4],1e-20,tickleFreq);
    % pulseSet = addFrictionPulse(pulseSet,[(numSteps+2)*5e-4 (5*(numSteps+2)+0.5)*1e-4],3.041e-20 * 6.1);
    finalCloud = evolveCloud(ionCloud,pulseSet,(5*(numSteps+1)+0.5)*1e-4); 
    % animated2dTrajectory(finalCloud,'z','x',[0 50e-4],'newTickleScan.gif');
    plotEnergy(finalCloud);

% Default for just Sr+: 50 kHz      
% 2Sr+ & 1 molecule 104: 48.50 kHz   Shift: 1.50 kHz
% 2Sr+ & 1 molecule 96:  49.25 kHz   Shift: 0.75 kHz
% 2Sr+ & 1 molecule 120: 47.00 kHz   Shift: 3.00 kHz

% Increasing the axial frequency to 100 kHz, trying to increase the 25Hz gap
% 2Sr+ & 1 molecule 104: 96.675 kHz  Shift: 3.325 kHz
% 2Sr+ & 1 molecule 120: 94.200 kHz  Shift: 5.800 kHz 
% 2Sr+ & 1 molecule 96:  98.475 kHz  Shift: 1.525 kHz
% 2Sr+ & 1 molecule 89:  99.764 kHz  Shift: 0.236 kHz

% Increasing again to 500 kHz
% 2Sr+ & 1 molecule 89:   kHz  Shift:  kHz




function findFreq()
firstCloud = initializeCloud(1,0,'customZ',100,100.1,50e3);
pulseSet1 = addThermalize({},[0 0.4e-3],0.05,75,-75e6);
pulseSet1{end+1} = absImpulse('x',[29e-3 30e-3],1e-19,0);
ionCloud = evolveCloud(firstCloud,pulseSet1,30e-3);
finalCloud = evolveCloud(ionCloud,{},2.5e-3);
% 
% xSr = (88 / (88+89)) * pullTrajectory(finalCloud,'z1');
% xMol = (89 / (88+89)) * pullTrajectory(finalCloud,'z2');
% 
% x = xSr + xMol;

x = pullTrajectory(finalCloud,'x1');


times = finalCloud.times;
L = length(x);
if mod(L,2) == 1
    L = L - 1;
    x = x(1:L);
    times = times(1:L);
end
y = fft(x);
Y = abs(y);
p = Y(1:L/2 +1);
p(2:end-1) = 2*p(2:end-1);
deltaT = times(2)-times(1);
deltaF = 1/deltaT;
f = deltaF*(0:(L/2))/L;
plot(f,p);
xlim([0,5e6]);

%secular Sr x freq = 1.120e5 Hz apparently
%secular Molecule x freq = 0.936e5 Hz, y freq = 1.108e5 Hz


function testOptical()
ionCloud = initializeCloud(2,1,'tightZ');
pulseSet = addThermalize({},[0 4e-4],0.05,2);
%pulseSet = {};
pulseSet{end+1} = absOpticalPulse('x',[4e-4 10e-4],1e5);
finalCloud = evolveCloud(ionCloud,pulseSet,12e-4);
plotFields(finalCloud);
%plotTrajectories(finalCloud,[7e-4 8e-4],[8e-4 10e-4],[10e-4 12e-4]);

function fullSequenceScan()
numSteps = 40;
finalTemps = zeros(1,numSteps);
parfor i = 1:numSteps
    ionCloud = initializeCloud(2,1,'tightZ');
    pulseSet = addThermalize({},[0 4e-4],0.05,6);
    pulseSet{end+1} = absOpticalPulse('x',[3.15e-3 4.15e-3],1e5);
    pulseSet{end+1} = squeezePulse(ionCloud.atoms{2},'x',(0.0 * pi),[4.2e-3 4.9e-3],6);
    pulseSet{end+1} = squeezePulse(ionCloud.atoms{1},'x',(0.0 * pi),[6.9e-3 7.9e-3],32);
    finalCloud = evolveCloud(ionCloud,pulseSet,8.0e-3,2.05e-3,1e4,5.45e-3,5.5e4);
    finalTemps(i) = (finalCloud.energies(end)/1.3e-23);
end
figure('Position',[100   100   500   500],'Name',sprintf('Final Temperatures for 10 Runs'));
semilogy(finalTemps,'b*');
xlabel('Run Number');
ylabel('Final Temperature');

function animationTesting()
ionCloud = initializeCloud(2,1,'tightZ');
pulseSet = addThermalize({},[0 4e-4],0.05,1.5);
finalCloud = evolveCloud(ionCloud,pulseSet,10e-4,6e-4,3e4,8e-4,5e4);
%plotFields(finalCloud);
%plotTrajectories(finalCloud,[5e-4 Inf]);
%plotJustTrajectories(finalCloud,[9e-4 Inf]);
animated2dTrajectory(finalCloud,'z','x',[9e-4 10e-4]);

function recooling()
ionCloud = initializeCloud(4,3,'tightZ');
pulseSet = addThermalize({},[0 4e-4],0.05,1.75);
pulseSet = addThermalize(pulseSet,[7e-4 13e-4],0.05,0);
pulseSet = addFrictionPulse(pulseSet,[12e-4 17e-4]);
finalCloud = evolveCloud(ionCloud,pulseSet,17e-4);
%plotFluorescence(finalCloud);
animated2dTrajectory(finalCloud,'z','x',[0e-3 1.7e-3],'recooling.gif');

function recoolingScan()
initCloud = initializeCloud(4,3,'tightZ');
pulseSet1 = addThermalize({},[0 4e-4],0.05,1.75);
cooledCloud = evolveCloud(initCloud,pulseSet1,7e-4);
initZ = getZPosition(cooledCloud);
initZ
sameSpot = 0;
numTries = 12;
pulseSet2 = addThermalize({},[0 4e-4],0.005,0);
pulseSet2 = addFrictionPulse(pulseSet2,[3e-4 8e-4]);
parfor i = 1:numTries
    finalCloud = evolveCloud(cooledCloud,pulseSet2,8e-4);
    finZ = getZPosition(finalCloud);
    if finZ == initZ
        sameSpot = sameSpot + 1;
    end
    finZ
end
sameFraction = sameSpot / numTries;
sameFraction

function showSlapShot()
ionCloud = initializeCloud(3,2,'tightZ');
pulseSet = addThermalize({},[0 4e-4],0.05,4,9e-22);
pulseSet = slapShot(pulseSet,[26.00e-4 26.00005e-4],15e-16,0,-2.5e-16,[3]);
pulseSet = slapShot(pulseSet,[26.025e-4 26.02505e-4],0,15e-16,-15e-16,[3]);
pulseSet = addFrictionPulse(pulseSet,[26e-4 1000e-4],9e-22);
finalCloud = evolveCloud(ionCloud,pulseSet,100e-3);
animated2dTrajectory(finalCloud,'z','x',[2.4e-3 3.5e-3],'slapImpulseCoolMicroMotions.gif');
animated2dTrajectory(finalCloud,'z','x',[98e-3 100e-3],'slapImpulseLongCoolMicroMotions.gif');
%animated2dTrajectory(finalCloud,'x','y',[2.4e-3 3.5e-3],'slapImpulseCoolXY.gif');
plotFluorescence(finalCloud,[2.0e-3 100e-3]);
%plotTrajectories(finalCloud,[0e-4 5e-4],[5e-4 6.1e-4],[6e-4 6.5e-4],[9e-4 10e-4]);
%plotIndividualEnergies(finalCloud);

function flipRate()
ionCloud = initializeCloud(4,3,'customZ',1e5,1.05e5,2e4);
pulseSet1 = addThermalize({},[0 4e-4],0.05,0);
intermediateCloud = evolveCloud(ionCloud,pulseSet1,4e-4);
numTries = 20;
numZFlips = zeros(1,numTries);
finalTemps = zeros(1,numTries);
parfor i = 1:numTries
    pulseSet2 = addFrictionPulse({},[0 (7e-4 + i*0.05e-4)]);
    finalCloud = evolveCloud(intermediateCloud,pulseSet2,1e-2);
    initZPos = getZPosition(finalCloud,finalCloud.trajectories(1,:));
    for j = (round((i*1e-4)/1.697688e-8) + 1):length(finalCloud.trajectories)
        currentZPos = getZPosition(finalCloud,finalCloud.trajectories(j,:));
        if currentZPos ~= initZPos
            numZFlips(i) = numZFlips(i) + 1;
        end
        initZPos = currentZPos;
    end
    numZFlips(i) = numZFlips(i) * (1/(1e-2 - (i*1e-4)));
    finalTemps(i) = (finalCloud.energies(end)/1.3e-23);
end
plot(finalTemps,numZFlips,'*','MarkerSize',10,'MarkerEdgeColor','b');
xlabel('Final Temperature');
ylabel('Flips/Second');
plotFinalTemps(tickleFrequencies,finalTemps);

function showpulseLaser()
ionCloud = initializeCloud(1,1,'tightZ');
pulseSet = addThermalize({},[0 4e-4],0.05,0);
pulseSet = addFrictionPulse(pulseSet,[5e-4 25e-4]);
pulseSet{end+1} = pulseLaser([15e-4 22e-4],1,1050e-9,525e-9,2e-39,0,0,0.5e-4); %FIX THIS
pulseSet = slapShot(pulseSet,[16.0e-4 16.00005e-4],0,0,6.0e-16,[1]);
finalCloud = evolveCloud(ionCloud,pulseSet,2.5e-3);
animated2dTrajectory(finalCloud,'z','x',[1.5e-3 1.9e-3],'pulseLaserSlap.gif',0,1);

function overdamping()
ionCloud = initializeCloud(6,6,'customZ',100,100,2e4);
pulseSet = addThermalize({},[0 4e-4],0.05,1,3.3354e-19);
pulseSet = slapShot(pulseSet,[1.00e-4 1.00005e-4],5e-16,0,5e-16,[2]);
%pulseSet = addFrictionPulse(pulseSet,[4e-4 5e-4],9e-20);
finalCloud = evolveCloud(ionCloud,pulseSet,0.2e-3);
animated2dTrajectory(finalCloud,'z','x',[0e-3 0.2e-3],'overdamped6Ions.gif');

function pulseLaserSlapshot()
ionCloud = initializeCloud(4,4,'customZ',100,100.1,6e4);
pulseSet = addThermalize({},[0 4e-4],0.05,1);
pulseSet = addFrictionPulse(pulseSet,[4e-4 35e-4],9e-22);
pulseSet{end+1} = pulseLaser(33e-4,4e-3,532e-9,10e-6,2e-39,0,0,0e-6);
finalCloud = evolveCloud(ionCloud,pulseSet,3.500000e-3);
%finalVelocity = sqrt((finalCloud.vector(4))^2 + (finalCloud.vector(5))^2 + (finalCloud.vector(6))^2)
animated2dTrajectory(finalCloud,'z','x',[3.2e-3 3.5e-3],'PulseZigZagZagZigUnderdamped.gif');

function showCWLaser()
ionCloud = initializeCloud(1,1,'tightZ');
pulseSet = addThermalize({},[0 4e-4],0.05,1);
pulseSet = addFrictionPulse(pulseSet,[4e-4 34e-4],9e-20);
pulseSet{end+1} = CWLaser([33e-4 33.00005e-4],100,532e-9,2e-6,2e-39,0,0,2e-6);
finalCloud = evolveCloud(ionCloud,pulseSet,3.300005e-3);
finalVelocity = sqrt((finalCloud.vector(4))^2 + (finalCloud.vector(5))^2 + (finalCloud.vector(6))^2)
%animated2dTrajectory(finalCloud,'z','x',[3.2e-3 3.4e-3],'GaussianBeamSlapshot.gif');

function showBufferGas()
ionCloud = initializeCloud(4,4,'tightZ');
ionCloud.bufferGasBool = true;
ionCloud.bufferGasTime = [0e-4 inf];
ionCloud.bufferGasDensity = 1e14; %per cubic centimeter
ionCloud.bufferGasTemp = 4;
pulseSet = addThermalize({},[0 4e-4],0.05,5,3.041e-20 * 1);
finalCloud = evolveCloud(ionCloud,pulseSet,10e-4);
finalCloud.numCollisions
animated2dTrajectory(finalCloud,'z','x',[5e-4 10e-4],'BufferGasDemo.gif');
plotTrueTemp(finalCloud);

function bufferGasFlipRate()
ionCloud = initializeCloud(4,3,'tightZ');
ionCloud.bufferGasBool = true;
ionCloud.bufferGasTime = [0e-4 inf];
ionCloud.bufferGasDensity = 1e16; %per cubic centimeter
ionCloud.bufferGasTemp = 4;
pulseSet1 = addThermalize({},[0 4e-4],0.05,1);
intermediateCloud = evolveCloud(ionCloud,pulseSet1,4e-4);
intermediateCloud.bufferGasBool = true;
intermediateCloud.bufferGasTime = [0e-4 inf];
intermediateCloud.bufferGasDensity = 1e16; %per cubic centimeter
intermediateCloud.bufferGasTemp = 4;
numTries = 16;
numZFlips = zeros(1,numTries);
laserIntensity = zeros(1,numTries);
numCollisions = zeros(1,numTries);
parfor i = 1:numTries
    pulseSet2 = addFrictionPulse({},[0 1e-2], 3.041e-20 * 0.5 * i);
    finalCloud = evolveCloud(intermediateCloud,pulseSet2,1e-2);
    initZVector = getZVector(finalCloud,finalCloud.trajectories(1,:));
    for j = 1:length(finalCloud.trajectories)
        currentZVector = getZVector(finalCloud,finalCloud.trajectories(j,:));
        tempZFlips = 0;
        for k = 1:finalCloud.numAtoms
            if currentZVector(k) ~= initZVector(k)
                tempZFlips = tempZFlips + 1;
            end
        end
        numZFlips(i) = numZFlips(i) + tempZFlips;
        initZVector = currentZVector;
    end
    numZFlips(i) = numZFlips(i) * (1/(1e-2));
    laserIntensity(i) = i * 0.5;
    numCollisions(i) = finalCloud.numCollisions;
end
numCollisions
plot(laserIntensity,numZFlips,'*','MarkerSize',10,'MarkerEdgeColor','b');
xlabel('Laser Intensity (Multiple of Saturation Intensity)');
ylabel('Flips/Second');

function bufferGasFinalTemp
ionCloud = initializeCloud(4,3,'tightZ');
ionCloud.bufferGasBool = true;
ionCloud.bufferGasTime = [0e-4 inf];
ionCloud.bufferGasDensity = 1e16; %per cubic centimeter
ionCloud.bufferGasTemp = 4;
pulseSet1 = addThermalize({},[0 4e-4],0.05,1);
intermediateCloud = evolveCloud(ionCloud,pulseSet1,4e-4);
intermediateCloud.bufferGasBool = true;
intermediateCloud.bufferGasTime = [0e-4 inf];
intermediateCloud.bufferGasDensity = 1e16; %per cubic centimeter
intermediateCloud.bufferGasTemp = 4;
numTries = 16;
numZFlips = zeros(1,numTries);
laserIntensity = zeros(1,numTries);
numCollisions = zeros(1,numTries);
parfor i = 1:numTries
    pulseSet2 = addFrictionPulse({},[0 1e-2], 3.041e-20 * 0.5 * i);
    finalCloud = evolveCloud(intermediateCloud,pulseSet2,1e-2);
    initZVector = getZVector(finalCloud,finalCloud.trajectories(1,:));
    for j = 1:length(finalCloud.trajectories)
        currentZVector = getZVector(finalCloud,finalCloud.trajectories(j,:));
        tempZFlips = 0;
        for k = 1:finalCloud.numAtoms
            if currentZVector(k) ~= initZVector(k)
                tempZFlips = tempZFlips + 1;
            end
        end
        numZFlips(i) = numZFlips(i) + tempZFlips;
        initZVector = currentZVector;
    end
    numZFlips(i) = numZFlips(i) * (1/(1e-2));
    laserIntensity(i) = i * 0.5;
    numCollisions(i) = finalCloud.numCollisions;
end
numCollisions
plot(laserIntensity,numZFlips,'*','MarkerSize',10,'MarkerEdgeColor','b');
xlabel('Laser Intensity (Multiple of Saturation Intensity)');
ylabel('Flips/Second');

function eraseDarkIons()
ionCloud = initializeCloud(2,1,'customZ',100,100.1,50e3);
pulseSet = addThermalize({},[0 4e-4],0.05,1);
pulseSet = addFrictionPulse(pulseSet,[4e-4 5e-4],3.041e-20 * 6.1 * 2);
pulseSet{end+1} = absImpulse('x',[4.5e-4 4.75e-4],1e-17,0.936e5); 
pulseSet{end+1} = absImpulse('y',[4.75e-4 5.0e-4],1e-17,0.936e5);
pulseSet = addFrictionPulse(pulseSet,[5e-4 56e-4],3.041e-20 * 6.1 * 0.5);
finalCloud = evolveCloud(ionCloud,pulseSet,56e-4);
%animated2dTrajectory(finalCloud,'z','x',[0 Inf],'eraseDarkIons.gif');
plotJustTrajectories(finalCloud,[0e-4 4.5e-4],[4.5e-4 4.75e-4],[4.75e-4 5e-4],[5e-4 6e-4],[6e-4 Inf]);

function showChirp()
% MAKE SURE EVERYTHING IS ROBUST; SHOULD BE ABLE TO CHANGE PARAMETERS BY
% ~20%-30% AND STILL HAVE SEQUENCE WORK
initCloud = initializeCloud(2,1,'customZ',100,100.1,5e4);
pulseSet1 = addThermalize({},[0 4e-4],0.05,100,1.5,-100e6);
ionCloud1 = evolveCloud(initCloud,pulseSet1,40e-3);
pulseSet2 = addFrictionPulse({},[0 20e-3],1.5,-25e6);
ionCloud = evolveCloud(ionCloud1,pulseSet2,20e-3);
% Above is just cooling the ions down

pulseSetChirp = addFrictionPulse({},[0 7.5e-3],1.5,-5e6);
pulseSetChirp{end+1} = chirpPulse('z',[0 7.5e-3],[45e3 55e3],1e-21);

pulseSetOrbit = addFrictionPulse({},[0 1e-4],3.0,-100e6);
pulseSetOrbit{end+1} = absImpulse('x',[0.5e-4 0.75e-4],-1e-19,0.936e5); 
pulseSetOrbit{end+1} = absImpulse('y',[0.75e-4 1.0e-4],-1e-19,0.936e5);

intermediateCloud = evolveCloud(ionCloud,pulseSetOrbit,1e-4);

pulseSet3 = addFrictionPulse({},[0 60e-3],1.5,-25e6);
ionCloud2 = evolveCloud(intermediateCloud,pulseSet3,60e-3);


animated2dTrajectory(ionCloud2,'x','y',[59.5e-3 60.0e-3],'eraseDarkIons.gif');

totalFluor = idealFluorescence(ionCloud2);
Ftimes = linspace(0,0.2,length(totalFluor));
plot(Ftimes,totalFluor);

% Chirp with just Sr

% finalCloud1 = evolveCloud(ionCloud,pulseSetChirp,7.5e-3);
% finalCloud2 = evolveCloud(ionCloud2,pulseSetChirp,7.5e-3);


% fluorescence1 = idealFluorescence(finalCloud1,[0 Inf]);
% freqs1 = linspace(45,55,length(fluorescence1));
% plot(freqs1,fluorescence1);

% fluorescence2 = idealFluorescence(finalCloud2,[0 Inf]);
% freqs2 = linspace(45,55,length(fluorescence2));
% plot(freqs2,fluorescence2);



function chirpTickleScan()
initCloud = initializeCloud(1,1,'customZ',100,100.1,5e4);
pulseSet1 = addThermalize({},[0 4e-4],0.05,1.25);
ionCloud = evolveCloud(initCloud,pulseSet1,5e-4);
numSteps = 4;
freqRange = [90e3 140e3];

parfor i = 1:numSteps
    pulseSet = {}
    pulseSet = addFrictionPulse(pulseSet,[1e-4 75e-4],3.041e-20 * 1);
    pulseSet{end+1} = chirpPulse('x',[1e-4 75e-4],freqRange,1e-18);
    finalCloud = evolveCloud(ionCloud,pulseSet,75e-4);
    fluorescence = idealFluorescence(finalCloud);
end
totalFluor = zeros(1,length(fluorescence(1)));
for i = 1:numSteps
    totalFluor = totalFluor + fluorescence(i);
end
totalFluor = totalFluor / numSteps;
freqs = linspace(freqRange(1),freqRange(2),length(totalFluor));
plot(freqs,totalFluor);

%sequence: do a chirp of just bright ions, then a chirp of bright and dark
%ions



function pulseSet = makePulseSet(cloud)
%first make a pulseList describing all the fields.  This can be compiled
%into actual fields
    pulseSet = {};
    pulseSet = addThermalize(pulseSet,[0 4e-4],0.1);
    pulseSet{end+1} = squeezePulse(cloud.atoms{1},'x',0,[5e-4 1.6e-3],20);
 
function pulseSet = addThermalize(pulseSet,timeSpan,T,frictionMult,frictionAmp,delta)
if nargin < 4
    frictionMult = 1;
end
if nargin < 5
    frictionAmp = 1;     
end
if nargin < 6
    delta = 10e6;
end

noiseAmp = 4e-22* sqrt(T/10);
    pulseSet{end+1} = absFrictionPulse('x',(frictionMult*timeSpan),frictionAmp,delta);
    pulseSet{end+1} = absFrictionPulse('y',(frictionMult*timeSpan),frictionAmp,delta);
    pulseSet{end+1} = absFrictionPulse('z',(frictionMult*timeSpan),frictionAmp,delta);
    pulseSet{end+1} = absNoisePulse('x',timeSpan,noiseAmp );
    pulseSet{end+1} = absNoisePulse('y',timeSpan,noiseAmp);
    pulseSet{end+1} = absNoisePulse('z',timeSpan,noiseAmp);
    
function pulseSet = addFrictionPulse(pulseSet,timeSpan,absAmp,delta)
if nargin < 3
    absAmp = 1;
end
if nargin < 4
    delta = -10e6;
end
pulseSet{end+1} = absFrictionPulse('x',timeSpan,absAmp,delta);
pulseSet{end+1} = absFrictionPulse('y',timeSpan,absAmp,delta);
pulseSet{end+1} = absFrictionPulse('z',timeSpan,absAmp,delta);

function pulseSet = slapShot(pulseSet,timeSpan,xAmp,yAmp,zAmp,atoms)
pulseSet{end+1} = slap('x',timeSpan,xAmp,atoms);
pulseSet{end+1} = slap('y',timeSpan,yAmp,atoms);
pulseSet{end+1} = slap('z',timeSpan,zAmp,atoms);



function pulse = pulseLaser(time,pulseEnergy,lambda,waist,polarizability,xCenterPulse,yCenterPulse,zCenterPulse)
if nargin < 6
    xCenterPulse = 0;
    yCenterPulse = 0;
    zCenterPulse = 0;
end
pulse.pulseType = sprintf('pulseLaser');
pulse.timeSpan = [(time) (time + 5e-9)];
pulse.lambdaPulse = lambda;
pulse.waistPulse = waist;
intensity = (pulseEnergy / 5e-9) * 2 / (pi * waist * waist);
pulse.EzeroPulse = sqrt(intensity / (3e8 * 8.8e-12)); %3e8 = c, 8.8 e-12 = epsilon0
pulse.xRPulse = (pi * waist * waist) / (lambda);
pulse.polarizabilityPulse = polarizability;
pulse.xCenterPulse = xCenterPulse;
pulse.yCenterPulse = yCenterPulse;
pulse.zCenterPulse = zCenterPulse;

function pulse = chirpPulse(dimension,timeSpan,frequencySpan,absAmp)
pulse.pulseType = sprintf('chirp%c',dimension);
pulse.timeSpan = timeSpan;
pulse.frequencySpan = frequencySpan;
pulse.amp = absAmp;

function pulse = CWLaser(timespan,power,lambda,waist,polarizability,xCenterCW,yCenterCW,zCenterCW)
if nargin < 6
    xCenterCW = 0;
    yCenterCW = 0;
    zCenterCW = 0;
end
pulse.pulseType = sprintf('CWLaser');
pulse.timeSpan = timespan;
pulse.lambdaCW = lambda;
pulse.waistCW = waist;
intensity = power * 2 / (pi * waist * waist);
pulse.EzeroCW = sqrt(intensity / (3e8 * 8.8e-12)); %3e8 = c, 8.8 e-12 = epsilon0
pulse.xRCW = (pi * waist * waist) / (lambda);
pulse.polarizabilityCW = polarizability;
pulse.xCenterCW = xCenterCW;
pulse.yCenterCW = yCenterCW;
pulse.zCenterCW = zCenterCW;

function pulse = slap(dimension,timeSpan,absAmp,atoms,freq)
if nargin < 5
    freq = 0;
end
pulse.pulseType = sprintf('slap%c',dimension);
pulse.timeSpan = timeSpan;
pulse.amp = absAmp;
pulse.atoms = atoms;
pulse.freq = freq;

function pulse = absOpticalPulse(dimension,timeSpan,freq)
pulse.pulseType = sprintf('optical%c',dimension);
pulse.timeSpan = timeSpan;
pulse.phase = 0;
pulse.freq = freq;
pulse.amp = 8e-22* sqrt(.005);

function pulse = absImpulse(dimension,timeSpan,absAmp,freq,phase)
if nargin < 4
    freq = 0;
end
if nargin < 5
    phase = 0;
end
pulse.pulseType = sprintf('E%c',dimension);
pulse.timeSpan = timeSpan;
pulse.phase = phase;
pulse.freq = freq;
pulse.amp = absAmp;

function pulse = absNoisePulse(dimension,timeSpan,absAmp)
pulse.pulseType = sprintf('noise%c',dimension);
pulse.timeSpan = timeSpan;
pulse.phase = 0;
pulse.freq = 0;
pulse.amp = absAmp;

function pulse = absFrictionPulse(dimension,timeSpan,absAmp,delta)
if nargin < 4
   delta = -10e6
end
pulse.pulseType = sprintf('friction%c',dimension);
pulse.timeSpan = timeSpan;
pulse.phase = 0;
pulse.freq = 0;
pulse.amp = absAmp;
pulse.delta = delta;

function pulse = squeezePulse(atom,dimension,phase,timeSpan,squeezeFactor)
pot = potFromAtom(atom,dimension);
pulse.pulseType = sprintf('E2%c',dimension);
pulse.timeSpan = timeSpan;
pulse.phase = phase;
pulse.freq = pot.freq * 2;
numCycles = (timeSpan(2) - timeSpan(1))*pot.freq; 
ampPerCycle = (squeezeFactor^(1/numCycles) - 1);
pulse.amp = 1 * ampPerCycle * pot.forceK;  %work out from squeezeFactor. units of V/m^2




function pot = potFromAtom(atom,dimension)
switch dimension
    case 'x'
        pot = atom.potentialx;
    case 'y'
        pot = atom.potentialy;
    case 'z'
        pot = atom.z;
end

function fieldSet = compileFields(pulseSet,times)
times = doubleTimes(times);
fieldSet = blankSet(length(times));
fieldSet.xAtoms = [];
fieldSet.yAtoms = [];
fieldSet.zAtoms = [];

fieldSet.waistPulse = 1;
fieldSet.EzeroPulse = 0;
fieldSet.xRPulse = 1;
fieldSet.polarizabilityPulse = 0;
fieldSet.OPxCenterPulse = 0;
fieldSet.OPyCenterPulse = 0;
fieldSet.OPzCenterPulse = 0;

fieldSet.waistCW = 1;
fieldSet.EzeroCW = 0;
fieldSet.xRCW = 1;
fieldSet.polarizabilityCW = 0;
fieldSet.OPxCenterCW = 0;
fieldSet.OPyCenterCW = 0;
fieldSet.OPzCenterCW = 0;

fieldSet.chirpx = zeros(1,length(times));
fieldSet.chirpy = zeros(1,length(times));
fieldSet.chirpz = zeros(1,length(times));

fieldSet.delta = 0;

for i = 1:length(pulseSet)
    thisPulse = pulseSet{i};
    fieldSet.(thisPulse.pulseType) = fieldSet.(thisPulse.pulseType) + wiggleFromPulse(thisPulse,times);
    if (strcmp(thisPulse.pulseType,'frictionx')) || (strcmp(thisPulse.pulseType,'frictiony')) || (strcmp(thisPulse.pulseType,'frictionz'))
        fieldSet.delta = thisPulse.delta;
    end
    if strcmp(thisPulse.pulseType,'pulseLaser')
        fieldSet.waistPulse = thisPulse.waistPulse;
        fieldSet.EzeroPulse = thisPulse.EzeroPulse;
        fieldSet.xRPulse = thisPulse.xRPulse;
        fieldSet.polarizabilityPulse = thisPulse.polarizabilityPulse;
        fieldSet.OPxCenterPulse = thisPulse.xCenterPulse;
        fieldSet.OPyCenterPulse = thisPulse.yCenterPulse;
        fieldSet.OPzCenterPulse = thisPulse.zCenterPulse;
        OpticalDepthPulse = fieldSet.polarizabilityPulse * (fieldSet.EzeroPulse * fieldSet.EzeroPulse) / 1.38e-23;
        if OpticalDepthPulse ~= 0
            fprintf('Depth of Pulse Laser Potential: %4.2f K \n',OpticalDepthPulse); %Fix this to handle multiple OPs
        end
    end
    if strcmp(thisPulse.pulseType,'CWLaser')
        fieldSet.waistCW = thisPulse.waistCW;
        fieldSet.EzeroCW = thisPulse.EzeroCW;
        fieldSet.xRCW = thisPulse.xRCW;
        fieldSet.polarizabilityCW = thisPulse.polarizabilityCW;
        fieldSet.OPxCenterCW = thisPulse.xCenterCW;
        fieldSet.OPyCenterCW = thisPulse.yCenterCW;
        fieldSet.OPzCenterCW = thisPulse.zCenterCW;
        OpticalDepthCW = fieldSet.polarizabilityCW * (fieldSet.EzeroCW * fieldSet.EzeroCW) / 1.38e-23;
        if OpticalDepthCW ~= 0
            fprintf('Depth of CW Laser Potential: %4.5f K \n',OpticalDepthCW); %Fix this to handle multiple OPs
        end
    end
    if strcmp(thisPulse.pulseType,'slapx')
        fieldSet.xAtoms = thisPulse.atoms;
    end
    if strcmp(thisPulse.pulseType,'slapy')
        fieldSet.yAtoms = thisPulse.atoms;
    end
    if strcmp(thisPulse.pulseType,'slapz')
        fieldSet.zAtoms = thisPulse.atoms;
    end
end
fieldSet.Ex = fieldSet.Ex + fieldSet.noisex + fieldSet.chirpx;
fieldSet.Ey = fieldSet.Ey + fieldSet.noisey + fieldSet.chirpy;
fieldSet.Ez = fieldSet.Ez + fieldSet.noisez + fieldSet.chirpz;

fieldSet.times = times;

function wiggle = wiggleFromPulse(pulse,times)
switch pulse.pulseType
    case{'noisex','noisey','noisez'}
        dt = times(2) - times(1);
        wiggle = randn(size(times)) * pulse.amp / sqrt(dt);
        
    case{'slapx','slapy','slapz'}
        nArray = pulse.amp * cos(times * pulse.freq * 2 * pi); %ones(size(times));
        wiggle = nArray;
        
    case{'pulseLaser','CWLaser'}
        wiggle = ones(size(times));

    case{'opticalx','opticaly','opticalz'}
        dt = times(2) - times(1);
        t = times(end) - times(1);
        numFlips = round(t * pulse.freq);
        avgFlipInterval = ceil(length(times) / numFlips);
        noiseArray = randn(size(times)) * pulse.amp / sqrt(dt);
        wiggle = noiseArray;
        flipInterval = zeros(1,(numFlips + 1));
        for l = 1:(numFlips + 1)
            flipInterval(l) = ceil(avgFlipInterval + (50 * randn));
        end
        j = 1;
        k = 1;
        for i = 2:length(times)
            if mod(k,flipInterval(j)) ~= 0
                wiggle(i) = wiggle(i-1);
                k = k+1;
            else
                j = j+1;
                k = 1;
            end
        end
        
    case{'chirpx','chirpy','chirpz'}
        intStart = round(pulse.timeSpan(1) / 5e-9);
        if intStart == 0
            intStart = 1;
        end
        intEnd = round(pulse.timeSpan(2) / 5e-9);
        wiggle = zeros(1,length(times));
        chirpTimes = times(intStart:intEnd) - times(intStart);
        wiggle(intStart:intEnd) = pulse.amp * chirp(chirpTimes,pulse.frequencySpan(1),chirpTimes(end),pulse.frequencySpan(2));
        
    otherwise
        wiggle = pulse.amp * cos((times * pulse.freq * 2 * pi) + pulse.phase);
    
        
end
wiggle(times < pulse.timeSpan(1)) = 0;
wiggle(times > pulse.timeSpan(2)) = 0;



function fieldSet = blankSet(n)
fieldSet.Ex = zeros(1,n);
fieldSet.Ey = zeros(1,n);
fieldSet.Ez = zeros(1,n);
fieldSet.noisex = zeros(1,n);
fieldSet.noisey = zeros(1,n);
fieldSet.noisez = zeros(1,n);
fieldSet.opticalx = zeros(1,n);
fieldSet.opticaly = zeros(1,n);
fieldSet.opticalz = zeros(1,n);
fieldSet.E2x = zeros(1,n);
fieldSet.E2y = zeros(1,n);
fieldSet.E2z = zeros(1,n); 
fieldSet.frictionx = zeros(1,n);
fieldSet.frictiony = zeros(1,n);
fieldSet.frictionz = zeros(1,n);
fieldSet.slapx = zeros(1,n);
fieldSet.slapy = zeros(1,n);
fieldSet.slapz = zeros(1,n);
fieldSet.pulseLaser = zeros(1,n);
fieldSet.CWLaser = zeros(1,n);


function moleculeZ = getZPosition(cloud,vector)
if nargin < 2
    vector = cloud.vector;
end
zPositions = zeros(1,cloud.numAtoms);
for i = 1:cloud.numAtoms
    zPositions(i) = vector(3 + (6 * (i-1)));
end
molZPos = zPositions(end);
sortedPositions = sort(zPositions);
for j = 1:cloud.numAtoms
    if sortedPositions(j) == molZPos
        moleculeZ = j;
    end
end

function zVector = getZVector(cloud,vector)
if nargin < 2
    vector = cloud.vector;
end
zPositions = zeros(1,cloud.numAtoms);
for i = 1:cloud.numAtoms
    zPositions(i) = vector(3 + (6 * (i-1)));
end
zVector = sort(zPositions);

function newTimes = doubleTimes(times)
dt = (times(2) - times(1))/ 2;
newLength = length(times) * 2;
newTimes = times(1) + ((0:(newLength-1))*dt);

function plotFinalTemps(tickleFrequencies,finalTemps)
figure('Position',[100   100   500   500],'Name',sprintf('Final Temperatures vs Frequency'));
plot(tickleFrequencies,finalTemps,'b');
xlabel('Tickle Frequency');
ylabel('Final Temperature');

function plotJustTrajectories(cloud,lateTimes,lateTimes2,lateTimes3,lateTimes4,lateTimes5)
doeslateTimes2 = 1;
doeslateTimes3 = 1;
doeslateTimes4 = 1;
doeslateTimes5 = 1;
isMixed = 0;
if strcmp(cloud.isMixed,'mixed')
    isMixed = 1;
end
if strcmp(cloud.isMixed,'mlcls')
    isMixed = 2;
end
if nargin < 2
    lateTimes = [0.5e-3 Inf];
end
if nargin < 3
    doeslateTimes2 = 0;
end
if nargin < 4
    doeslateTimes3 = 0;
end
if nargin < 5
    doeslateTimes4 = 0;
end
if nargin < 6
    doeslateTimes5 = 0;
end
width = (2 + doeslateTimes2 + doeslateTimes3 + doeslateTimes4 + doeslateTimes5) * 300;
figure('Position',[82   107   width   630],'Name',sprintf('%d atom cloud',cloud.numAtoms));
% 
% str1 = " Sr - (" + (cloud.atoms{1}.potentialx.freq / 1000) + " kHz," + (cloud.atoms{1}.potentialy.freq / 1000) + " kHz," + (cloud.atoms{1}.potentialz.freq / 1000)+ " kHz) ";
% if (isMixed == 2)
%     str1 = " ";
% end
% 
% str2 = " ";
% if ((isMixed == 1) || (isMixed == 2))
%     str2 = " SrO - (" + (cloud.atoms{2}.potentialx.freq / 1000) + " kHz," + (cloud.atoms{2}.potentialy.freq / 1000) + " kHz," + (cloud.atoms{2}.potentialz.freq / 1000)+ " kHz) ";
% end
% annotation('textbox', [0.1, 0.99, 0.75, 0.01], 'String', "Atoms: " + cloud.numIons + " Molecules: " + cloud.numMolecules + "  Frequencies: " + str1 + str2,'FitBoxToText','on');

if doeslateTimes5
    subplot(2,6,1);
    plot2dTrajectory(cloud,'x','y');
    subplot(2,6,2);
    plot2dTrajectory(cloud,'x','y',lateTimes);
    subplot(2,6,3);
    plot2dTrajectory(cloud,'x','y',lateTimes2);
    subplot(2,6,4);
    plot2dTrajectory(cloud,'x','y',lateTimes3);
    subplot(2,6,5);
    plot2dTrajectory(cloud,'x','y',lateTimes4);
    subplot(2,6,6);
    plot2dTrajectory(cloud,'x','y',lateTimes5);
    
    subplot(2,6,7);
    plot2dTrajectory(cloud,'z','x');
    subplot(2,6,8);
    plot2dTrajectory(cloud,'z','x',lateTimes);
    subplot(2,6,9);
    plot2dTrajectory(cloud,'z','x',lateTimes2);
    subplot(2,6,10);
    plot2dTrajectory(cloud,'z','x',lateTimes3);
    subplot(2,6,11);
    plot2dTrajectory(cloud,'z','x',lateTimes4);
    subplot(2,6,12);
    plot2dTrajectory(cloud,'z','x',lateTimes5);
else
    if doeslateTimes4
        subplot(2,5,1);
        plot2dTrajectory(cloud,'x','y');
        subplot(2,5,2);
        plot2dTrajectory(cloud,'x','y',lateTimes);
        subplot(2,5,3);
        plot2dTrajectory(cloud,'x','y',lateTimes2);
        subplot(2,5,4);
        plot2dTrajectory(cloud,'x','y',lateTimes3);
        subplot(2,5,5);
        plot2dTrajectory(cloud,'x','y',lateTimes4);
        
        subplot(2,5,6);
        plot2dTrajectory(cloud,'z','x');
        subplot(2,5,7);
        plot2dTrajectory(cloud,'z','x',lateTimes);
        subplot(2,5,8);
        plot2dTrajectory(cloud,'z','x',lateTimes2);
        subplot(2,5,9);
        plot2dTrajectory(cloud,'z','x',lateTimes3);
        subplot(2,5,10);
        plot2dTrajectory(cloud,'z','x',lateTimes4);
    else
        if doeslateTimes3
            subplot(2,4,1);
            plot2dTrajectory(cloud,'x','y');
            subplot(2,4,2);
            plot2dTrajectory(cloud,'x','y',lateTimes);
            subplot(2,4,3);
            plot2dTrajectory(cloud,'x','y',lateTimes2);
            subplot(2,4,4);
            plot2dTrajectory(cloud,'x','y',lateTimes3);
            
            subplot(2,4,5);
            plot2dTrajectory(cloud,'z','x');
            subplot(2,4,6);
            plot2dTrajectory(cloud,'z','x',lateTimes);
            subplot(2,4,7);
            plot2dTrajectory(cloud,'z','x',lateTimes2);
            subplot(2,4,8);
            plot2dTrajectory(cloud,'z','x',lateTimes3);
        else
            if doeslateTimes2    
                subplot(2,3,1);
                plot2dTrajectory(cloud,'x','y');
                subplot(2,3,2);
                plot2dTrajectory(cloud,'x','y',lateTimes);
                subplot(2,3,3);
                plot2dTrajectory(cloud,'x','y',lateTimes2);
                
                subplot(2,3,4);
                plot2dTrajectory(cloud,'z','x');
                subplot(2,3,5);
                plot2dTrajectory(cloud,'z','x',lateTimes);
                subplot(2,3,6);
                plot2dTrajectory(cloud,'z','x',lateTimes2);
            else
                subplot(2,2,1);
                plot2dTrajectory(cloud,'x','y');
                subplot(2,2,2);
                plot2dTrajectory(cloud,'x','y',lateTimes);
                
                subplot(2,2,3);
                plot2dTrajectory(cloud,'z','x');
                subplot(2,2,4);
                plot2dTrajectory(cloud,'z','x',lateTimes);
            end
        end
    end
end


function plotTrajectories(cloud,lateTimes,lateTimes2,lateTimes3,lateTimes4,lateTimes5)
doeslateTimes2 = 1;
doeslateTimes3 = 1;
doeslateTimes4 = 1;
doeslateTimes5 = 1;
isMixed = 0;
if strcmp(cloud.isMixed,'mixed')
    isMixed = 1;
end
if strcmp(cloud.isMixed,'mlcls')
    isMixed = 2;
end
if nargin < 2
    lateTimes = [0.5e-3 Inf];
end
if nargin < 3
    doeslateTimes2 = 0;
end
if nargin < 4
    doeslateTimes3 = 0;
end
if nargin < 5
    doeslateTimes4 = 0;
end
if nargin < 6
    doeslateTimes5 = 0;
end
figure('Position',[82   107   1414   839],'Name',sprintf('%d atom cloud',cloud.numAtoms));
subplot(3,1,1);
plotXYZ(cloud,1)
plotXYZ(cloud,2)
% plotXYZ(cloud,0)
a= dbstack;
title([a(2).name '()']);

subplot(4,1,2);
plotEnergy(cloud);
%subplot(3,1,2);
%plotXYZ(cloud,2)

str1 = " Sr - (" + (cloud.atoms{1}.potentialx.freq / 1000) + " kHz," + (cloud.atoms{1}.potentialy.freq / 1000) + " kHz," + (cloud.atoms{1}.potentialz.freq / 1000)+ " kHz) ";
str2 = " ";
if (isMixed == 2)
    str1 = " ";
    str2 = " SrO - (" + (cloud.atoms{1}.potentialx.freq / 1000) + " kHz," + (cloud.atoms{1}.potentialy.freq / 1000) + " kHz," + (cloud.atoms{1}.potentialz.freq / 1000)+ " kHz) ";
end

if (isMixed == 1)
    str2 = " SrO - (" + (cloud.atoms{2}.potentialx.freq / 1000) + " kHz," + (cloud.atoms{2}.potentialy.freq / 1000) + " kHz," + (cloud.atoms{2}.potentialz.freq / 1000)+ " kHz) ";
end
annotation('textbox', [0.1, 0.99, 0.75, 0.01], 'String', "Atoms: " + cloud.numIons + " Molecules: " + cloud.numMolecules + "  Frequencies: " + str1 + str2,'FitBoxToText','on');

if doeslateTimes2
    subplot(4,6,13);
    plot2dTrajectory(cloud,'x','y');
    subplot(4,6,14);
    plot2dTrajectory(cloud,'x','y',lateTimes);
    subplot(4,6,15);
    plot2dTrajectory(cloud,'x','y',lateTimes2);
else
    subplot(4,4,9);
    plot2dTrajectory(cloud,'x','y');
    subplot(4,4,10);
    plot2dTrajectory(cloud,'x','y',lateTimes);
end
if doeslateTimes3
    subplot(4,6,16);
    plot2dTrajectory(cloud,'x','y',lateTimes3);
end
if doeslateTimes4
    subplot(4,6,17);
    plot2dTrajectory(cloud,'x','y',lateTimes4);
end
if doeslateTimes5
    subplot(4,6,18);
    plot2dTrajectory(cloud,'x','y',lateTimes5);
end

if doeslateTimes2
    subplot(4,6,19);
    plot2dTrajectory(cloud,'z','x');
    subplot(4,6,20);
    plot2dTrajectory(cloud,'z','x',lateTimes);
    subplot(4,6,21);
    plot2dTrajectory(cloud,'z','x',lateTimes2);
else
    subplot(4,4,13);
    plot2dTrajectory(cloud,'z','x');
    subplot(4,4,14);
    plot2dTrajectory(cloud,'z','x',lateTimes);
end
if doeslateTimes3
    subplot(4,6,22);
    plot2dTrajectory(cloud,'z','x',lateTimes3);
end
if doeslateTimes4
    subplot(4,6,23);
    plot2dTrajectory(cloud,'z','x',lateTimes4);
end
if doeslateTimes5
    subplot(4,6,24);
    plot2dTrajectory(cloud,'z','x',lateTimes5);
end


function animated2dTrajectory(cloud,dim1,dim2,timeWindow,filename,narrowRegion,permTrails)
if nargin < 5
    filename = 'animatedTrajectories.gif';
end
if nargin < 6
    narrowRegion = 0;
end
if nargin < 7
    permTrails = 0;
end
fig = figure('Position',[100   100   1414   839],'Name',sprintf('%d atom cloud',cloud.numAtoms));
if nargin < 4
    timeWindow = [0 Inf];
    %titlestring = 'all times';
%else
%    titlestring = sprintf('%3.0f usec < t < %3.0f usec',timeWindow(1)*1e6,timeWindow(2)*1e6);
end
dim1Array = cell(1,cloud.numAtoms);
dim2Array = cell(1,cloud.numAtoms);
for i = 1:cloud.numAtoms
    dim1Array{i} = pullTrajectory(cloud,sprintf('%c%d',dim1,i),timeWindow);
    dim2Array{i} = pullTrajectory(cloud,sprintf('%c%d',dim2,i),timeWindow);
end

cmCoords1 = zeros(length(dim1Array{1}),1);
cmCoords2 = zeros(length(dim2Array{1}),1);
totalMass = 0;
for i = 1:cloud.numIons
    cmCoords1 = cmCoords1 + (88 * dim1Array{i});
    cmCoords2 = cmCoords2 + (88 * dim2Array{i});
    totalMass = totalMass + 88;
end
for j = (cloud.numIons + 1):(cloud.numIons + cloud.numMolecules)
    cmCoords1 = cmCoords1 + (89 * dim1Array{j});
    cmCoords2 = cmCoords2 + (89 * dim2Array{j});
    totalMass = totalMass + 89;
end
cmCoords1 = cmCoords1 / totalMass;
cmCoords2 = cmCoords2 / totalMass;
stepLength = 100;
for i = 1:(floor(length(dim1Array{1}) / stepLength))
    currentTime = round(((i * stepLength) * 1.697688e-8) * 1e5);
    currentTime = currentTime / 100;
    for j = 1:cloud.numAtoms
        if permTrails
            x = 1;
        else
            if i < 2
                x = 1;
            else
                x = (i - 1) * stepLength;
            end
        end
    end
    for k = 1:cloud.numIons
        plot(dim1Array{k}(x:(i * stepLength)),dim2Array{k}(x:(i * stepLength)),dim1Array{k}((i * stepLength)),dim2Array{k}((i * stepLength)),'.b','MarkerSize',25)
        hold on
    end
    for k = (cloud.numIons + 1):(cloud.numIons + cloud.numMolecules)
        plot(dim1Array{k}(x:(i * stepLength)),dim2Array{k}(x:(i * stepLength)),'-r',dim1Array{k}((i * stepLength)),dim2Array{k}((i * stepLength)),'.r','MarkerSize',25)
        hold on
    end
    plot(cmCoords1(i * stepLength),cmCoords2(i * stepLength),'+','MarkerSize',15,'MarkerEdgeColor','k')
    hold on
    axis([-10e-4 10e-4 -10e-4 10e-4]) 
    forceSquare();
    xlabel(dim1);
    ylabel(dim2);
    if narrowRegion ~= 0
        line([narrowRegion(1) narrowRegion(1)],[-1 1],'Color','red');
        line([narrowRegion(2) narrowRegion(2)],[-1 1],'Color','red');
    end
    trueCurrentTime = (timeWindow(1)*1e3) + currentTime;
    title("t = " + trueCurrentTime + " milliseconds");
    drawnow
    frame = getframe(fig);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1
        imwrite(imind,cm,filename,'gif','DelayTime',0.05,'Loopcount',inf); 
    else 
        imwrite(imind,cm,filename,'gif','DelayTime',0.05,'WriteMode','append'); 
    end 
    hold off
end



function plot2dTrajectory(cloud,dim1,dim2,timeWindow)
if nargin < 4
    timeWindow = [0 Inf];
    titlestring = 'all times';
else
    titlestring = sprintf('%3.0f usec < t < %3.0f usec',timeWindow(1)*1e6,timeWindow(2)*1e6);
end
for i = 1:cloud.numAtoms
    x1 = pullTrajectory(cloud,sprintf('%c%d',dim1,i),timeWindow);
    y1 = pullTrajectory(cloud,sprintf('%c%d',dim2,i),timeWindow);
    plot(x1,y1);
    hold all;
end
forceSquare();
xlabel(dim1);
ylabel(dim2);
title(titlestring);x = randn;
        


function forceSquare()
a = xlim;
b = ylim;
extent = max(abs([a b]));
xlim([-extent extent]);
ylim([-extent extent]);

function x = cauchy_dist(location_parameter, scale_parameter)
p_cdf = rand(); %uniform random from 0->1, since cdf by definition 0->1
x = location_parameter + scale_parameter*tan(pi*(p_cdf-0.5)); %solve cdf eqn for x

function plotFluorescence(cloud,timeSpan,region)
if nargin < 2
    timeSpan = [0 Inf];
end
if nargin < 3
    region = [-1e-3 1e-3];
end
figure('Position',[600   200   1000   600]);
c = makeConstants;
indivEnergyChanges = zeros(cloud.numIons,length(cloud.individualEnergies));
for i = 1:(length(cloud.individualEnergies)-1)
    indivEnergyChanges(:,i) = cloud.individualEnergies(1:cloud.numIons,i) - cloud.individualEnergies(1:cloud.numIons,i+1);
end
indivEnergyChanges(:,end) = indivEnergyChanges(:,end-1);
for i = 1:length(indivEnergyChanges)
    for j = 1:cloud.numIons
       if indivEnergyChanges(j,i) < 0
           indivEnergyChanges(j,i) = 0;
       end
    end
end
energyTrajectories = zeros(2*cloud.numIons,length(indivEnergyChanges));
energyTrajectories(1:cloud.numIons,:) = indivEnergyChanges;
for i = 1:cloud.numIons
    for j = 1:length(indivEnergyChanges)
        energyTrajectories(cloud.numIons + i,j) = cloud.trajectories(j,6 * (i-1) + 3);
    end
end
actualPhotons = zeros(1,length(indivEnergyChanges));
currentEnergy = 0;
for i = 1:length(indivEnergyChanges)
    for j = 1:cloud.numIons
         if region(1) < energyTrajectories(cloud.numIons + j,i) && region(2) > energyTrajectories(cloud.numIons + j,i)
             currentEnergy = currentEnergy + energyTrajectories(j,i);
         end
    end
    while currentEnergy > 0
        randEnergyChange = -1;
        while randEnergyChange < 0
            randAtomFreq = cauchy_dist(7.1116297286e14,2e7);
%             randLaserFreq = normrnd((7.1116297286e14 - 1e7),2.34e5);
            randEnergyChange = c.H * (randAtomFreq - (7.1116297286e14 - 1e7));
        end
        if currentEnergy > randEnergyChange
            currentEnergy = currentEnergy - randEnergyChange;
            actualPhotons(i) = actualPhotons(i) + 1;
        end
        if currentEnergy < randEnergyChange && currentEnergy > 0
            random = rand;
            if random <= (currentEnergy/randEnergyChange)
                actualPhotons(i) = actualPhotons(i) + 1;
            end
            currentEnergy = currentEnergy - randEnergyChange;
        end
    end
    currentEnergy = 0;
end
numBins = ceil(length(actualPhotons)/300);
fluorescence = zeros(1,numBins);
for i = 1:(numBins - 1)
    for j = 1:300
        fluorescence(i) = fluorescence(i) + actualPhotons(300*(i-1) + j);
    end
    fluorescence(i) = fluorescence(i) / 300;
end
for i = ((numBins-1)*300):length(actualPhotons)
    fluorescence(numBins) = fluorescence(numBins) + actualPhotons(i);
end
fluorescence(numBins) = fluorescence(numBins) / (length(actualPhotons) - ((numBins-1) * 300));
fluorescence = (fluorescence / (1.697688e-8));
fluorescence = round(0.003 * fluorescence);
newTimes = (1:numBins) * 300 * 1.697688e-8;
indexSpan = round(timeSpan / (300 * 1.697688e-8));
if indexSpan(1) == 0
    indexSpan(1) = 1;
end
if indexSpan(2) == Inf
    indexSpan(2) = numBins;
end
plot(newTimes(indexSpan(1):indexSpan(2)),fluorescence(indexSpan(1):indexSpan(2)));
xlabel('time');
ylabel('photons per second');
title(sprintf('Fluorescence from the region %g m < z < %g m',region(1),region(2)));

function fluor = idealFluorescence(finalCloud,timeSpan)
if nargin < 2
    timeSpan = [0 Inf];
end
fieldSet = finalCloud.fieldSet;
const = makeConstants;

l1 = length(pullTrajectory(finalCloud,'vx1'));
l2 = length(fieldSet.frictionx);

fluorescence = zeros(1,min(l1,l2));


delta = fieldSet.delta; %  detuning: + => blue, - => red   CHANGE THIS ONE TO ALTER LASER
        
omega0 = 7.1116297286e14; % angular frequency of Sr transition
omega = omega0 + delta; % angular frequency of cooling laser
k = omega / const.C; % wavenumber of cooling laser
gamma = 2 * pi * 21e6;

for i = 1:finalCloud.numIons
    VX = pullTrajectory(finalCloud,sprintf('vx%d',i));
    VY = pullTrajectory(finalCloud,sprintf('vy%d',i));
    VZ = pullTrajectory(finalCloud,sprintf('vz%d',i));
    
    if length(fieldSet.frictionx) > length(VX)
        frictionx = fieldSet.frictionx(1:length(VX));
        frictiony = fieldSet.frictiony(1:length(VY));
        frictionz = fieldSet.frictionz(1:length(VZ));
        vx = VX;
        vy = VY;
        vz = VZ;
    elseif length(fieldSet.frictionx) < length(VX)
        vx = VX(1:length(fieldSet.frictionx));
        vy = VY(1:length(fieldSet.frictiony));
        vz = VZ(1:length(fieldSet.frictionz));
        frictionx = fieldSet.frictionx;
        frictiony = fieldSet.frictiony;
        frictionz = fieldSet.frictionz;
    else
        vx = VX;
        vy = VY;
        vz = VZ;
        frictionx = fieldSet.frictionx;
        frictiony = fieldSet.frictiony;
        frictionz = fieldSet.frictionz;
    end
    Rvx = zeros(1,length(vx));
    Rvy = zeros(1,length(vy));
    Rvz = zeros(1,length(vz));
    for j = 1:length(vx)
        Rvx(j) = ((gamma / 2) * frictionx(j)) / (1 + frictionx(j) + ((4 / (gamma^2)) * (delta - (k * vx(j) ) ).^2 ) );
        Rvy(j) = ((gamma / 2) * frictiony(j)) / (1 + frictiony(j) + ((4 / (gamma^2)) * (delta - (k * vy(j) ) ).^2 ) );
        Rvz(j) = ((gamma / 2) * frictionz(j)) / (1 + frictionz(j) + ((4 / (gamma^2)) * (delta - (k * vz(j) ) ).^2 ) );
    end
    fluorescence = fluorescence + Rvx + Rvy + Rvz;
end

fluorescence = round(0.003 * fluorescence);
indexSpan = round(timeSpan / (1.697688e-8));
if indexSpan(1) == 0
    indexSpan(1) = 1;
end
if indexSpan(2) == inf
    indexSpan(2) = length(fluorescence);
end
fluor = fluorescence(indexSpan(1):indexSpan(2));

function plotIndividualEnergies(cloud)
figure('Position',[600   200   600   400]);
times = cloud.times;
individualTemperatures = cloud.individualEnergies/1.3e-23;
for i = 1:cloud.numAtoms
    plot(times,individualTemperatures(i,:));
    hold on
end

xlabel('time');
ylabel('Temperature of Each Ion');

function plotEnergy(cloud,lineType)
if nargin < 2
    lineType = 'b';
end
times = cloud.times;
% plot(times,cloud.fracErrors,'m');+6
% energies = max(cloud.fracErrors) * cloud.energies / max(cloud.energies);
% hold all;
% plot(times,energies,'b');

temperatures = cloud.energies/1.3e-23;
errorString = sprintf('minimum temperature %3.3f K\nmaximum integration error %3.2e',min(temperatures),max(cloud.fracErrors));
plot(times,temperatures,lineType);        atom.friction = 1;

xlabel('time');
ylabel('Temperature');
if nargin == 1
    addText(errorString);
end

function plotTrueTemp(cloud,lineType)
figure('Position',[600   200   600   400]);
if nargin < 2
    lineType = 'b';
end
times = cloud.times;
temperatures = (cloud.energies)/(cloud.numAtoms * 1.3e-23);
plot(times,temperatures,lineType);
xlabel('time');
ylabel('Temperature');


function p = plotJustZ(cloud,n,whichChannel)
if n > cloud.numAtoms
    return
end
if nargin < 3
    whichChannel = 'z';
end
times = cloud.times;
if n == 0
    zArray = cell(1,cloud.numAtoms);
    for i = 1:cloud.numAtoms
        zArray{i} = pullTrajectory(cloud,sprintf('%c%d','z',i));
    end
    cmCoordsz = zeros(length(times),1);
    totalMass = 0;
    for i = 1:cloud.numIons
        cmCoordsz = cmCoordsz + (88 * zArray{i});
        totalMass = totalMass + 88;
    end
    for i = (cloud.numIons + 1):(cloud.numIons + cloud.numMolecules)
        cmCoordsz = cmCoordsz + (89 * zArray{i});
        totalMass = totalMass + 89;
    end
    cmCoordsz = cmCoordsz / totalMass;
    hold all;
    plot(times,cmCoordsz,'k');
    legend('COMz');
else
    times = cloud.times;
    labelz = sprintf('%c%d',whichChannel,n);
    z = pullTrajectory(cloud,labelz);
    p = plot(times,z,'k:');
    a = legend;
    oldLegend = a.String;
    oldLegend{end} = labelz;
    legend(oldLegend);
end

function plotXYZ(cloud,n)
if n > cloud.numAtoms
    return
end
times = cloud.times;
if n == 0
    xArray = cell(1,cloud.numAtoms);
    yArray = cell(1,cloud.numAtoms);
    zArray = cell(1,cloud.numAtoms);
    for i = 1:cloud.numAtoms
        xArray{i} = pullTrajectory(cloud,sprintf('%c%d','x',i));
        yArray{i} = pullTrajectory(cloud,sprintf('%c%d','y',i));
        zArray{i} = pullTrajectory(cloud,sprintf('%c%d','z',i));
    end
    cmCoordsx = zeros(length(times),1);
    cmCoordsy = zeros(length(times),1);
    cmCoordsz = zeros(length(times),1);
    totalMass = 0;
    for i = 1:cloud.numIons
        cmCoordsx = cmCoordsx + (88 * xArray{i});
        cmCoordsy = cmCoordsy + (88 * yArray{i});
        cmCoordsz = cmCoordsz + (88 * zArray{i});
        totalMass = totalMass + 88;
    end
    for i = (cloud.numIons + 1):(cloud.numIons + cloud.numMolecules)
        cmCoordsx = cmCoordsx + (89 * xArray{i});
        cmCoordsy = cmCoordsy + (89 * yArray{i});
        cmCoordsz = cmCoordsz + (89 * zArray{i});
        totalMass = totalMass + 89;
    end
    cmCoordsx = cmCoordsx / totalMass;
    cmCoordsy = cmCoordsy / totalMass;
    cmCoordsz = cmCoordsz / totalMass;
    hold all;
    plot(times,cmCoordsx,'r',times,cmCoordsy,'b',times,cmCoordsz,'k');
    legend('COMx','COMy','COMz');
else
    times = cloud.times;
    labelx = sprintf('x%d',n);
    labely = sprintf('y%d',n);
    labelz = sprintf('z%d',n);
    x = pullTrajectory(cloud,labelx);
    y = pullTrajectory(cloud,labely);
    z = pullTrajectory(cloud,labelz);
    plot(times,x,'r');
    hold all;
    plot(times,y,'b');
    plot(times,z,'k');
    legend(labelx,labely,labelz);
end

function coord = pullTrajectory(cloud,whichVar,timeWindow)
if nargin < 3
    timeWindow = [0 Inf];
end
firsti = find(cloud.times >= timeWindow(1),1,'first');
lasti = find(cloud.times <= timeWindow(2),1,'last');
for i = 1:length(cloud.varOrder)
    if strcmp(whichVar,cloud.varOrder{i})
        coord = cloud.trajectories(:,i);
    end
end
coord = coord(firsti:lasti);

function plotFields(cloud)
%whichToPlot = whichFields(fieldSet);
%numToPlot = length(whichToPlot);
fieldSet = cloud.fieldSet;
zPotentialTime = cloud.zPotentialTime;

figure('Position',[200   100   600   500]);
t = fieldSet.times;

subplot((zPotentialTime.zChanges + 4),1,1);
plot(t,fieldSet.frictionx * 1e20);
hold all;
plot(t,fieldSet.frictiony * 1e20);
plot(t,fieldSet.frictionz * 1e20);
legend('Friction x','Friction y','Friction z');

subplot((zPotentialTime.zChanges + 4),1,2);
plot(t,fieldSet.Ex * 1e20);
hold all;
plot(t,fieldSet.Ey * 1e20);
plot(t,fieldSet.Ez * 1e20);
legend('Ex','Ey','Ez');

subplot((zPotentialTime.zChanges + 4),1,3);
plot(t,fieldSet.opticalx * 1e20);
hold all;
plot(t,fieldSet.opticaly * 1e20);
plot(t,fieldSet.opticalz * 1e20);
legend('opticalx','opticaly','opticalz');

subplot((zPotentialTime.zChanges + 4),1,4);
plot(t,fieldSet.E2x * 1e20);
hold all;
plot(t,fieldSet.E2y * 1e20);
plot(t,fieldSet.E2z * 1e20);
legend('squeezex','squeezey','squeezez');

if zPotentialTime.zChanges
    subplot(5,1,5);
    plot(zPotentialTime.zTimes,zPotentialTime.zPot);
    hold all;
    legend('Z Potential');
end

function which = whichFields(fieldSet)
which = {};
fields = fieldnames(S);
for i = 1:length(fields)
    if strcmp(fields{i},'times') == 0
        thisChannel = fieldSet.(fields{i});
        if max(abs(thisChannel)) > 0
            which{end+1} = fields{i};
        end
    end
end

function finalCloud = evolveCloud(ionCloud,pulseSet,t,tStartChange,newPot1,tSecondChange,newPot2)
if nargin < 4
    doesChange = 0;
    doesChangeTwice = 0;
end
if nargin == 5 
    doesChange = 1;
    doesChangeTwice = 0;
end
if nargin > 5
    doesChange = 1;
    doesChangeTwice = 1;
end
%no dynamic timing yet
%variable names h,k1,k2,k3,k4 picked to align with  https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
changeTime = 1.0e-4;
changeIndex = round(changeTime / 1.697688e-8);
h = ionCloud.dt;
numTimes = round(t / h);
fieldSet = compileFields(pulseSet,0:h:t);
zPotentialTime.zPot = [];
zPotentialTime.zTimes = [];
switch doesChange
    case 1
        zPotentialTime.zChanges = 1;
        iStartChange = round(tStartChange / 1.697688e-8);
    case 0
        zPotentialTime.zChanges = 0;
end
if doesChangeTwice
    iSecondChange = round(tSecondChange / 1.697688e-8);
else
    iSecondChange = 0;
end
trajectories = zeros(numTimes,length(ionCloud.vector));
t = fieldSet.times(1);
y = ionCloud.vector;

finalCloud = ionCloud;
energies = zeros(1,numTimes);
individualEnergies = zeros(ionCloud.numAtoms,numTimes);
fracErrors = energies;
simTimes = energies;
energies(1) = energy(y,ionCloud);
initEnergies = individualEnergy(y,ionCloud);
for i = 1:ionCloud.numAtoms
    individualEnergies(i,1) = initEnergies(i);
end
simTimes(1) = t;
trajectories(1,:) = y;

h = h * 1.697688;
i = 2;
if ionCloud.atoms{1}.friction == 1
    originalPot = ionCloud.potentialSet{3}.freq;
else
    originalPot = (ionCloud.potentialSet{3}.freq / 0.912);
end
if doesChange
    atomPotStep = ((newPot1 - originalPot) / changeIndex);
    moleculePotStep = (((newPot1*0.912) - (originalPot*0.912)) / changeIndex);
end
while t < fieldSet.times(end)
    if doesChange
        if i == iStartChange
            fprintf('About to start changing trap\n');
        end
        if (iStartChange < i) && (i < iStartChange + changeIndex + 1)
            k = i - iStartChange;
            for j = 1:ionCloud.numAtoms
                if ionCloud.atoms{j}.friction == 1
                    newPot = defaultPotentialZ(originalPot + (atomPotStep * k),ionCloud.atoms{j}.mass);
                else
                    newPot = defaultPotentialZ((originalPot * 0.912) + (moleculePotStep * k),ionCloud.atoms{j}.mass);  
                end
                ionCloud.potentialSet{(j*3)} = newPot;
            end
        end
        if i == iStartChange + changeIndex + 1
            fprintf('Finished changing trap\n');
        end
        
        if (i == iSecondChange) && doesChangeTwice
            fprintf('About to start changing trap\n');
            if ionCloud.atoms{1}.friction == 1
                updatedPot = ionCloud.potentialSet{3}.freq;
            else
                updatedPot = (ionCloud.potentialSet{3}.freq / 0.912);
            end
        end
        if (iSecondChange < i) && (i < iSecondChange + changeIndex + 1) && doesChangeTwice
            k = i - iSecondChange;
            atomPotStep2 = ((newPot2 - updatedPot) / changeIndex);
            moleculePotStep2 = (((newPot2*0.912) - (updatedPot*0.912)) / changeIndex);
            for j = 1:ionCloud.numAtoms
                if ionCloud.atoms{j}.friction == 1
                    newPot = defaultPotentialZ(updatedPot + (atomPotStep2 * k),ionCloud.atoms{j}.mass);
                else
                    newPot = defaultPotentialZ((updatedPot * 0.912) + (moleculePotStep2 * k),ionCloud.atoms{j}.mass);  
                end
                ionCloud.potentialSet{(j*3)} = newPot;
            end
        end
        if (i == iSecondChange + changeIndex + 1) && doesChangeTwice
            fprintf('Finished changing trap\n');
            if ionCloud.atoms{1}.friction == 1
                updatedPot = ionCloud.potentialSet{3}.freq;
            else
                updatedPot = (ionCloud.potentialSet{3}.freq / 0.912);
            end
        end
    end
    currentPot = ionCloud.potentialSet{3}.freq;
    [fracError] = rkError(y,t,h,ionCloud,fieldSet);
    [newy,newt,collided] = rkStep(y,t,h,ionCloud,fieldSet,0);  
    y = newy;
    t = newt;
    for j = 1:ionCloud.numAtoms
        finalCloud.numCollisions = finalCloud.numCollisions + collided(j);
    end
    if i > length(simTimes)
        trajectories = doubleLength(trajectories);
        energies = doubleLength(energies);
        individualEnergies = doubleLength(individualEnergies);
        simTimes = doubleLength(simTimes);
        fracErrors = doubleLength(fracErrors);
    end
    if mod(i,5000) == 0
        fprintf('time %f/%f, %d steps\n',t,max(fieldSet.times),i);
    end
    
    trajectories(i,:) = y;
    energies(i) = energy(y,ionCloud);
    currentEnergies = individualEnergy(y,ionCloud);
    for j = 1:ionCloud.numAtoms
        individualEnergies(j,i) = currentEnergies(j);
    end
    simTimes(i) = t;
    fracErrors(i) = fracError;
    zPotentialTime.zPot(i-1) = currentPot;
    for j = 1:ionCloud.numAtoms
        zForceK = ionCloud.potentialSet{((j-1)*3) + 3}.forceK;
        Vx = ionCloud.potentialSet{((j-1)*3) + 1}.voltage;
        Vy = ionCloud.potentialSet{((j-1)*3) + 2}.voltage;
        ionCloud.potentialSet{((j-1)*3) + 1} = defaultPotential(Vx,t,ionCloud.atoms{j}.mass,'x',zForceK);
        ionCloud.potentialSet{((j-1)*3) + 2} = defaultPotential(Vy,t,ionCloud.atoms{j}.mass,'y',zForceK);
    end
    i = i+1;
end
for i = 1:(length(zPotentialTime.zPot))
    zPotentialTime.zTimes(end + 1) = (i / length(zPotentialTime.zPot)) * t;
end

trajectories = trajectories(1:i-1,:);  %trim down to actual length
energies = energies(1:i-1);  %trim down to actual length
individualEnergies = individualEnergies(:,1:i-1);
fracErrors = fracErrors(1:i-1);
simTimes = simTimes(1:i-1);   %these horrible 8 lines preallocate arrays a bit..
1;

finalCloud.vector = y;
finalCloud.fracErrors = fracErrors;
finalCloud.energies = energies;
finalCloud.individualEnergies = individualEnergies;
finalCloud.times = simTimes;
finalCloud.trajectories = trajectories;
finalCloud.fieldSet = fieldSet;
finalCloud = updateAtoms(finalCloud);
finalCloud.zPotentialTime = zPotentialTime;

function [fracError] = rkError(y,t,h,ionCloud,fieldSet)
%returns fractional energy error 
energyStart= energy(y,ionCloud);
[newy,newt,collided] = rkStep(y,t,h,ionCloud,fieldSet,1);
energyEnd = energy(newy,ionCloud);
delE = abs(energyStart-energyEnd);
fracError = delE / energyStart;

function [newy,newt,collided] = rkStep(y,t,h,ionCloud,fieldSet,conservative)
collided = zeros(1,ionCloud.numAtoms);
if ionCloud.bufferGasBool
    for i = 1:ionCloud.numAtoms
        vecPos = (i-1) * 6;
        const = makeConstants();
        mHe = 4.0026 * const.AMU;
        tempHe = ionCloud.bufferGasTemp; %Kelvin
        d = 580e-12; %meters
        atomSpeedSquared = (y(vecPos+4))^2 + (y(vecPos+5))^2 + (y(vecPos+6))^2;
        meanFreePathTime = 1 / (((atomSpeedSquared + (3 * const.BOLTZMANN * tempHe / mHe))^(1/2)) * pi * d * d * ionCloud.bufferGasDensity * 1e6);
        
        if (t > ionCloud.bufferGasTime(1)) && (t < ionCloud.bufferGasTime(2))
            p = 1.697688e-8 / meanFreePathTime;
            n = rand;
            if n <= p                
                collided(i) = 1;
            end
        end
    end
end
    
    k1 = h * ydot(t,y,ionCloud,fieldSet,conservative,collided);  %use only k1 for zeros order R-K
    k2 = h * ydot(t+h/2,y+k1/2,ionCloud,fieldSet,conservative,collided);
    k3 = h * ydot(t+h/2,y+k2/2,ionCloud,fieldSet,conservative,collided);
    k4 = h * ydot(t+h,y+k3,ionCloud,fieldSet,conservative,collided);    
    newy = y + ((k1 + (2*k2) + (2*k3) + k4) /6); %this line is RK4 
    newt = t + h;
    
    
function fourier(cloud,dim)   
x = pullTrajectory(cloud,dim);
times = cloud.times;
L = length(x);
if mod(L,2) == 1
    L = L - 1;
    x = x(1:L);
    times = times(1:L);
end
y = fft(x);
Y = abs(y);
p = Y(1:L/2 +1);
p(2:end-1) = 2*p(2:end-1);
deltaT = times(2)-times(1);
deltaF = 1/deltaT;
f = deltaF*(0:(L/2))/L;
figure;
plot(f,p);
xlim([0,5e6]);

function v = doubleLength(v)
sz = size(v);
newlength = length(v) * 2;
if sz(1) > sz(2)
    v(newlength,1) = 0;
else
    v(1,newlength) = 0;
end

    
function i = findIndex(times,t)
dt = times(2) - times(1);
i = (t - times(1))/dt;
i = round(i+1);
i = min(i,length(times));

function ydotVals = ydot(t,y,cloud,fieldSet,conservative,collided) %see https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
%change this to reflect a potential? or simply hardwire (F = k_e q1 q2) / r^2) but ok for now
ydotVals = y * 0;
potentialSet = cloud.potentialSet;
ii = findIndex(fieldSet.times,t);
const = makeConstants();
for i = 1:cloud.numAtoms  %each set is [x y z vx vy vz]
    mass = cloud.atoms{i}.mass;
    vecPos = (i-1) * 6;
    potPos = (i-1) * 3;
    forcex = -y(vecPos+1) * potentialSet{potPos+1}.trueForceK;  %trap potential
    forcey = -y(vecPos+2) * potentialSet{potPos+2}.trueForceK; 
    forcez = -y(vecPos+3) * potentialSet{potPos+3}.forceK;
    
    spotSizePulse = fieldSet.waistPulse * (sqrt(1 + ((y(vecPos+1) - fieldSet.OPxCenterPulse) / fieldSet.xRPulse)^2 ));
    EfieldPulse = fieldSet.EzeroPulse * (fieldSet.waistPulse / spotSizePulse) * exp(-( ((y(vecPos+2) - fieldSet.OPyCenterPulse)^2) + ((y(vecPos+3) - fieldSet.OPzCenterPulse)^2) ) / (spotSizePulse * spotSizePulse));
    twoAlphaEPulse = 2 * fieldSet.polarizabilityPulse * EfieldPulse;
    
    forcex = forcex - (fieldSet.pulseLaser(ii) * twoAlphaEPulse * 2 * ((y(vecPos+1) - fieldSet.OPxCenterPulse) / (fieldSet.xRPulse^2)) * (spotSizePulse^(-2)) * EfieldPulse * ( ( ((y(vecPos+2) - fieldSet.OPyCenterPulse)^2) + ((y(vecPos+3) - fieldSet.OPzCenterPulse)^2) ) - (0.5 * (fieldSet.waistPulse^2))));
    forcey = forcey - (fieldSet.pulseLaser(ii) * twoAlphaEPulse * 2 * (y(vecPos+2) - fieldSet.OPyCenterPulse) * (spotSizePulse^(-2)) * EfieldPulse);
    forcez = forcez - (fieldSet.pulseLaser(ii) * twoAlphaEPulse * 2 * (y(vecPos+3) - fieldSet.OPzCenterPulse) * (spotSizePulse^(-2)) * EfieldPulse);
    
    spotSizeCW = fieldSet.waistCW * (sqrt(1 + ((y(vecPos+1) - fieldSet.OPxCenterCW) / fieldSet.xRCW)^2 ));
    EfieldCW = fieldSet.EzeroCW * (fieldSet.waistCW / spotSizeCW) * exp(-( ((y(vecPos+2) - fieldSet.OPyCenterCW)^2) + ((y(vecPos+3) - fieldSet.OPzCenterCW)^2) ) / (spotSizeCW * spotSizeCW));
    twoAlphaECW = 2 * fieldSet.polarizabilityCW * EfieldCW;
    
    forcex = forcex - (fieldSet.CWLaser(ii) * twoAlphaECW * 2 * ((y(vecPos+1) - fieldSet.OPxCenterCW) / (fieldSet.xRCW^2)) * (spotSizeCW^(-2)) * EfieldCW * ( ( ((y(vecPos+2) - fieldSet.OPyCenterCW)^2) + ((y(vecPos+3) - fieldSet.OPzCenterCW)^2) ) - (0.5 * (fieldSet.waistCW^2))));
    forcey = forcey - (fieldSet.CWLaser(ii) * twoAlphaECW * 2 * (y(vecPos+2) - fieldSet.OPyCenterCW) * (spotSizeCW^(-2)) * EfieldCW);
    forcez = forcez - (fieldSet.CWLaser(ii) * twoAlphaECW * 2 * (y(vecPos+3) - fieldSet.OPzCenterCW) * (spotSizeCW^(-2)) * EfieldCW);
    
    
    isSlappedX = 0;
    isSlappedY = 0;
    isSlappedZ = 0;
    for j = 1:length(fieldSet.xAtoms)
        if i == fieldSet.xAtoms(j)
            isSlappedX = 1;
        end
    end
    for j = 1:length(fieldSet.yAtoms)
        if i == fieldSet.yAtoms(j)
            isSlappedY = 1;
        end
    end
    for j = 1:length(fieldSet.zAtoms)
        if i == fieldSet.zAtoms(j)
            isSlappedZ = 1;
        end
    end
    
    if conservative == 0  %add in drives
        forcex = forcex + (fieldSet.Ex(ii));  %drive field and noise
        forcex = forcex - (y(vecPos+1) * fieldSet.E2x(ii)); %squeezing or trap scaling
        if isSlappedX
            forcex = forcex + fieldSet.slapx(ii);
        end
        
        forcey = forcey + (fieldSet.Ey(ii));
        forcey = forcey - (y(vecPos+2) * fieldSet.E2y(ii)); %squeezing or trap scaling
        if isSlappedY
            forcey = forcey + fieldSet.slapy(ii);
        end
            
        forcez = forcez + (fieldSet.Ez(ii));
        forcez = forcez - (y(vecPos+3) * fieldSet.E2z(ii)); %squeezing or trap scaling
        if isSlappedZ
            forcez = forcez + fieldSet.slapz(ii);
        end
        
    end
    if cloud.atoms{i}.friction == 1
        
        delta = fieldSet.delta; %  detuning: + => blue, - => red   CHANGE THIS ONE TO ALTER LASER
        
        omega0 = 7.1116297286e14; % angular frequency of Sr transition
        omega = omega0 + delta; % angular frequency of cooling laser
        k = omega / const.C; % wavenumber of cooling laser
        gamma = 2 * pi * 21e6; % linewidth
        
        % fieldSet.frictionx(ii) is the fraction I / I-sat
        % I think I need to make the frequency of something a spread
        % instead of a delta function
        
        Rvx = ((gamma / 2) * fieldSet.frictionx(ii)) / (1 + fieldSet.frictionx(ii) + ((4 / (gamma^2)) * (delta - (k * y(vecPos+4) ) )^2 ) );
        Rvy = ((gamma / 2) * fieldSet.frictiony(ii)) / (1 + fieldSet.frictiony(ii) + ((4 / (gamma^2)) * (delta - (k * y(vecPos+5) ) )^2 ) );
        Rvz = ((gamma / 2) * fieldSet.frictionz(ii)) / (1 + fieldSet.frictionz(ii) + ((4 / (gamma^2)) * (delta - (k * y(vecPos+6) ) )^2 ) );
        
        theta = pi * rand();
        phi = 2 * pi * rand();
        
        forcex = forcex + (const.HBAR * k * Rvx) + (const.HBAR * k * sqrt(Rvx) * sin(theta) * cos(phi));
        forcey = forcey + (const.HBAR * k * Rvy) + (const.HBAR * k * sqrt(Rvy) * sin(theta) * sin(phi));
        forcez = forcez + (const.HBAR * k * Rvz) + (const.HBAR * k * sqrt(Rvz) * cos(theta));
        
        
    end
    if cloud.atoms{i}.friction == 0
        forcex = forcex + (fieldSet.opticalx(ii));
        forcey = forcey + (fieldSet.opticaly(ii));
        forcez = forcez + (fieldSet.opticalz(ii));
    end
    
    if collided(i)
        const = makeConstants();
        mHe = 4.0026 * const.AMU;
        tempHe = cloud.bufferGasTemp; %Kelvin
        mAtom = cloud.atoms{i}.mass;
        
        vHeX = normrnd(-y(vecPos+4),((2 * const.BOLTZMANN * tempHe) / mHe)^(1/2));
        vHeY = normrnd(-y(vecPos+5),((2 * const.BOLTZMANN * tempHe) / mHe)^(1/2));
        vHeZ = normrnd(-y(vecPos+6),((2 * const.BOLTZMANN * tempHe) / mHe)^(1/2));
                
        vAtomNewX = (((2*mHe)/(mHe + mAtom))*vHeX) - (((mHe-mAtom)/(mHe+mAtom))*y(vecPos+4));
        vAtomNewY = (((2*mHe)/(mHe + mAtom))*vHeY) - (((mHe-mAtom)/(mHe+mAtom))*y(vecPos+5));
        vAtomNewZ = (((2*mHe)/(mHe + mAtom))*vHeZ) - (((mHe-mAtom)/(mHe+mAtom))*y(vecPos+6));
                
        deltaPx = mAtom*(vAtomNewX - y(vecPos+4));
        deltaPy = mAtom*(vAtomNewY - y(vecPos+5));
        deltaPz = mAtom*(vAtomNewZ - y(vecPos+6));
        deltaT = 1.697688e-8;
                
        forcex = forcex + (deltaPx / deltaT);
        forcey = forcey + (deltaPy / deltaT);
        forcez = forcez + (deltaPz / deltaT);
    end
    
    ydotVals(vecPos+1:vecPos+3) = y(vecPos+4:vecPos+6);
    ydotVals(vecPos+4) = forcex/mass; %dv_x/dt = f_x/m
    ydotVals(vecPos+5) = forcey/mass; %dv_x/dt = f_x/m
    ydotVals(vecPos+6) = forcez/mass; %dv_x/dt = f_x/m
end
%now put in interactions
if cloud.interacting
    for i = 1:cloud.numAtoms
        for j = i+1:cloud.numAtoms
            ivecPos = (i-1) * 6;
            jvecPos = (j-1) * 6;
            iPos = y(ivecPos+(1:3));
            jPos = y(jvecPos+(1:3));
            F = couloumbForce(iPos,jPos);
            acci = F/cloud.atoms{i}.mass;
            accj = -F/cloud.atoms{j}.mass;
            ydotVals(ivecPos+(4:6)) = ydotVals(ivecPos+(4:6)) + acci;
            ydotVals(jvecPos+(4:6)) = ydotVals(jvecPos+(4:6)) + accj;
        end
    end
end

function e = energy(y,cloud)
%returns just energy from trap itself, not applied fields
e = 0;
potentialSet = cloud.potentialSet;
for i = 1:cloud.numAtoms  %each set is [x y z vx vy vz]
    mass = cloud.atoms{i}.mass;
    
    vecPos = (i-1) * 6;
    potPos = (i-1) * 3;

    potx = 0.5 * y(vecPos+1) * y(vecPos+1) * potentialSet{potPos+1}.forceK;  %trap potential
    poty = 0.5 * y(vecPos+2) * y(vecPos+2) * potentialSet{potPos+2}.forceK;  %trap potential
    potz = 0.5 * y(vecPos+3) * y(vecPos+3) * potentialSet{potPos+3}.forceK;  %trap potential
    
    kinx = 0.5 * mass * y(vecPos+4) * y(vecPos+4);
    kiny = 0.5 * mass * y(vecPos+5) * y(vecPos+5);
    kinz = 0.5 * mass * y(vecPos+6) * y(vecPos+6);
    
    e = e + potx + poty + potz + kinx + kiny + kinz;
end
if cloud.interacting
    for i = 1:cloud.numAtoms
        for j = i+1:cloud.numAtoms
            ivecPos = (i-1) * 6;
            jvecPos = (j-1) * 6;
            iPos = y(ivecPos+(1:3));
            jPos = y(jvecPos+(1:3));
            U = couloumbPot(iPos,jPos);
            e = e + U;
        end
    end
end

function e = individualEnergy(y,cloud)
%returns just energy from trap itself, not applied fields
e = zeros(1,cloud.numAtoms);
potentialSet = cloud.potentialSet;
for i = 1:cloud.numAtoms  %each set is [x y z vx vy vz]
    mass = cloud.atoms{i}.mass;
    
    vecPos = (i-1) * 6;
    potPos = (i-1) * 3;

    potx = 0.5 * y(vecPos+1) * y(vecPos+1) * potentialSet{potPos+1}.forceK;  %trap potential
    poty = 0.5 * y(vecPos+2) * y(vecPos+2) * potentialSet{potPos+2}.forceK;  %trap potential
    potz = 0.5 * y(vecPos+3) * y(vecPos+3) * potentialSet{potPos+3}.forceK;  %trap potential
    
    kinx = 0.5 * mass * y(vecPos+4) * y(vecPos+4);
    kiny = 0.5 * mass * y(vecPos+5) * y(vecPos+5);
    kinz = 0.5 * mass * y(vecPos+6) * y(vecPos+6);
    
    e(i) = potx + poty + potz + kinx + kiny + kinz;
    if cloud.interacting
        for j = i+1:cloud.numAtoms
            ivecPos = (i-1) * 6;
            jvecPos = (j-1) * 6;
            iPos = y(ivecPos+(1:3));
            jPos = y(jvecPos+(1:3));
            U = couloumbPot(iPos,jPos);
            e(i) = e(i) + U;
        end
    end
end

function U = couloumbPot(iPos,jPos)
    Q = 1.602176462e-19; %from https://en.wikipedia.org/wiki/Elementary_charge
    Ke = 8.987551787e9; %from https://en.wikipedia.org/wiki/Coulomb%27s_constant
    
    delR = jPos - iPos;
    r = norm(delR);
    U = (Ke * Q * Q / r);
    
function F = couloumbForce(iPos,jPos)
    Q = 1.602176462e-19; %from https://en.wikipedia.org/wiki/Elementary_charge
    Ke = 8.987551787e9; %from https://en.wikipedia.org/wiki/Coulomb%27s_constant
    
    delR = jPos - iPos;
    r = norm(delR);
    unitr = delR/r;
    forcemag = (Ke * Q * Q / r^2);
    F = -unitr * forcemag;


function cloud = initializeCloud(numAtoms,numFriction,cloudType,xPot,yPot,zPot)
if nargin < 2
    numFriction = numAtoms;
    cloudType = 'tightZ';
end
if numFriction == numAtoms
    isMixed = 'atoms';
else
    if numFriction ~= numAtoms && numFriction ~= 0
        isMixed = 'mixed';
    end
    if numFriction == 0
        isMixed = 'mlcls';
    end
end
atoms = cell(1,numAtoms);
if strcmp(isMixed,'atoms')
    for i = 1:numAtoms
        switch cloudType
            case 'tightZ'
                atoms{i} = initializeTightAtom(1);
            case 'looseZ'
                atoms{i} = initializeLooseAtom(1);
            case 'customZ'
                atoms{i} = initializeCustomAtom(1,xPot,yPot,zPot);
        end
    end     
end

if strcmp(isMixed,'mlcls')
    for i = 1:numAtoms
        switch cloudType
            case 'tightZ'
                atoms{i} = initializeTightAtom(2);
            case 'looseZ'
                atoms{i} = initializeLooseAtom(2);
            case 'customZ'
                atoms{i} = initializeCustomAtom(2,xPot,yPot,zPot);
        end
    end     
end

if strcmp(isMixed,'mixed')
    for i = 1:numFriction
        switch cloudType
            case 'tightZ'
                atoms{i} = initializeTightAtom(1);
            case 'looseZ'
                atoms{i} = initializeLooseAtom(1);
            case 'customZ'
                atoms{i} = initializeCustomAtom(1,xPot,yPot,zPot);
        end
    end
    for i = (numFriction+1):numAtoms
      switch cloudType
           case 'tightZ'
               atoms{i} = initializeTightAtom(2);
           case 'looseZ'
               atoms{i} = initializeLooseAtom(2);
           case 'customZ'
               atoms{i} = initializeCustomAtom(2,xPot,yPot,zPot);
      end
    end      
end

cloud.numAtoms = numAtoms;
cloud.numIons = numFriction;
cloud.numMolecules = numAtoms - numFriction;
cloud.atoms = atoms;
cloud.dt = 1e-8;
cloud.interacting = 1;
cloud = updateCloud(cloud);
cloud = tickleCloud(cloud);
cloud.isMixed = isMixed;
cloud.bufferGasBool = false;
cloud.bufferGasTime = [0 0];
cloud.bufferGasDensity = 0;
cloud.bufferGasTemp = 0;
cloud.numCollisions = 0;


function cloud = tickleCloud(cloud,r)
if nargin < 2
    r = 1e-4;
end
switch cloud.numAtoms
    case 1
        cloud.vector = r * [1 0.5 0.1 0 0 0];
        cloud = updateAtoms(cloud);
    case 2
        cloud.vector = r * [1 0.5 0.1 0 0 0 0.6 0.5 -0.2 0 0 0];
        cloud = updateAtoms(cloud);
    case 3
        cloud.vector = r * [1 0.5 0.1 0 0 0 0.6 0.5 -0.1 0 0 0 -0.6 0.1 0.4 0 0 0];
        cloud = updateAtoms(cloud);
    case 4
        cloud.vector = r * [1 0.5 0.1 0 0 0 0.6 0.5 -0.1 0 0 0 -0.6 0.1 0.4 0 0 0 0.1 0.1 0.4 0 0 0];
        cloud = updateAtoms(cloud);
    case 5
        cloud.vector = r * [1 0.5 0.1 0 0 0 0.6 0.5 -0.1 0 0 0 -0.6 0.1 0.4 0 0 0 0.1 0.1 0.4 0 0 0 0.1 0.5 1.5 0 0 0];
        cloud = updateAtoms(cloud);
    case 6
        cloud.vector = r * [1 0.5 0.1 0 0 0 0.6 0.5 -0.1 0 0 0 -0.6 0.1 0.4 0 0 0 0.1 0.1 0.4 0 0 0 0.1 0.5 1.5 0 0 0 -0.5 0 -1 0 0 0];
        cloud = updateAtoms(cloud);
    case 7
        cloud.vector = r * [1 0.5 0.1 0 0 0 0.6 0.5 -0.1 0 0 0 -0.6 0.1 0.4 0 0 0 0.1 0.1 0.4 0 0 0 0.1 0.5 1.5 0 0 0 -0.5 0 -1 0 0 0 -0.25 -0.25 -0.25 0 0 0];
        cloud = updateAtoms(cloud);
    case 8
        cloud.vector = r * [1 0.5 0.1 0 0 0 0.6 0.5 -0.1 0 0 0 -0.6 0.1 0.4 0 0 0 0.1 0.1 0.4 0 0 0 0.1 0.5 1.5 0 0 0 -0.5 0 -1 0 0 0 -0.25 -0.25 -0.25 0 0 0 0.5 -0.25 0.5 0 0 0];
        cloud = updateAtoms(cloud);
    case 9
        cloud.vector = r * [1 0.5 0.1 0 0 0 0.6 0.5 -0.1 0 0 0 -0.6 0.1 0.4 0 0 0 0.1 0.1 0.4 0 0 0 0.1 0.5 1.5 0 0 0 -0.5 0 -1 0 0 0 -0.25 -0.25 -0.25 0 0 0 0.5 -0.25 0.5 0 0 0 -0.1 0.5 1 0 0 0];
        cloud = updateAtoms(cloud);
    case 10
        cloud.vector = r * [1 0.5 0.1 0 0 0 0.6 0.5 -0.1 0 0 0 -0.6 0.1 0.4 0 0 0 0.1 0.1 0.4 0 0 0 0.1 0.5 1.5 0 0 0 -0.5 0 -1 0 0 0 -0.25 -0.25 -0.25 0 0 0 0.5 -0.25 0.5 0 0 0 -0.1 0.5 1 0 0 0 1 -0.5 1 0 0 0];
        cloud = updateAtoms(cloud);
    case 11
        cloud.vector = r * [1 0.5 0.1 0 0 0 0.6 0.5 -0.1 0 0 0 -0.6 0.1 0.4 0 0 0 0.1 0.1 0.4 0 0 0 0.1 0.5 1.5 0 0 0 -0.5 0 -1 0 0 0 -0.25 -0.25 -0.25 0 0 0 0.5 -0.25 0.5 0 0 0 -0.1 0.5 1 0 0 0 1 -0.5 1 0 0 0 0.2 1 -1 0 0 0];
        cloud = updateAtoms(cloud);
    case 12
        cloud.vector = r * [1 0.5 0.1 0 0 0 0.6 0.5 -0.1 0 0 0 -0.6 0.1 0.4 0 0 0 0.1 0.1 0.4 0 0 0 0.1 0.5 1.5 0 0 0 -0.5 0 -1 0 0 0 -0.25 -0.25 -0.25 0 0 0 0.5 -0.25 0.5 0 0 0 -0.1 0.5 1 0 0 0 1 -0.5 1 0 0 0 0.2 1 -1 0 0 0 -1 -1 0 0 0 0];
        cloud = updateAtoms(cloud);
    case 13
        cloud.vector = r * [1 0.5 0.1 0 0 0 0.6 0.5 -0.1 0 0 0 -0.6 0.1 0.4 0 0 0 0.1 0.1 0.4 0 0 0 0.1 0.5 1.5 0 0 0 -0.5 0 -1 0 0 0 -0.25 -0.25 -0.25 0 0 0 0.5 -0.25 0.5 0 0 0 -0.1 0.5 1 0 0 0 1 -0.5 1 0 0 0 0.2 1 -1 0 0 0 -1 -1 0 0 0 0 0.2 -0.4 -0.4 0 0 0];
        cloud = updateAtoms(cloud);
    case 14
        cloud.vector = r * [1 0.5 0.1 0 0 0 0.6 0.5 -0.1 0 0 0 -0.6 0.1 0.4 0 0 0 0.1 0.1 0.4 0 0 0 0.1 0.5 1.5 0 0 0 -0.5 0 -1 0 0 0 -0.25 -0.25 -0.25 0 0 0 0.5 -0.25 0.5 0 0 0 -0.1 0.5 1 0 0 0 1 -0.5 1 0 0 0 0.2 1 -1 0 0 0 -1 -1 0 0 0 0 0.2 -0.4 -0.4 0 0 0 0 0.45 0.5 0 0 0];
        cloud = updateAtoms(cloud);
    case 15
        cloud.vector = r * [1 0.5 0.1 0 0 0 0.6 0.5 -0.1 0 0 0 -0.6 0.1 0.4 0 0 0 0.1 0.1 0.4 0 0 0 0.1 0.5 1.5 0 0 0 -0.5 0 -1 0 0 0 -0.25 -0.25 -0.25 0 0 0 0.5 -0.25 0.5 0 0 0 -0.1 0.5 1 0 0 0 1 -0.5 1 0 0 0 0.2 1 -1 0 0 0 -1 -1 0 0 0 0 0.2 -0.4 -0.4 0 0 0 0 0.45 0.5 0 0 0 0.25 -0.5 0.25 0 0 0];
        cloud = updateAtoms(cloud);
    case 16
        cloud.vector = r * [1 0.5 0.1 0 0 0 0.6 0.5 -0.1 0 0 0 -0.6 0.1 0.4 0 0 0 0.1 0.1 0.4 0 0 0 0.1 0.5 1.5 0 0 0 -0.5 0 -1 0 0 0 -0.25 -0.25 -0.25 0 0 0 0.5 -0.25 0.5 0 0 0 -0.1 0.5 1 0 0 0 1 -0.5 1 0 0 0 0.2 1 -1 0 0 0 -1 -1 0 0 0 0 0.2 -0.4 -0.4 0 0 0 0 0.45 0.5 0 0 0 0.25 -0.5 0.25 0 0 0 -0.75 0.25 -0.05 0 0 0];
        cloud = updateAtoms(cloud);
    case 17
        cloud.vector = r * [1 0.5 0.1 0 0 0 0.6 0.5 -0.1 0 0 0 -0.6 0.1 0.4 0 0 0 0.1 0.1 0.4 0 0 0 0.1 0.5 1.5 0 0 0 -0.5 0 -1 0 0 0 -0.25 -0.25 -0.25 0 0 0 0.5 -0.25 0.5 0 0 0 -0.1 0.5 1 0 0 0 1 -0.5 1 0 0 0 0.2 1 -1 0 0 0 -1 -1 0 0 0 0 0.2 -0.4 -0.4 0 0 0 0 0.45 0.5 0 0 0 0.25 -0.5 0.25 0 0 0 -0.75 0.25 -0.05 0 0 0 -0.5 0.25 -0.6 0 0 0];
        cloud = updateAtoms(cloud);
    case 18
        cloud.vector = r * [1 0.5 0.1 0 0 0 0.6 0.5 -0.1 0 0 0 -0.6 0.1 0.4 0 0 0 0.1 0.1 0.4 0 0 0 0.1 0.5 1.5 0 0 0 -0.5 0 -1 0 0 0 -0.25 -0.25 -0.25 0 0 0 0.5 -0.25 0.5 0 0 0 -0.1 0.5 1 0 0 0 1 -0.5 1 0 0 0 0.2 1 -1 0 0 0 -1 -1 0 0 0 0 0.2 -0.4 -0.4 0 0 0 0 0.45 0.5 0 0 0 0.25 -0.5 0.25 0 0 0 -0.75 0.25 -0.05 0 0 0 -0.5 0.25 -0.6 0 0 0 -0.35 1.0 1.3 0 0 0];
        cloud = updateAtoms(cloud);
end


function cloud = updateAtoms(cloud)
    for i = 1:cloud.numAtoms
        vpoint = (6*(i-1));
    end
    cloud.atoms{i}.xyz = (cloud.vector(vpoint+1:vpoint+3));
    cloud.atoms{i}.vxvyvz = (cloud.vector(vpoint+4:vpoint+6));


function cloud = updateCloud(cloud)
%eventually include more 
%vector order is: [x1 y2 z2 vx2 vy2 vz x2 y2 z2 vx2 vy2 vz2...];
    
    cloud.varOrder = varOrder(cloud.numAtoms);
    cloud.vector = [];
    for i = 1:cloud.numAtoms
        cloud.vector = [cloud.vector cloud.atoms{i}.xyz cloud.atoms{i}.vxvyvz];
    end
    cloud.vector = cloud.vector';
    
    potentialSet = {}; %potentialSet = {x1 y1 z1 x2 y2 z2...}
    for i = 1:cloud.numAtoms
        potentialSet{end+1} = cloud.atoms{i}.potentialx; 
        potentialSet{end+1} = cloud.atoms{i}.potentialy;
        potentialSet{end+1} = cloud.atoms{i}.potentialz;
    end
    cloud.potentialSet = potentialSet; 
   
   
    
function order = varOrder(n)
    order = cell(1,6*n);
    for i = 1:n
        order{6*(i-1) + 1} = sprintf('x%d',i);
        order{6*(i-1) + 2} = sprintf('y%d',i);
        order{6*(i-1) + 3} = sprintf('z%d',i);
        order{6*(i-1) + 4} = sprintf('vx%d',i);
        order{6*(i-1) + 5} = sprintf('vy%d',i);
        order{6*(i-1) + 6} = sprintf('vz%d',i);
    end
    
function atom = initializeCustomAtom(type,xPot,yPot,zPot)
if nargin == 0
    type = 'Sr';
end
switch type
    case 'Sr'
        type = 1;
    case 'SrO'
        type = 2;
end
modType = mod(type,2);
atom.constants = makeConstants();
atom.xyz = [0 0 0];
atom.vxvyvz = [0 0 0];
atom.type = type;
switch modType
    case 1
        atom.mass = 88 * atom.constants.AMU;
        atom.charge = 1 * atom.constants.E;
        atom.potentialz = defaultPotentialZ(zPot,atom.mass);
        atom.potentialx = defaultPotential(xPot,0,atom.mass,'x',atom.potentialz.forceK);
        atom.potentialy = defaultPotential(yPot,0,atom.mass,'y',atom.potentialz.forceK);
        atom.friction = 1;
    case 0
        atom.mass = 89 * atom.constants.AMU;
        atom.charge = 1 * atom.constants.E;
        atom.potentialz = defaultPotentialZ((zPot * 0.994366),atom.mass); %    NEED TO CHANGE THIS TO VARY MASS: ZPOT * SQRT(88 / MASS OF MOLECULE)
        atom.potentialx = defaultPotential(xPot,0,atom.mass,'x',atom.potentialz.forceK);
        atom.potentialy = defaultPotential(yPot,0,atom.mass,'y',atom.potentialz.forceK);
        atom.friction = 0;
end
    
function atom = initializeTightAtom(type)
if nargin == 0
    type = 'Sr';
end
switch type
    case 'Sr'
        type = 1;
    case 'SrO'
        type = 2;
end
modType = mod(type,2);
atom.constants = makeConstants();
atom.xyz = [0 0 0];
atom.vxvyvz = [0 0 0];
atom.type = type;
switch modType
    case 1
        atom.mass = 88 * atom.constants.AMU;
        atom.charge = 1 * atom.constants.E;
        atom.potentialz = defaultPotentialZ(5e4,atom.mass);
        atom.potentialx = defaultPotential(100,0,atom.mass,'x',atom.potentialz.forceK);
        atom.potentialy = defaultPotential(100.1,0,atom.mass,'y',atom.potentialz.forceK);
        atom.friction = 1;
    case 0
        atom.mass = 89 * atom.constants.AMU; %CHANGED THIS TO TEST TICKLESCAN OF DARK ION MASS
        atom.charge = 1 * atom.constants.E;
        atom.potentialz = defaultPotentialZ(4.972e4,atom.mass);
        atom.potentialx = defaultPotential(100,0,atom.mass,'x',atom.potentialz.forceK);
        atom.potentialy = defaultPotential(100.1,0,atom.mass,'y',atom.potentialz.forceK);
        atom.friction = 0;
end

function atom = initializeLooseAtom(type)
if nargin == 0
    type = 'Sr';
end
switch type
    case 'Sr'
        type = 1;
    case 'SrO'
        type = 2;
end
modType = mod(type,2);
atom.constants = makeConstants();
atom.xyz = [0 0 0];
atom.vxvyvz = [0 0 0];
atom.type = type;
switch modType
    case 1
        atom.mass = 88 * atom.constants.AMU;
        atom.charge = 1 * atom.constants.E;
        atom.potentialz = defaultPotentialZ(1e4,atom.mass);
        atom.potentialx = defaultPotential(100,0,atom.mass,'x',atom.potentialz.forceK);
        atom.potentialy = defaultPotential(100.1,0,atom.mass,'y',atom.potentialz.forceK);
        atom.friction = 1;
    case 0
        atom.mass = 89 * atom.constants.AMU;
        atom.charge = 1 * atom.constants.E;
        atom.potentialz = defaultPotentialZ(0.912e4,atom.mass);
        atom.potentialx = defaultPotential(100,0,atom.mass,'x',atom.potentialz.forceK);
        atom.potentialy = defaultPotential(100.1,0,atom.mass,'y',atom.potentialz.forceK);
        atom.friction = 0;
end

function pot = defaultPotential(V,t,m,axis,zForceK,RF)
if nargin < 6
    RF = 2e6;
end
RF = RF * 2 * pi;
c = makeConstants();
q = 1 * c.E;
pot.angFreq = sqrt((q * V / ( (2^(1/2)) * m * RF * (2.8956e-3)^2))^2 - (0.5 * zForceK / m));
pot.freq = pot.angFreq / (2 * pi);  %in Hz
pot.mass = m;
pot.forceK = pot.mass * pot.angFreq^2;
% try to just input the radial average secular frequency, then move to RF
% from there, don't bother with calculating it from voltage
if axis == 'x'
    pot.trueForceK = (q * V * cos(RF * t) / ((2.8956e-3)^2)) - (0.5 * zForceK);
end
if axis == 'y'
    pot.trueForceK = (-1 * q * V * cos(RF * t) / ((2.8956e-3)^2)) - (0.5 * zForceK);
end
pot.axis = axis;
pot.voltage = V;
pot.descriptor = sprintf('%s potential, f = %3.1f MHz',pot.axis,pot.freq/1e6);
%function H = defaultHamiltonian()

function pot = defaultPotentialZ(f,m)
pot.freq = f;  %in Hz
pot.angFreq = pot.freq * 2 * pi;
pot.mass = m;
pot.forceK = pot.mass * pot.angFreq^2;
pot.axis = 'z';
pot.descriptor = sprintf('%s potential, f = %3.1f MHz',pot.axis,pot.freq/1e6);
%function H = defaultHamiltonian()

function c = makeConstants
%example: twospins(1,0.5);
c.H = 6.626068766e-34;                  %in SI
c.HBAR = c.H/(2 * pi);                  %in SI
c.C = 2.99792458e8;                     %in m/s
c.E = 1.602176462e-19;                  %in Coulombs
c.BOHRMAGNETON = 9.274009994e-24;       %in Joules/Tesla
c.BOHRMAGNETONHZ = c.BOHRMAGNETON/c.H;  %in Hz/Tesla
c.BOLTZMANN = 1.380650e-23;             %in SI
c.AMU = 1.66e-27;


function addText(tstring,ABOVE)
if (nargin == 1)
    ABOVE = 0;
end
%adds text to the current plot.
if ~isempty(tstring)
    a = xlim;
    b = ylim;
    xstart = a(1) + 0.1 * (a(2) - a(1));
    ystart = b(1) + 0.55 * (b(2) - b(1));
    if ABOVE == 1
        ystart = b(1) + 0.8 * (b(2) - b(1));
    end
    if ABOVE == -1
        ystart = b(1) + 0.2 * (b(2) - b(1));
       
    end
    if ABOVE == -2
        ystart = b(1) + 0.05 * (b(2) - b(1));
       
    end
        text(xstart,ystart,tstring,'FontSize',10,'FontWeight','normal');
       
end


