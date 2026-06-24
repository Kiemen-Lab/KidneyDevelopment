function train_segmentation_model(pthDL)
% pthDL: path to save the model, should already contain model metadata
datafile=[pthDL,'net.mat'];
load(datafile,'sxy','classNames','nm');

variableInfo = who('-file',datafile);
if ismember('net', variableInfo)
    error([' A network has already been trained for model ',nm,...
        '. Choose a new model name to retrain.'])
end

disp(' ')
disp('Starting model training...')
disp(' ')

classes=1:length(classNames);

nmim='im\';
nmlabel='label\';

pthTrain=[pthDL,'training\'];
pthVal=[pthDL,'validation\'];

% 1 make training data
TrainHE=[pthTrain,nmim];
Trainlabel=[pthTrain,nmlabel];
imdsTrain = imageDatastore(TrainHE);
pxdsTrain = pixelLabelDatastore(Trainlabel,classNames,classes);
pximdsTrain = pixelLabelImageDatastore(imdsTrain,pxdsTrain); %'DataAugmentation',augmenter);
tbl = countEachLabel(pxdsTrain);

% make validation data
ValHE=[pthVal,nmim];
Vallabel=[pthVal,nmlabel];
imdsVal = imageDatastore(ValHE);
pxdsVal = pixelLabelDatastore(Vallabel,classNames,classes);
pximdsVal = pixelLabelImageDatastore(imdsVal,pxdsVal); %'DataAugmentation',augmenter);

options = trainingOptions('adam',...    % stochastic gradient descent solver
    'MaxEpochs',8,...                   % maximum amount of times we let the model see each training tile
    'MiniBatchSize',4,...               % datapoints per 'mini-batch' - ideally a small power of 2 (32, 64, 128, or 256)
    'Shuffle','every-epoch',...         % reallocate mini-batches each epoch (so min-batches are new mixtures of data)
    'ValidationData',pximdsVal,...      % image datastore object of validation data
    'ValidationPatience',7,...          % stop training when validation data doesn't improve for __ iterations 5
    'InitialLearnRate',0.0005,...       % amount that network is allowed to update its weights in each minibatch
    'LearnRateSchedule','piecewise',... % drop learning rate during training to prevent overfitting
    'LearnRateDropPeriod',1,...         % drop learning rate every _ epochs
    'LearnRateDropFactor',0.75,...      % multiply learning rate by this factor to drop it
    'ValidationFrequency',128,...       % initial loss should be -ln( 1 / # classes )
    'ExecutionEnvironment','gpu',...    % train on gpu
    'Plots','training-progress');        % view progress while training
    %'OutputFcn', @(info)SaveTrainingPlot(info,pthDL)); % save training progress as image

% Design network
numclass = numel(classes);
imageFreq = tbl.PixelCount ./ tbl.ImagePixelCount;
classWeights = median(imageFreq) ./ imageFreq;

lgraph = deeplabv3plusLayers([sxy sxy 3],numclass,"resnet50");
pxLayer = pixelClassificationLayer('Name','labels','Classes',tbl.Name,'ClassWeights',classWeights);
lgraph = replaceLayer(lgraph,"classification",pxLayer);

% train
[net, info] = trainNetwork(pximdsTrain,lgraph,options);
datafile=[pthDL,'net.mat'];
if exist(datafile,'file')
    save(datafile,'net','info','-append');
else
    save(datafile,'net','info');
end

end

function stop=SaveTrainingPlot(info,pthSave)
    stop=false;  %prevents this function from ending trainNetwork prematurely
    if info.State == "done"
        currentfig = findall(groot,'Type','Figure');
        savefig(currentfig,[pthSave,'training_process.fig'])
        %exportapp(currentfig,[pthSave,'training_process.png'])
    end
end
