%% Simulation Constants
S_rate = 8000;  %sampling rate 8kHz
S_per_frm = 80; %samples per frame 80
T = 30;         %time window for scope 30s

%% Filter Design
band_pass = designfilt('bandpassiir', ...       % Response type
       'StopbandFrequency1',30, ...    % Frequency constraints 30Hz-3.5kHz
       'PassbandFrequency1',40, ...
       'PassbandFrequency2',3000, ...
       'StopbandFrequency2',3500, ...
       'StopbandAttenuation1',10, ...   % Magnitude constraints
       'PassbandRipple',1, ...
       'StopbandAttenuation2',50, ...
       'DesignMethod','ellip', ...      % Design method
       'MatchExactly','passband', ...   % Design method options
       'SampleRate',S_rate)  ;             % Sample rate


%% Audio Input

audioSource = audioDeviceReader('OutputDataType', 'single', ...     %reading audio from mic
                                'NumChannels', 1, ...               %mono
                                'SamplesPerFrame', S_per_frm, ...   %80 samples per frame
                                'SampleRate', S_rate);              %8kHz sample rate
%% Display Settings
scope = dsp.TimeScope(2, 'SampleRate', [S_rate/S_per_frm, S_rate], ...  %frames per second, sample rate
                         'BufferLength', S_rate*T, ...                  % buffer length set to show full data set
                         'YLimits', [-0.3 1.1], ...
                         'ShowGrid', true, ...
                         'Title','Decision speech and speech data', ...
                         'TimeSpanOverrunAction','Scroll',...           %Continues to scroll after time span is surpassed
                         'TimeSpan',T);
%% Initialization                 
VAD_cst_param = vadInitCstParams;% Initialize VAD parameters
clear vadG729
%% Implementation

numTSteps = T*100;% Run for T seconds
while(numTSteps)
  
  speech = audioSource(); % Retrieve 10 ms of speech data from the audio recorder
  
  speech = filter(band_pass, speech);%filter the data
 
  decision = vadG729(speech, VAD_cst_param); % Call the VAD algorithm
  
  scope(decision, speech); % Plot speech frame and decision: 1 for speech, 0 for silence
  %numTSteps = numTSteps - 1; %include line if you want to limit simulation
end
release(scope);