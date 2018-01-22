function vad_flag = vadG729(speech, VAD_cst_param)
% VADG729 Implement the Voice Activity Detection Algorithm.
% Note that although G.729 VAD operates on pre-processed speech data, this
% function is a standalone version, i.e. the pre-processing (highpass
% filtering and linear predictive analysis) is also included.
%
% This function is in support of the 'G.729 Voice Activity Detection'
% example and may change in a future release

% Copyright 2015 The MathWorks, Inc.

%% Algorithm Components Initialization
persistent HPF AC LEV1 LEV2 LPC2LSF ZCD VAD_var_param

if isempty(HPF)
    % Create a IIR digital filter used for pre-processing
    HPF = dsp.BiquadFilter('SOSMatrix', VAD_cst_param.HPF_sos);                    
                    
    % Create an autocorrelator and set its properties to compute the lags 
    % in the range [0:NP].
    AC = dsp.Autocorrelator('MaximumLagSource', 'Property', ...
                            'MaximumLag', VAD_cst_param.NP);
                        
    % Create a Levinson solver which compute the reflection coefficients 
    % from auto-correlation function using the Levinson-Durbin recursion. 
    % The first object is configured to output polynomial coefficients and 
    % the second object is configured to output reflection coefficients.
    LEV1 = dsp.LevinsonSolver('AOutputPort', true, ...
                             'KOutputPort', false);

    LEV2 = dsp.LevinsonSolver('AOutputPort', false, ...
                              'KOutputPort', true);
                         
    % Create a converter from linear prediction coefficients (LPC) to line 
    % spectral frequencies (LSF)
    LPC2LSF = dsp.LPCToLSF;    
    
    % Create a zero crossing detector
    ZCD = dsp.ZeroCrossingDetector;
    
    % initialize variable parameters
    VAD_var_param.window_buffer = single(zeros(VAD_cst_param.L_WINDOW - VAD_cst_param.L_FRAME, 1));
    VAD_var_param.frm_count = single(0);
    VAD_var_param.MeanLSF = single(zeros(VAD_cst_param.M, 1));
    VAD_var_param.MeanSE = single(0);
    VAD_var_param.MeanSLE = single(0);
    VAD_var_param.MeanE = single(0);
    VAD_var_param.MeanSZC = single(0);
    VAD_var_param.count_sil = single(0);
    VAD_var_param.count_update = single(0);
    VAD_var_param.count_ext = single(0);
    VAD_var_param.less_count = single(0);
    VAD_var_param.flag = single(1);
    VAD_var_param.prev_markers = single([1, 1]);
    VAD_var_param.prev_energy = single(0);
    VAD_var_param.Prev_Min = single(Inf);
    VAD_var_param.Next_Min = single(Inf);
    VAD_var_param.Min_buffer = single(Inf * ones(1, VAD_cst_param.No));    
end

%% Constants Initialization
NOISE = single(0);
VOICE = single(1);
v_flag = single(0);%#ok
VAD_var_param.frm_count = single(VAD_var_param.frm_count + 1);

%% Pre-processing
% Filter the speech frame: this high-pass filter serves as a precaution 
% against undesired low-frequency components.
speech_hp = step(HPF, 32768*speech);

% Store filtered data to the pre-processed speech buffer
speech_buf = [VAD_var_param.window_buffer; speech_hp];

% LPC analysis

% Windowing of signal
speech_win = VAD_cst_param.hamwindow .* speech_buf;

% Autocorrelation
r = step(AC, speech_win) .* VAD_cst_param.lagwindow;

% LSF
a = step(LEV1, r(1:VAD_cst_param.M+1));
LSF = step(LPC2LSF, a) / (2 * pi);

% Reflection coefficients
rc = step(LEV2, r(1:3));

%% VAD starts here

%% Parameters extraction

% Full-band energy
Ef = 10 * log10(r(1) / VAD_cst_param.L_WINDOW);

% Low-band energy
El = r(1) * VAD_cst_param.lbf_corr(1) + 2 * sum(r(2:end) .* VAD_cst_param.lbf_corr(2:end));
El = 10 * log10(El / VAD_cst_param.L_WINDOW);

% Spectral Distorsion
SD = sum((LSF-VAD_var_param.MeanLSF).^2);

% Zero-crossing rate
idx = VAD_cst_param.L_WINDOW - VAD_cst_param.L_NEXT - VAD_cst_param.L_FRAME + 1;
ZC = double(step(ZCD, speech_buf(idx:idx+VAD_cst_param.L_FRAME)))/VAD_cst_param.L_FRAME;

% Long-term minimum energy
VAD_var_param.Next_Min = min(Ef, VAD_var_param.Next_Min);
Min = min(VAD_var_param.Prev_Min, VAD_var_param.Next_Min);
if (mod(VAD_var_param.frm_count, single(8)) == 0)
    VAD_var_param.Min_buffer = [VAD_var_param.Min_buffer(2:end), VAD_var_param.Next_Min];
    VAD_var_param.Prev_Min = min(VAD_var_param.Min_buffer);
    VAD_var_param.Next_Min = single(Inf);
end

if (VAD_var_param.frm_count < VAD_cst_param.Ni)
    %% Initialization of running averages if frame number is less than 32
    if (Ef < 35)
        VAD_var_param.less_count = VAD_var_param.less_count + 1;
        marker = NOISE;
    else
        % include only the frames that have an energy Ef greater than 21
        marker = VOICE;
        NE = (VAD_var_param.frm_count - 1) - VAD_var_param.less_count;
    
        VAD_var_param.MeanE = (VAD_var_param.MeanE * NE + Ef) / (NE+1);
        VAD_var_param.MeanSZC = (VAD_var_param.MeanSZC * NE + ZC) / (NE+1);
        VAD_var_param.MeanLSF = (VAD_var_param.MeanLSF * NE + LSF) / (NE+1);
    end

else 
    %% Start calculating the chararcteristic energies of background noise
    if (VAD_var_param.frm_count == VAD_cst_param.Ni)
        VAD_var_param.MeanSE = VAD_var_param.MeanE - 10;
        VAD_var_param.MeanSLE = VAD_var_param.MeanE - 12;
    end
    
    % Difference measures between current frame parameters and running
    % averages of background noise characteristics
    dSE = VAD_var_param.MeanSE - Ef;
    dSLE = VAD_var_param.MeanSLE - El;
    dSZC = VAD_var_param.MeanSZC - ZC;

    %% Initial VAD decision
    if (Ef < 21)
        marker = NOISE;
    else
        marker = vad_decision(dSLE, dSE, SD, dSZC);
    end

    v_flag = single(0);
    
    %% Voice activity decision smoothing
    % from energy considerations and neighbouring past frame decisions
    
    % Step 1
    if ((VAD_var_param.prev_markers(1) == VOICE) && (marker == NOISE) ...
        && (Ef > VAD_var_param.MeanSE + 2) && (Ef > 21))
        marker = VOICE;
        v_flag = single(1);
    end
    
    % Step 2
    if (VAD_var_param.flag == 1)
        if ((VAD_var_param.prev_markers(2) == VOICE) ...
            && (VAD_var_param.prev_markers(1) == VOICE) ...
            && (marker == NOISE) ...
            && (abs(Ef - VAD_var_param.prev_energy) <= 3))
        
            VAD_var_param.count_ext = VAD_var_param.count_ext + 1;
            marker = VOICE;
            v_flag = single(1);
           
            if (VAD_var_param.count_ext <= 4)
                VAD_var_param.flag = single(1);
            else
                VAD_var_param.count_ext = single(0);
                VAD_var_param.flag = single(0);
            end
        end
    else
        VAD_var_param.flag = single(1);
    end

    if (marker == NOISE)
        VAD_var_param.count_sil = VAD_var_param.count_sil + 1;
    end
  
    % Step 3    
    if ((marker == VOICE) && (VAD_var_param.count_sil > 10) ...
        && (Ef - VAD_var_param.prev_energy <= 3))
        marker = NOISE;
        VAD_var_param.count_sil = single(0);
    end
    
    if (marker == VOICE)
        VAD_var_param.count_sil = single(0);
    end

    % Step 4
    if ((Ef < VAD_var_param.MeanSE + 3) && (VAD_var_param.frm_count > VAD_cst_param.No) ...
        && (v_flag == 0) && (rc(2) < 0.6))
        marker = NOISE;
    end
    
    %% Update running averages only in the presence of background noise
    
    if ((Ef < VAD_var_param.MeanSE + 3) && (rc(2) < 0.75) && (SD < 0.002532959))
        VAD_var_param.count_update = VAD_var_param.count_update + 1;
        % Modify update speed coefficients
        if (VAD_var_param.count_update < VAD_cst_param.INIT_COUNT)
            COEF = single(0.75);
            COEFZC = single(0.8);
            COEFSD = single(0.6);
        elseif (VAD_var_param.count_update < VAD_cst_param.INIT_COUNT + 10)
            COEF = single(0.95);
            COEFZC = single(0.92);
            COEFSD = single(0.65);
        elseif (VAD_var_param.count_update < VAD_cst_param.INIT_COUNT + 20)
            COEF = single(0.97);
            COEFZC = single(0.94);
            COEFSD = single(0.70);
        elseif (VAD_var_param.count_update < VAD_cst_param.INIT_COUNT + 30)
            COEF = single(0.99);
            COEFZC = single(0.96);
            COEFSD = single(0.75);
        elseif (VAD_var_param.count_update < VAD_cst_param.INIT_COUNT + 40)
            COEF = single(0.995);
            COEFZC = single(0.99);
            COEFSD = single(0.75);
        else
            COEF = single(0.995);
            COEFZC = single(0.998);
            COEFSD = single(0.75);
        end

        % Update mean of parameters LSF, SE, SLE, SZC
        VAD_var_param.MeanSE = COEF * VAD_var_param.MeanSE + (1-COEF) * Ef;
        VAD_var_param.MeanSLE = COEF * VAD_var_param.MeanSLE + (1-COEF) * El;
        VAD_var_param.MeanSZC = COEFZC * VAD_var_param.MeanSZC + (1-COEFZC) * ZC;
        VAD_var_param.MeanLSF = COEFSD * VAD_var_param.MeanLSF + (1-COEFSD) * LSF;        
    end

    if ((VAD_var_param.frm_count > VAD_cst_param.No) && ((VAD_var_param.MeanSE < Min) &&  (SD < 0.002532959)) || (VAD_var_param.MeanSE > Min + 10 ))
        VAD_var_param.MeanSE = Min;
        VAD_var_param.count_update = single(0);
    end
end

%% Update parameters for next frame
VAD_var_param.prev_energy = Ef;
VAD_var_param.prev_markers = [marker, VAD_var_param.prev_markers(1)];

idx = VAD_cst_param.L_FRAME + 1;
VAD_var_param.window_buffer = speech_buf(idx:end);

%% Return final decision
vad_flag = marker;

return

function dec = vad_decision(dSLE, dSE, SD, dSZC)

% Active voice decision using multi-boundary decision regions in the space
% of the 4 difference measures
a = single([0.00175, -0.004545455, -25, 20, 0, ...
     8800, 0, 25, -29.09091, 0, ...
     14000, 0.928571, -1.5, 0.714285]);

b = single([0.00085, 0.001159091, -5, -6, -4.7, ...
     -12.2, 0.0009, -7.0, -4.8182, -5.3, ...
     -15.5, 1.14285, -9, -2.1428571]);

dec = single(0);

if SD > a(1)*dSZC+b(1)
    dec = single(1);
    return;
end

if SD > a(2)*dSZC+b(2)
    dec = single(1);
    return;
end

if dSE < a(3)*dSZC+b(3)
    dec = single(1);
    return;
end

if dSE < a(4)*dSZC+b(4)
    dec = single(1);
    return;
end

if dSE < b(5)
    dec = single(1);
    return;
end
    
if dSE < a(6)*SD+b(6)
    dec = single(1);
    return;
end

if SD > b(7)
    dec = single(1);
    return;
end

if dSLE < a(8)*dSZC+b(8)
    dec = single(1);
    return;
end

if dSLE < a(9)*dSZC+b(9)
    dec = single(1);
    return;
end

if dSLE < b(10)
    dec = single(1);
    return;
end

if dSLE < a(11)*SD+b(11)
    dec = single(1);
    return;
end

if dSLE > a(12)*dSE+b(12)
    dec = single(1);
    return
end

if dSLE < a(13)*dSE+b(13)
    dec = single(1);
    return;
end

if dSLE < a(14)*dSE+b(14)
    dec = single(1);
    return;
end

return
