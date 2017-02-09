%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%   SKA Africa                                                                %
%   http://www.kat.ac.za                                                      %
%   Copyright (C) 2013 Andrew Martens (andrew@ska.ac.za)                      %
%                                                                             %
%   This program is free software; you can redistribute it and/or modify      %
%   it under the terms of the GNU General Public License as published by      %
%   the Free Software Foundation; either version 2 of the License, or         %
%   (at your option) any later version.                                       %
%                                                                             %
%   This program is distributed in the hope that it will be useful,           %
%   but WITHOUT ANY WARRANTY; without even the implied warranty of            %
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             %
%   GNU General Public License for more details.                              %
%                                                                             %
%   You should have received a copy of the GNU General Public License along   %
%   with this program; if not, write to the Free Software Foundation, Inc.,   %
%   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.               %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function to find channel response of casper pfb chain
function[dout, result] = casper_pfb_fir_async_generic_channel_response_dsim_data(varargin)

  %load('cwg0_out_208984375.mat');
  %load('cwg1_out_208929882.mat');
  %load('cwg1_out_208983375.mat');
  %load('cwg1_out_53508358.mat');
  %load('cwg1_out_54926318.36_1e6.mat');
  load('cwg0_100e6_length16e6.mat');
  dsim_signal = CWG0;
  clear CWG0;
  
  %load('cwg0_out_54753906.25_1e6.mat');
  %load('cwg1_out_54840112.3_1e6.mat');
  %load('cwg1_out_804552000_1e6.mat');
  %load('cwg1_out_54552000_1e6.mat');
  
  %cwg1_out = cwg1(392:end);
  %clear cwg1;
   
  % Import Data
  %load('cwg_428074893.mat');
  %load('cwg_428000401.mat');

  %cwg1_out = cwg1_out(392:end);
  clear cwg0_out;
  
  offset = 0; %500;
  %dsim_signal = cwg1_out;
  
  %tic  
  warning off Simulink:Engine:OutputNotConnected
  ut_log_group = 'casper_pfb_fir_generic_channel_response_debug';

  utlog('entering', {'trace', ut_log_group});

  result = -1; 
  data = []; 

  defaults = { ...
    'n_spectra', 2^0, ...             % the amount of averaging to perform   (was 8)
    'n_taps', 8, ...                % number of taps in pfb_fir
    'fft_stages', 13, ...           % size of FFT
    'window_type', 'hamming', ...   % window type
    'fwidth', 1, ...                % relative width of main lobe in channel response
    'plots', 'on', ...             % output a plot of the response when done
    'offset', 0, ...                % offset of start frequency in fft bins
    'freq_step', 1/2, ...          % frequency step for each iteration in fraction of an fft_bin     1/15
    'span', 4, ...                  % frequency span over which to iterate
    'amp', 0.5, ...                 % amplitude of sinusoid
    'noise_power', 1/2^15, ...      % power in guassian white noise to add to input sinusoid
    'dc_power', 0, ...              % dc power to add to input sinusoid
  };

  args = {varargin{:}, 'defaults', defaults};
  [n_taps, temp, results(1)]                 = utpar_get({args, 'n_taps'});
  [fft_stages, temp, results(2)]             = utpar_get({args, 'fft_stages'});
  [window_type, temp, results(3)]            = utpar_get({args, 'window_type'});
  [fwidth, temp, results(4)]                 = utpar_get({args, 'fwidth'});
  [offset, temp, results(5)]                 = utpar_get({args, 'offset'});
  [fs, temp, results(6)]                     = utpar_get({args, 'freq_step'});
  [span, temp, results(7)]                   = utpar_get({args, 'span'});
  [amp, temp, results(8)]                    = utpar_get({args, 'amp'});
  [noise_power, temp, results(9)]            = utpar_get({args, 'noise_power'});
  [dc_power, temp, results(10)]              = utpar_get({args, 'dc_power'});
  [plots, temp, results(11)]                 = utpar_get({args, 'plots'});
  [n_spectra, temp, results(12)]             = utpar_get({args, 'n_spectra'});

  if ~isempty(find(results ~= 0)),
      utlog('error getting parameters from varargin',{'error', ut_log_group});
      return;
  end

  %%%%%%%%%%%%%%%%%%%
  % system settings %
  %%%%%%%%%%%%%%%%%%%
  
  fft_shift             = 2^fft_stages-1;
  n_inputs              = 3;
  n_bits_in             = 10;
  n_bits_out            = 18;

  %%%%%%%%%%%%%%%%%%%%
  % pfb_fir settings %   
  %%%%%%%%%%%%%%%%%%%%

  pfb_fir_model_name    = 'pfb_fir_generic_async_test_sim4';
  pfb_fir_block_name    = 'pfb_fir_generic';
  pfb_fir_type          = 'default';
  pfb_fir_block_version = '0'; 

  pfb_size              = fft_stages;
  pfb_n_bits_in         = n_bits_in;
  pfb_n_bits_out        = n_bits_out;
  pfb_coeff_bit_width   = 18;
  pfb_quantization      = 'Round  (unbiased: Even Values)';
  pfb_async             = 'on';
  
  

  %%%%%%%%%%%%%%%%
  % fft settings %   
  %%%%%%%%%%%%%%%%

  fft_model_name        = 'fft_async_test_sim_FD';
  
  fft_block_name        = 'fft_wideband_real';
  fft_type              = fft_block_name;
  fft_block_version     = '1';

  fft_n_bits_in         = 18;
  fft_coeff_bit_width   = 18;
  fft_quantization      = 'Round  (unbiased: Even Values)';
  fft_overflow          = 'Wrap';
  fft_async             = 'on';

  %%%%%%%%%%%%%%%%%%%%%%
  % pfb_fir parameters %
  %%%%%%%%%%%%%%%%%%%%%%

  pfb_fir_parameters = { ...
      'TotalTaps', n_taps, ...
      'CoeffBitWidth', pfb_coeff_bit_width, ...
      'BitWidthIn', pfb_n_bits_in, ...
      'BitWidthOut', pfb_n_bits_out, ...
      'n_inputs', n_inputs, ...
      'PFBSize', fft_stages, ...
      'quantization', pfb_quantization, ...
      'async', pfb_async, ...
  };

  pfb_fir_block = { ...
    'name', pfb_fir_block_name, ...
    'type', pfb_fir_type, ...
    'version', pfb_fir_block_version, ...
    'parameters', pfb_fir_parameters, ...
  }; %block

  [pfb_fir_model, result] = utpfb_fir_generic_bbox_modelify( ...
    'name', pfb_fir_model_name, ...
    'block', pfb_fir_block);
  if result ~= 0,
      utlog('error creating black box pfb_fir model', {'error', ut_log_group});
      return;
  end

  %erect the model if it doesn't already exist
  sys = find_system('type', 'block_diagram', 'name', pfb_fir_model_name);
  if isempty(sys),
    result = utmodel_erect(pfb_fir_model{:});
    if result ~= 0,
      utlog('error erecting pfb_fir model', {'error', ut_log_group});
      return;
    end
  else %update if it exists
  %TODO should use utmodel_update

  end

  %%%%%%%%%%%%%%%%%%
  % fft parameters %
  %%%%%%%%%%%%%%%%%%

  fft_parameters = { ...
      'coeff_bit_width', fft_coeff_bit_width, ...
      'input_bit_width', fft_n_bits_in, ...
      'bin_pt_in', fft_n_bits_in-1, ...
      'n_inputs', n_inputs, ...
      'FFTSize', fft_stages, ...
      'quantization', fft_quantization, ...
      'overflow', fft_overflow, ...
      'async', fft_async,...
  };

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %construct black box model for erecting
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  fft_block = { ...
    'name', fft_block_name, ...
    'type', fft_type, ...
    'version', fft_block_version, ...
    'parameters', fft_parameters, ...
  }; %block

  [fft_model, result] = utfft_bbox_modelify( ...
    'name', fft_model_name, ...
    'block', fft_block);
  if result ~= 0,
      utlog('error creating fft black box model', {'error', ut_log_group});
      return;
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % erect the model if it doesn't already exist %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  sys = find_system('type', 'block_diagram', 'name', fft_model_name);
  if isempty(sys),
    result = utmodel_erect(fft_model{:});
    if result ~= 0,
      utlog('error erecting fft model', {'error', ut_log_group});
      return;
    end
  else, %update if exists
  %TODO should use utmodel_update

  end
  
  response = zeros(2^(fft_stages-1), (span/fs));

  
  start_idx = 1;
  end_idx = 0;
  
  %get response
  for index = 1:1:floor(length(dsim_signal)/(2^fft_stages*n_taps)),
    tic
    index
    
    if index == 20
        a = 1;
    elseif index == 40
        a = 1;
    elseif index == 60
        a = 1;
    elseif index == 80
        a = 1; 
    end
     
    start_idx = end_idx+1;
    end_idx = index*2^fft_stages*n_taps;  
    
    if (length(dsim_signal)>=(2^fft_stages*n_taps))
        din(index,:) = dsim_signal(start_idx:end_idx)';    
    else
        error('dsim_signal too short')
    end

    dvalid = ones(1,2^fft_stages*n_taps*n_spectra);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % push data through pfb_fir %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [pfb_fir_dout, pfb_fir_dvalid, result] = utpfb_fir_generic_bbox_sim( pfb_fir_model{:}, ...
        'name', pfb_fir_block_name, 'din', din(index,:)', 'dvalid', dvalid, 'debug', 'off');     %MOD - dvalid
    
    if result ~= 0,
      utlog('error black box simulating pfb_fir', {'error', ut_log_group});
      return;
    end 
   
    [r,c] = size(pfb_fir_dout);
    pfb_fir_dout_rs = reshape(pfb_fir_dout, r*c, 1);
    
    pfb_fir(:,index) = pfb_fir_dout_rs(1:2^fft_stages*n_spectra);

    %Temp
    pfb_fir_dvalid = ones(2^fft_stages*n_spectra,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % push data through fft %
    %%%%%%%%%%%%%%%%%%%%%%%%%
   
    %utlog('simulating fft response', {ut_log_group});
    %[fft_dout, result] = utfft_bbox_sim( fft_model{:}, ...
    %    'name', fft_block_name, 'fft_shift', fft_shift, 'din', pfb_fir_dout, 'dvalid', pfb_fir_dvalid, 'debug', 'on');
    %if result ~= 0,
    %  utlog('error black box simulation of fft', {'error', ut_log_group});
    %  return;
    %end 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perform averaging on the result %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %fft_dout_sum = sum(abs(fft_dout(:,1:n_spectra)),2); %sum along rows
 
    %response(:, (freq_index-offset/fs)+1) = fft_dout_sum(1:2^(fft_stages-1),1)./n_spectra;
    %response(:, index) = fft_dout_sum(1:2^(fft_stages-1),1);
    
    toc
  end %for

  %dout = response'; %flip so that column contains fft bin response 

  result = 0;

  
  %[SNR,SFDR] = Adj_Channel_suppression(dout,span);
  %[SNR,SFDR,ch_deviation,SFDR_locations,variance_PFB] = PFB_statistics(dout,span);
  
  if strcmp(plots, 'on'),
    %plot results
    figure(2);
    plot([offset:fs:(offset+span)-fs], 20*log10(dout));
    title('pfb chanel response');
    
    %Add SNR to plot 
    str_dB = 'dB';
    str_SNR = 'SNR = ';
    str_SFDR = 'SFDR = ';
    xlabel('Channel');
    ylabel(str_dB);
    
    %strmin = [num2str(ceil(SNR))];
    %text((0:length(SNR)-1),(-10*ones(1,length(SNR))),strcat(str_SNR,strmin,str_dB),'HorizontalAlignment','left');
    
    %Add SFDR to plot 
    %strSFDR = [num2str(ceil(SFDR))];
    %text((0:length(SFDR)-1),(-15*ones(1,length(SFDR))),strcat(str_SFDR,strSFDR,str_dB),'HorizontalAlignment','left');
  end
 
%     if strcmp(plots, 'on'),
%         %plot results
%         figure(3);
%         plot([offset:fs:(offset+span)-fs], 20*log10(dout));
%         title('pfb channel response');
% 
%         hold on;
%         plot([offset:fs:(offset+span)-fs], 20*log10(SFDR_locations),'*');
% 
% 
%         %Add SNR to plot 
%         str_dB = 'dB';
%         str_SNR = 'SNR = ';
%         str_SFDR = 'SFDR = ';
%         str_Var = 'Var = ';
%         str_ChDev = 'Ch Dev = ';
%         xlabel('Channel');
%         ylabel(str_dB);
% 
%         %Add SNR data to plot 
%         strmin = [num2str(ceil(SNR))];
%         text((0:length(SNR)-1),(-10*ones(1,length(SNR))),strcat(str_SNR,strmin,str_dB),'HorizontalAlignment','left');
% 
%         %Add SFDR data to plot 
%         strSFDR = [num2str(ceil(SFDR))];
%         text((0:length(SFDR)-1),(-14*ones(1,length(SFDR))),strcat(str_SFDR,strSFDR,str_dB),'HorizontalAlignment','left');
% 
%         %Add variance info to the figure
%         strVariance = [num2str(variance_PFB)];
%         text((0:length(variance_PFB)-1),(-6*ones(1,length(variance_PFB))),strcat(str_Var,strVariance),'HorizontalAlignment','left');
% 
%         %Add channel deviation info to the figure
%         strChDev = [num2str(ch_deviation)];
%         text((0:length(strChDev)-1),(-2*ones(1,length(strChDev))),strcat(str_ChDev,strChDev,str_dB),'HorizontalAlignment','left');
%     end

  utlog('exiting', {'trace', ut_log_group});
  
  toc
end %function
