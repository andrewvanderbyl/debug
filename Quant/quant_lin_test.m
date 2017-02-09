% This m-file works inconjuction with the Silumlink test model
% 'dsim_quant_1'. Run the model then run this script to see the results.
% The model will divide the coefficient scaling by 2 every 4096 ticks for
% a total of 10 interations. This script will split the output into a 
% matrix where each column is ticks_per_scale in length and will compute 
% the mean and plot it on a log scale.


% -------------- %
% Linearity Test %
% -------------- %

% Test multiplier output in quant
% -------------------------------
initial_latency = 607;
ticks_per_scale = 4096;
plot_iterations = 8;
% hist_8b = 128;
% hist_10b = 2^7;
% hist_54b = 2^12;

%vec_mult = q_m_10b_2.signals.values(initial_latency+1:initial_latency + ticks_per_scale*(floor(length(q_m_10b_2.signals.values)/ticks_per_scale)-1));
vec_mult = br_a_2.signals.values(initial_latency+1:initial_latency + ticks_per_scale*(floor(length(br_a_2.signals.values)/ticks_per_scale)-1));


mat_split = reshape(vec_mult,ticks_per_scale,floor(length(br_a_2.signals.values)/ticks_per_scale) - 1);
lin_plot_mult = mean(abs(mat_split));

% Test quant output. This is quantized to 8 bits
% ----------------------------------------------
vec_quant =  q_10b_8b_2.signals.values(initial_latency+1:initial_latency + ticks_per_scale*(floor(length(q_10b_8b_2.signals.values)/ticks_per_scale)-1));
mat_split = reshape(vec_quant,ticks_per_scale,floor(length(q_10b_8b_2.signals.values)/ticks_per_scale) - 1);
lin_plot_quant = mean(abs(mat_split));



% Creat Linear line in DB
% -----------------------
div_vec = zeros(length(lin_plot_mult),1);
initial_vec = (repmat(1,length(lin_plot_mult),1));
for i=0:length(lin_plot_mult)-1
   div_vec(i+1) = 2^(i);    
end

linlog = initial_vec./div_vec;



% Plot results
% ------------
figure(1);
plot(20*log10(lin_plot_mult));
hold on;
plot(20*log10(lin_plot_quant),'r');
%plot(20*log10(linlog),'g');
xlabel('Scaled input iteration')
ylabel('dB')
title('Quantizer Linearity test')
%legend('Lin Mult','Lin Quant','Linear Ref')
legend('Lin Mult','Lin Quant')

diff = 20*log10(lin_plot_mult(1:plot_iterations)) - 20*log10(lin_plot_quant(1:plot_iterations));
lin_diff = zeros(length(diff)+1,1);
lin_diff(1) = nan;
lin_diff(2:end) = diff;
figure(2);
plot(lin_diff);
hold on;
xlabel('Scaled input iteration')
ylabel('dB')
title('Quantizer Linearity test: Difference')
legend('Lin Mult - Lin Quant')



% % ------------- %
% % 'Batman' Test %
% % ------------- %
% 
% % Output of CWG: 54 bits
% figure(2)
% hist(cwg_float_1_1.signals.values(initial_latency:end),hist_54b)
% title('Histogram: Output of CWG: 54bit')
% xlabel('Bins')
% ylabel('Bin hits')
% 
% % Output of CWG: 10 bits
% figure(3)
% hist(cwg_float_2_1.signals.values(initial_latency:end),hist_10b)
% title('Histogram Output of CWG: 10/18bit')
% xlabel('Bins')
% ylabel('Bin hits')
% 
% 
% % Input to Quant Multiplier: 54 bits
% figure(3)
% hist(q_a_54b_in_2.signals.values(initial_latency:end),hist_54b)
% title('Histogram: Input of Quant multiplier: 54bit')
% xlabel('Bins')
% ylabel('Bin hits')
% 
% % Input to Quant Multiplier: 10 bits
% figure(4)
% hist(q_a_10b_in_2.signals.values(initial_latency:end),hist_10b)
% title('Histogram: Input of Quant multiplier: 10/18bit')
% xlabel('Bins')
% ylabel('Bin hits')
% 
% 
% % Output of Quant Multiplier: 54 bits
% figure(5)
% hist(q_m_54b_2.signals.values(initial_latency:end),hist_54b)
% title('Histogram: Output of Quant multiplier: 54bit')
% xlabel('Bins')
% ylabel('Bin hits')
% 
% Output of Quant Multiplier: 10 bits
% figure(6)
% hist(q_m_10b_2.signals.values(initial_latency:end),hist_10b)
% title('Histogram: Output of Quant multiplier: 10/18bit')
% xlabel('Bins')
% ylabel('Bin hits')


% % Output of Quant: 54b in, 8 bits out
% figure(7)
% hist(q_54b_8b_2.signals.values(initial_latency:end),hist_8b)
% title('Histogram: Output of Quant: 54 bit in, 8bit out')
% xlabel('Bins')
% ylabel('Bin hits')

% Output of Quant: 10b in, 8 bits out
% figure(8)
% hist(q_10b_8b_2.signals.values(initial_latency:end),hist_8b)
% title('Histogram: Output of Quant: 10/18 bit in, 8bit out')
% xlabel('Bins')
% ylabel('Bin hits')


% Output of Quant: 10b in, 8 bits out
% figure(9)
% hist(br_a_2.signals.values(initial_latency:end),hist_8b)
% title('Histogram: Output of Quant: 10/18 bit in, 8bit out')
% xlabel('Bins')
% ylabel('Bin hits')

% input = br_a_2.signals.values(initial_latency:end); 
% out = input;
% for i=1:length(input)
%     if input(i) < (-1+1/(2^7))
%         out(i) = (-1+1/(2^7)); 
%     end
% end


% Test multiplier output in quant
% -------------------------------
initial_latency = 607;
ticks_per_scale = 4096*4;
hist_8b = 128-1;  % The -1 compensates for the bias correction as the values are not full scale, so 1 fewer bins needed.
hist_10b = 2^7;
hist_54b = 2^12;


% Input to Quant: 8 bits out
figure(10)
hist(br_a_2.signals.values(initial_latency:end),hist_8b)
title('Histogram: Input to Quant: 8bit out')
xlabel('Bins')
ylabel('Bin hits')

% Output of Quant: 8 bits out
figure(11)
hist(bc_out_2.signals.values(initial_latency+4:end),hist_8b)
title('Histogram: Output of Quant: 8bit out')
xlabel('Bins')
ylabel('Bin hits')


% Output of Quant: 10b in, 8 bits out
figure(12)
hist(mux_out_2.signals.values(initial_latency+4+2:end),hist_8b)
title('Histogram: Output of Quant: 10/18 bit in, 8bit out')
xlabel('Bins')
ylabel('Bin hits')
