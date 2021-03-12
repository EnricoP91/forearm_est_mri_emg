function output = prepare_figure(semg_selected, dimensions, muscles, emg_time, win_time)
%% ------------------------------------------------------------------------
%     Creation and initialization of the plot. The plot is composed of 3
%     subplots:
%     i) polar plot representing the profile of EMG around the forearm
%     ii) bar plot representing the instantaneous activation of muscles
%     iii) time evolution plot representing the portion of EMG analyzed on
%     a random selected channel
% -------------------------------------------------------------------------
 if numel(dimensions) > 2
    error("The dimension value in input are more    than 2.")
 else
    if mod(dimensions(1),1)==0 && mod(dimensions(2),1) == 0
        n_electrode = dimensions(1);
        n_muscle = dimensions(2);
    end
 end
 
 % initialization 
 output.heightWin = abs(max(semg_selected(1,:)))+abs(min(semg_selected(1,:)));
 theta = linspace(0, 2*pi, n_electrode+1);
 
 % figure
 output.h1 = figure('Units','normalized','OuterPosition',[0 0 1 1]);
 
 % polarplot 
 output.ax1 = subplot(2,2,1);
 output.psubreal = polarplot(theta,zeros(size(theta)), 'ro--', 'Linewidth',2);
 hold on
 output.psub = polarplot(theta,zeros(size(theta)),'bo-');
 hold off
 legend('Measured EMG', 'Estimated Reconstruction')
 title("Voltage Profile [mV]")
 rlim([0, 300])
 
 % bar plot    
 output.ax2 = subplot(2,2,2);
 output.bsub = bar(muscles,zeros(1,n_muscle));
 ylabel("Current [A]")
 title("Muscles current estimation")
 ylim([-1e-3,3e-3])
 grid on
 
 % time evolution     
 output.ax3 = subplot(2,2,[3 4]);
 hold on
 plot(emg_time, semg_selected(1,:),'Linewidth',0.5)
 grid on
 xlim([win_time-1,win_time+1])
 ylim([min(semg_selected(1,:)), max(semg_selected(1,:))])

 output.rsub = rectangle('Position',[0, min(semg_selected(1,:)), win_time, output.heightWin], 'EdgeColor', 'r',...
                        'Linewidth',2);
 xlabel("Time [s]")
 ylabel("Voltage [mV]")
 title("EMG Profile (one channel)")
                    