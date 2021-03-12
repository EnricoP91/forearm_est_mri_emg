function update_plot(figure_struct, time_emg, win_start, semg_selected, win_time, data_est )

set(figure_struct.rsub,'Position', [time_emg(win_start), min(semg_selected(1,:)), win_time, figure_struct.heightWin])
figure_struct.ax3.XLim = [time_emg(win_start)-1, time_emg(win_start)+win_time+1];
set(figure_struct.bsub, 'YData', data_est.Iest)
set(figure_struct.psub, 'RData',data_est.Voutarm_f) 
set(figure_struct.psubreal,'RData',data_est.rms_semgreal_f)
drawnow