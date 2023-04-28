figure, loglog(FP_x,PYm_x,'b','linewidth',2), hpsx=gcf;
xlabel('f (Hz)'), ylabel('Power Spectrum in X (V^2 s)')
[0.097 0.17 0.27 0.43 0.57 0.91 1.23 2.11 3.03 5.07 8.77 19.9 36.7 83 170]


xline(0.097,':k'); text(0.097,0.5,'0.097','HorizontalAlignment','center') % 1 y
xline(0.17,':k'); text(0.17,0.5,'0.097','HorizontalAlignment','center')
xline(0.27,':k'); text(0.27,0.5,'0.097','HorizontalAlignment','center')
xline(0.43,':k'); text(0.097,0.5,'0.43','HorizontalAlignment','center') % 1 y
xline(0.57,':k'); text(0.097,0.5,'0.57','HorizontalAlignment','center')
xline(0.91,':k'); text(0.097,0.5,'0.91','HorizontalAlignment','center')
xline(1.23,':k'); text(0.097,0.5,'0.097','HorizontalAlignment','center') % 1 y
xline(2.11,':k'); text(0.097,0.5,'0.097','HorizontalAlignment','center')
xline(3.03,':k'); text(0.097,0.5,'0.097','HorizontalAlignment','center')
xline(5.07,':k'); text(0.097,0.5,'0.097','HorizontalAlignment','center') % 1 y
xline(8.77,':k'); text(0.097,0.5,'0.097','HorizontalAlignment','center')
xline(19.9,':k'); text(0.097,0.5,'0.097','HorizontalAlignment','center')
xline(36.7,':k'); text(0.097,0.5,'0.097','HorizontalAlignment','center') % 1 y
xline(83,':k'); text(0.097,0.5,'0.097','HorizontalAlignment','center')
xline(163,':k'); text(0.097,0.5,'0.097','HorizontalAlignment','center')