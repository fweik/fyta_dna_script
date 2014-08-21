

set xlabel "contour length [Ang]"
set ylabel "<R^2> [Ang^2]"
set key top left

g(x)=2*P*x*(1.-(P/x)*(1. - exp(-x/P)))

wlc(x,P)=2*P*x*(1.-(P/x)*(1. - exp(-x/P)))

f(x)=a*x+b

fit g(x) "avg_ete_lB_0.dat" u 4:2:3 via P
P_no_es=P
fit g(x) "avg_ete_lB_561.dat" u 4:2:3 via P
P_100mM=P
fit g(x) "avg_ete_lB_.dat" u 4:2:3 via P
P_10mM=P
#fit [800:] f(x) "avg_ete_lB_0.dat" u 4:2:3 via a,b

plot [] [-100000:] 'avg_ete_lB_0.dat' u 4:2:5:3 w xyerror t 'simulation data no electrostatics', wlc(x, P_no_es) t sprintf("worm-like chain for P=%.1f nm", 0.1*P_no_es), \
             'avg_ete_lB_561.dat' u 4:2:5:3 w xyerror t 'simulation data 0.1M', wlc(x, P_100mM) t sprintf("worm-like chain for P=%.1f nm", 0.1*P_100mM), \
             'avg_ete_lB_561.dat' u 4:2:5:3 w xyerror t 'simulation data 0.01M', wlc(x, P_10mM) t sprintf("worm-like chain for P=%.1f nm", 0.1*P_10mM)

set term postscript eps enhanced color "Helvetica" 12
set output 'fit_wormlike_chain.eps'

plot [] [0:] 'avg_ete_lB_0.dat' u 4:2:5:3 w xyerror t 'simulation data no electrostatics', wlc(x, P_no_es) t sprintf("worm-like chain for P=%.1f nm", 0.1*P_no_es), \
             'avg_ete_lB_561.dat' u 4:2:5:3 w yerror t 'simulation data 0.1M', wlc(x, P_100mM) t sprintf("worm-like chain for P=%.1f nm", 0.1*P_100mM)

!epstopdf 'fit_wormlike_chain.eps'

pause -1
