

set xlabel "contour length [Ang]"
set ylabel "<R^2> [Ang^2]"
set key top left

g(x)=2*P*x*(1.-(P/x)*(1. - exp(-x/P)))

wlc(x,P)=2*P*x*(1.-(P/x)*(1. - exp(-x/P)))

f(x)=a*x+b

fit [:400] g(x) "avg_ete_lambda_9.6.dat" u 4:2 via P
P_100mM=P
fit [:400] g(x) "avg_ete_lambda_30.dat" u 4:2 via P
P_10mM=P

plot [] [] 'avg_ete_lambda_9.6.dat' u 4:2:5:3 w xyerror t 'simulation data 0.1M', wlc(x, P_100mM) t sprintf("worm-like chain for P=%.1f nm", P_100mM), \
     'avg_ete_lambda_30.dat' u 4:2:5:3 w xyerror t 'simulation data 0.01M', wlc(x, P_10mM) t sprintf("worm-like chain for P=%.1f nm", P_10mM)

set term postscript eps enhanced color "Helvetica" 12
set output 'fit_wormlike_chain.eps'

plot [] [] 'avg_ete_lambda_9.6.dat' u 4:2:5:3 w xyerror t 'simulation data 0.1M', wlc(x, P_100mM) t sprintf("worm-like chain for P=%.1f nm", P_100mM), \
     'avg_ete_lambda_30.dat' u 4:2:5:3 w xyerror t 'simulation data 0.01M', wlc(x, P_10mM) t sprintf("worm-like chain for P=%.1f nm", P_10mM)

!epstopdf 'fit_wormlike_chain.eps'

pause -1
