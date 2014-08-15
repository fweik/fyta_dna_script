#
# (C) Copyright 2014, Damaris Holder, Florian Weik
#
# All rights reserved. 
#
# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#

source io.tcl
source analysis.tcl

# General MD parameters
set time_step 0.1
set total_int_steps 200000000
set steps_per_loop 10000
set skin 1.0

# Langevin parameters
# kT = 0.025eV ~ 300K (k_B = 8.617e-5 eV/K)
set kT 0.026
set gamma 1.0

# Output options
set vmd "no"
set vtf "yes"
set vtf_filename "/work/fweik/dna.vtf"

# Box geometry
# Length along the molecule
set box_z 2000
# Other directions
set box_xy 1000
# Shift along molecule axis
set zshift 8.5
# Shift in other directions
set center_xy [expr 0.5*$box_xy]

# Molecule
set n_basepairs 20
set configuration_filename "configurations/1000bp_conf.dat"
set sequence_filename "sequences/Sequence1.dat"
# Fix one end of molecule?
set fix_lower_end "yes"
# External forces on molecule
#Sheer force in +x direction
set ext_force_sheer 0.0
#Stretch force in +z direction
set ext_force_stretch 0.0

# Set up MD
setmd time_step $time_step
thermostat langevin $kT $gamma
setmd skin $skin
setmd box_l $box_xy $box_xy $box_z
cellsystem layered
set int_loops [expr $total_int_steps/$steps_per_loop]

# Read config
read_configuration $n_basepairs $configuration_filename

# Read sequence
set ladderlist [read_sequence $sequence_filename]

# Configure molecule

# setting charges
for { set i 0 } { $i <= [setmd max_part] } { incr i } {
      set type [part $i pr type]
      if { $type == 0 } {
	part $i q -1
      }
      if { $type == 2 } {
	part $i q -1
      }}

# masses

set massA 0.753
set massT 0.702
set massG 0.843
set massC 0.618

for { set i 0 } { $i <= [setmd max_part] } { incr i } {
  # AT
  # Sugar Left A
  if { $type == 1 } {
     if { [lindex $ladderlist [expr ($i+1)/4]] == 6 } {
     part $i mass $massA
     }
  }
  # Sugar Right T
  if { $type == 3 } {
     if { [lindex $ladderlist [expr ($i+1)/4]] == 6 } {
     part $i mass $massT
     }
  }
  # TA
  # Sugar Left T
  if { $type == 1 } {
     if { [lindex $ladderlist [expr ($i+1)/4]] == 6 } {
     part $i mass $massT
     }
  }
  # Sugar Right A
  if { $type == 3 } {
     if { [lindex $ladderlist [expr ($i+1)/4]] == 6 } {
     part $i mass $massA
     }
  }
  # GC
  # Sugar Left G
  if { $type == 1 } {
     if { [lindex $ladderlist [expr ($i+1)/4]] == 6 } {
     part $i mass $massG
     }
  }
  # Sugar Right C
  if { $type == 3 } {
     if { [lindex $ladderlist [expr ($i+1)/4]] == 6 } {
     part $i mass $massC
     }
  }
  # CG
  # Sugar Left C
  if { $type == 1 } {
     if { [lindex $ladderlist [expr ($i+1)/4]] == 6 } {
     part $i mass $massC
     }
  }
  # Sugar Right G
  if { $type == 3 } {
     if { [lindex $ladderlist [expr ($i+1)/4]] == 6 } {
     part $i mass $massG
     }
  }
}
      
# non-bonded IA between the sugars
inter forcecap 1
inter 0 0 lj-gen 1 1 2.5 0 4.976 -2 -4 18.773 -0.333
inter 2 2 lj-gen 1 1 2.5 0 4.976 -2 -4 18.773 -0.333
#inter 0 2 lj-gen 1 1 2.5 0 4.976 -2 -4 18.773 -0.333 

# electrostatic interactions

# set lB 0.00000000072
set lB 561
# 0.000000056
set kappa 9.6
set lambdaDB [expr 1./$kappa]
set oneoveralpha [expr 1/(-14.23)]
#set up of the electrostatic interactions
#inter coulomb lB dh κ rcut
# inter coulomb .78 dh 9.6 2
# inter coulomb .000000007 dh 9.6 2
#48
# inter coulomb bjerrum dh kappa r_cut eps_int eps r0 r1 1/alpha
inter coulomb $lB dh $lambdaDB [expr $lambdaDB*5] 78 3 4 13 $oneoveralpha
# inter coulomb $lB dh $lambdaDB [expr $lambdaDB*5]
# 300 K
# epsilon = 80
#dh:		debye hueckel specification
#kappa:		debye length
#r_cut:		cutoff radius (five times the debye length)
#eps_int:	dielectric constant in the interior of the DNA helix
#r0 and r1:	determine the boundaries between unscreened and screened electrostatic interactions
#alpha:		factor within exponential function to smooth the potential

set PI 3.14159

set phiA [expr ($PI/180)*54.53]
set phiT [expr ($PI/180)*55.93]
set phiG [expr ($PI/180)*52.69]
set phiC [expr ($PI/180)*54.87]
set sigmaA [expr ($PI/180)*18.67]
set sigmaT [expr ($PI/180)*16.82]
set sigmaG [expr ($PI/180)*25.57]
set sigmaC [expr ($PI/180)*22.43]
set r0AT 8.866
set r0GC 9.018
set r0sb 1.455
set E0AT -0.639
set E0GC -1.165
set E0sb -3.542
set kdAT 1.8
set kdGC 1.288
set lAT 0.703
set lGC 0.727
set alphaAT [expr 1/$lAT]
set alphaGC [expr 1/$lGC]
set alphasb [expr 1./0.4]
set f2 -0.132
set f3 +0.215

# adenyn and thymin #alpha=1/l
inter 211 cg_dna_basepair [list $r0AT $alphaAT $E0AT $kdAT $sigmaA $sigmaT $phiA $phiT $E0sb $r0sb $alphasb $f2 $f3 ]
# thymin and adenyn
inter 212 cg_dna_basepair [list $r0AT $alphaAT $E0AT $kdAT $sigmaT $sigmaA $phiT $phiA $E0sb $r0sb $alphasb $f2 $f3  ]
# guanin and cytosin
inter 221 cg_dna_basepair [list $r0GC $alphaGC $E0GC $kdGC $sigmaG $sigmaC $phiG $phiC $E0sb $r0sb $alphasb $f2 $f3  ]
# cytosin and guanin
inter 222 cg_dna_basepair [list $r0GC $alphaGC $E0GC $kdGC $sigmaC $sigmaG $phiC $phiG $E0sb $r0sb $alphasb $f2 $f3  ]

#TA-TA
inter 311 cg_dna_stacking {3.55 0.567 -0.492 0.028904 0.055176 0.058776 0.030718 0.059253 0.045492 0.027276 0. 0. 0. 0. 0. 0. 0.}
#AT-AT
inter 312 cg_dna_stacking {3.55 0.567 -0.492 0.028904 0.055176 0.058776 0.030718 0.059253 0.045492 0.027276 0. 0. 0. 0. 0. 0. 0.}
#CG-CG
inter 313 cg_dna_stacking {3.58 0.477 -0.49 0.088441 0.052882 0.026464 0.023318 0.054324 0.019645 0.021300 0. 0. 0. 0. 0. 0. 0.}
#GC-GC
inter 314 cg_dna_stacking {3.58 0.477 -0.49 0.088441 0.052882 0.026464 0.023318 0.054324 0.019645 0.021300 0. 0. 0. 0. 0. 0. 0.}
#AT-TA
inter 315 cg_dna_stacking {3.546 0.53 -0.494 -0.050096 0.026216 0.015713 -0.001645 -0.011908 -0.020130 -0.012014 0.047012 -0.027064 -0.006655 -0.014087 0.022053 0.023852 -0.008176}
#TA-AT
inter 316 cg_dna_stacking {3.67 0.513 -0.494 -0.050096 0.026216 0.015713 -0.001645 -0.011908 -0.020130 -0.012014 0.047012 -0.027064 -0.006655 -0.014087 0.022053 0.023852 -0.008176}
#GC-CG
inter 317 cg_dna_stacking {3.498 0.632 -0.492 -0.110166 -0.023904 0.002148 0.011226 -0.007213 0.002539 -0.009429 0.024314 -0.035707 0.000514 0.011428 0.038481 0.013361 0.013435}
#CG-GC
inter 318 cg_dna_stacking {3.537 0.563 -0.492 -0.110166 -0.023904 0.002148 0.011226 -0.007213 0.002539 -0.009429 0.024314 -0.035707 0.000514 0.011428 0.038481 0.013361 0.013435}
#GC-TA
inter 319 cg_dna_stacking {3.543 0.549 -0.501 -0.057278 0.025344 0.005640 0.010569 -0.018645 -0.003659 -0.008617 0.047173 -0.027305 -0.004606 -0.003417 0.031641 0.011550 0.004849}
#AT-CG
inter 320 cg_dna_stacking {3.543 0.549 -0.501 -0.057278 0.025344 0.005640 0.010569 -0.018645 -0.003659 -0.008617 0.047173 -0.027305 -0.004606 -0.003417 0.031641 0.011550 0.004849}
#TA-GC
inter 321 cg_dna_stacking {3.613 0.529 -0.50065 -0.057278 0.025344 0.005640 0.010569 -0.018645 -0.003659 -0.008617 -0.047173 0.027305 0.004606 0.003417 -0.031641 -0.011550 -0.004849}
#CG-AT
inter 322 cg_dna_stacking {3.613 0.529 -0.50065 -0.057278 0.025344 0.005640 0.010569 -0.018645 -0.003659 -0.008617 -0.047173 0.027305 0.004606 0.003417 -0.031641 -0.011550 -0.004849}
#GC-AT
inter 323 cg_dna_stacking {3.566 0.555 -0.501 0.015588 0.011820 0.022925 0.010727 0.044782 0.018715 0.010901 0.011257 -0.020028 0.003741 0.009697 -0.010209 -0.000221 -0.003437}
#TA-CG
inter 324 cg_dna_stacking {3.566 0.555 -0.501 0.015588 0.011820 0.022925 0.010727 0.044782 0.018715 0.010901 0.011257 -0.020028 0.003741 0.009697 -0.010209 -0.000221 -0.003437}
#AT-GC
inter 325 cg_dna_stacking {3.535 0.538 -0.500585 0.015588 0.011820 0.022925 0.010727 0.044782 0.018715 0.010901 -0.011257 0.020028 -0.003741 -0.009697 0.010209 0.000221 0.003437}
#CG-TA
inter 326 cg_dna_stacking {3.535 0.538 -0.500585 0.015588 0.011820 0.022925 0.010727 0.044782 0.018715 0.010901 -0.011257 0.020028 -0.003741 -0.009697 0.010209 0.000221 0.003437}

# backbone - backbone:
# Adenin
set kbss_A [expr 0.0001*(180/$PI)*(180/$PI)*3.489]
set theta_3_A [expr 94.17*$PI/180]
set theta_5_A [expr 61.36*$PI/180]
inter 111 angle_harmonic $kbss_A $theta_3_A
inter 112 angle_harmonic $kbss_A $theta_5_A

# Thymin
set kbss_T [expr 0.0001*(180/$PI)*(180/$PI)*4.689]
set theta_3_T [expr 92.67*$PI/180]
set theta_5_T [expr 68.40*$PI/180]
inter 121 angle_harmonic $kbss_T $theta_3_T
inter 122 angle_harmonic $kbss_T $theta_5_T

# Guanin
set kbss_G [expr 0.0001*(180/$PI)*(180/$PI)*4.165]
set theta_3_G [expr 90.30*$PI/180]
set theta_5_G [expr 63.69*$PI/180]
inter 131 angle_harmonic $kbss_G $theta_3_G
inter 132 angle_harmonic $kbss_G $theta_5_G

# Cytosin
set kbss_C [expr 0.0001*(180/$PI)*(180/$PI)*6.178]
set theta_3_C [expr 91.55*$PI/180]
set theta_5_C [expr 66.37*$PI/180]
inter 141 angle_harmonic $kbss_C $theta_3_C
inter 142 angle_harmonic $kbss_C $theta_5_C


for { set i 0 } { $i <= [setmd max_part] } { incr i } {
    set type [part $i pr type]
# AT
    # Sugar Left
    if { $type == 0 } {
    set s [lindex $ladderlist [expr ($i+5)/4-1]]
    set t [lindex $ladderlist [expr (($i+5)/4)]]
    puts "(i+5)/4: [expr ($i+5)/4]     ((i+5)/4)+1: [expr (($i+5)/4)+1]"
    puts "s: $s      t:$t"
    if { $s == 6 } {
    	# Basepair interaction
	# The order is s1 b1 b2 s2, different from the conf file
	part $i bond 211 [expr $i+1] [expr $i+3] [expr $i+2]
	# Check if first bp
	if { $i > 0 } {
	    # Angular potential with previous sugar and base
	    part $i bond 112 [expr $i+1] [expr $i - 4]
	}
	# Check if last bp
	if { $i <= [expr [setmd max_part] - 4] } {
	    # Angular potential with next sugar and base
	    part $i bond 111 [expr $i+4] [expr $i+1]
	    # Stacking interaction
	    # The order is s1 b1 b2 s2 s1 b1 b2 s2, different from the conf file
	    if { $t == 6 } {
		part $i bond 312 [expr $i+1] [expr $i+3] [expr $i+2] [expr $i+4] [expr $i+5] [expr $i+7] [expr $i+6]
	    }
	    if { $t == 7 } {
		part $i bond 315 [expr $i+1] [expr $i+3] [expr $i+2] [expr $i+4] [expr $i+5] [expr $i+7] [expr $i+6]
	    }
    	    if { $t == 8 } {
		part $i bond 325 [expr $i+1] [expr $i+3] [expr $i+2] [expr $i+4] [expr $i+5] [expr $i+7] [expr $i+6]
	    }
    	    if { $t == 9 } {
		part $i bond 320 [expr $i+1] [expr $i+3] [expr $i+2] [expr $i+4] [expr $i+5] [expr $i+7] [expr $i+6]
	    }
	}
    }}
    # Sugar Right
    if { $type == 2 } {
    set s [lindex $ladderlist [expr ($i+5)/4]]
    if { $s == 6 } {
	if { $i > 3 } {
	    # Angular potential with previous sugar and base
	    part $i bond 121 [expr $i+1] [expr $i - 4]
	}
	# Check if last bp
	if { $i <= [expr [setmd max_part] - 4] } {
	    # Angular potential with next sugar and base
	    part $i bond 122 [expr $i+4] [expr $i+1]
	}
    }}

# TA
    # Sugar Left
    if { $type == 0 } {
    set s [lindex $ladderlist [expr ($i+5)/4-1]]
    set t [lindex $ladderlist [expr (($i+5)/4)]]
    puts "(i+5)/4: [expr ($i+5)/4]     ((i+5)/4)+1: [expr (($i+5)/4)+1]"
    puts "s: $s      t:$t"
    if { $s == 7 } {
    	# Basepair interaction
	# The order is s1 b1 b2 s2, different from the conf file
	part $i bond 212 [expr $i+1] [expr $i+3] [expr $i+2]
	# Check if first bp
	if { $i > 0 } {
	    # Angular potential with previous sugar and base
	    part $i bond 122 [expr $i+1] [expr $i - 4]
	}
	# Check if last bp
	if { $i <= [expr [setmd max_part] - 4] } {
	    # Angular potential with next sugar and base
	    part $i bond 121 [expr $i+4] [expr $i+1]
	    # Stacking interaction
	    # The order is s1 b1 b2 s2 s1 b1 b2 s2, different from the conf file
	    if { $t == 6 } {
		part $i bond 316 [expr $i+1] [expr $i+3] [expr $i+2] [expr $i+4] [expr $i+5] [expr $i+7] [expr $i+6]
	    }
	    if { $t == 7 } {
		part $i bond 311 [expr $i+1] [expr $i+3] [expr $i+2] [expr $i+4] [expr $i+5] [expr $i+7] [expr $i+6]
	    }
    	    if { $t == 8 } {
		part $i bond 321 [expr $i+1] [expr $i+3] [expr $i+2] [expr $i+4] [expr $i+5] [expr $i+7] [expr $i+6]
	    }
    	    if { $t == 9 } {
		part $i bond 324 [expr $i+1] [expr $i+3] [expr $i+2] [expr $i+4] [expr $i+5] [expr $i+7] [expr $i+6]
	    }
	}
    }}
    # Sugar Right
    if { $type == 2 } {
    set s [lindex $ladderlist [expr ($i+5)/4]]
    if { $s == 7 } {
	if { $i > 3 } {
	    # Angular potential with previous sugar and base
	    part $i bond 111 [expr $i+1] [expr $i - 4]
	}
	# Check if last bp
	if { $i <= [expr [setmd max_part] - 4] } {
	    # Angular potential with next sugar and base
	    part $i bond 112 [expr $i+4] [expr $i+1]
	}
    }}

# GC
    # Sugar Left
    if { $type == 0 } {
    set s [lindex $ladderlist [expr ($i+5)/4-1]]
    set t [lindex $ladderlist [expr (($i+5)/4)]]
    puts "(i+5)/4: [expr ($i+5)/4]     ((i+5)/4)+1: [expr (($i+5)/4)+1]"
    puts "s: $s      t:$t"
    if { $s == 8 } {
    	# Basepair interaction
	# The order is s1 b1 b2 s2, different from the conf file
	part $i bond 221 [expr $i+1] [expr $i+3] [expr $i+2]
	# Check if first bp
	if { $i > 0 } {
	    # Angular potential with previous sugar and base
	    part $i bond 132 [expr $i+1] [expr $i - 4]
	}
	# Check if last bp
	if { $i <= [expr [setmd max_part] - 4] } {
	    # Angular potential with next sugar and base
	    part $i bond 131 [expr $i+4] [expr $i+1]
	    # Stacking interaction
	    # The order is s1 b1 b2 s2 s1 b1 b2 s2, different from the conf file
	    if { $t == 6 } {
		part $i bond 323 [expr $i+1] [expr $i+3] [expr $i+2] [expr $i+4] [expr $i+5] [expr $i+7] [expr $i+6]
	    }
	    if { $t == 7 } {
		part $i bond 319 [expr $i+1] [expr $i+3] [expr $i+2] [expr $i+4] [expr $i+5] [expr $i+7] [expr $i+6]
	    }
    	    if { $t == 8 } {
		part $i bond 314 [expr $i+1] [expr $i+3] [expr $i+2] [expr $i+4] [expr $i+5] [expr $i+7] [expr $i+6]
	    }
    	    if { $t == 9 } {
		part $i bond 317 [expr $i+1] [expr $i+3] [expr $i+2] [expr $i+4] [expr $i+5] [expr $i+7] [expr $i+6]
	    }
	}
    }}
    # Sugar Right
    if { $type == 2 } {
    set s [lindex $ladderlist [expr ($i+5)/4]]
    if { $s == 8 } {
	if { $i > 3 } {
	    # Angular potential with previous sugar and base
	    part $i bond 141 [expr $i+1] [expr $i - 4]
	}
	# Check if last bp
	if { $i <= [expr [setmd max_part] - 4] } {
	    # Angular potential with next sugar and base
	    part $i bond 142 [expr $i+4] [expr $i+1]
	}
    }}

# CG
    # Sugar Left
    if { $type == 0 } {
    set s [lindex $ladderlist [expr ($i+5)/4-1]]
    set t [lindex $ladderlist [expr (($i+5)/4)]]
    puts "(i+5)/4: [expr ($i+5)/4]     ((i+5)/4)+1: [expr (($i+5)/4)+1]"
    puts "s: $s      t:$t"
    if { $s == 9 } {
    	# Basepair interaction
	# The order is s1 b1 b2 s2, different from the conf file
	part $i bond 222 [expr $i+1] [expr $i+3] [expr $i+2]
	# Check if first bp
	if { $i > 0 } {
	    # Angular potential with previous sugar and base
	    part $i bond 142 [expr $i+1] [expr $i - 4]
	}
	# Check if last bp
	if { $i <= [expr [setmd max_part] - 4] } {
	    # Angular potential with next sugar and base
	    part $i bond 141 [expr $i+4] [expr $i+1]
	    # Stacking interaction
	    # The order is s1 b1 b2 s2 s1 b1 b2 s2, different from the conf file
	    if { $t == 6 } {
		part $i bond 322 [expr $i+1] [expr $i+3] [expr $i+2] [expr $i+4] [expr $i+5] [expr $i+7] [expr $i+6]
	    }
	    if { $t == 7 } {
		part $i bond 326 [expr $i+1] [expr $i+3] [expr $i+2] [expr $i+4] [expr $i+5] [expr $i+7] [expr $i+6]
	    }
    	    if { $t == 8 } {
		part $i bond 318 [expr $i+1] [expr $i+3] [expr $i+2] [expr $i+4] [expr $i+5] [expr $i+7] [expr $i+6]
	    }
    	    if { $t == 9 } {
		part $i bond 313 [expr $i+1] [expr $i+3] [expr $i+2] [expr $i+4] [expr $i+5] [expr $i+7] [expr $i+6]
	    }
	}
    }}
    # Sugar Right
    if { $type == 2 } {
    set s [lindex $ladderlist [expr ($i+5)/4]]
    if { $s == 9 } {
	if { $i > 3 } {
	    # Angular potential with previous sugar and base
	    part $i bond 131 [expr $i+1] [expr $i - 4]
	}
	# Check if last bp
	if { $i <= [expr [setmd max_part] - 4] } {
	    # Angular potential with next sugar and base
	    part $i bond 132 [expr $i+4] [expr $i+1]
	}
    }}
    
}

if { $fix_lower_end == "yes" } {
    part 0 fix
    part 2 fix

    puts [part 0]
    puts [part 2]
}

analyze_bps

for { set i 0 } { $i <= [setmd max_part] } { incr i } {
    puts [part $i]
}

# Prepare output

if { $vmd == "yes" } {
  prepare_vmd_connection "dna_bla"
  after 3000
  imd positions
}

if { $vtf == "yes" } {
    set vtfchan [open $vtf_filename w]

    writevsf $vtfchan
    writevcf $vtfchan
}

puts "Integrating $int_loops times $steps_per_loop steps."

set fo [open "averagesPOLYseq.dat" "w"]
set pers [open "persistencePOLYseq.dat" "w"]
# set forces [open "ForcesPOLYAT.dat" "w"]


set largest 0

for { set i 1 } { $i <= $int_loops } { incr i } {
    puts "Loop $i of $int_loops."
    puts [time { integrate $steps_per_loop }]

#     for { set k 1 } { $k <= [setmd max_part] } { incr k } {
#     puts $forces [part $k print f]
#     }

#     for { set k 1 } { $k <= [setmd max_part] } { incr k } {
#     set x [part $k print f]
#     set F [expr [lindex $x 0]*[lindex $x 0]+[lindex $x 1]*[lindex $x 1]+[lindex $x 2]*[lindex $x 2]]
#     if { $F > $largest} {
#     set largest $F
#     }
#     }
#     puts $forces "$largest"
    #puts $fo [analyze_bps]
    #puts $pers [analyze_pl]

    if { $vmd == "yes" } { 
#	after 10
	imd positions
    }

    if { $vtf == "yes" } {
	writevcf $vtfchan
    }
}

close $fo
