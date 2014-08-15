setmd time_step 0.1
# kT = 0.025eV ~ 300K (k_B = 8.617e-5 eV/K)
# thermostat langevin 0.025 0.0000001
# thermostat langevin 0.025 0.0001 - 11
# thermostat langevin 0.025 0.001  - 13
# thermostat langevin 0.025 0.01   - 18
# thermostat langevin 0.025 0.1    - 100
thermostat langevin 0.026 1
setmd skin 0.5
set vmd "no"

set vtf_filename "/work/fweik/dna.vtf"
set vtf "yes"

# Along the molecule
set box_z 2000
# Other directions
set box_xy 1000
set zshift 8.5
set center_xy [expr 0.5*$box_xy]

set n_basepairs 20

set total_int_steps 200000000
set steps_per_loop 10000

set fix_lower_end "yes"

#Sheer force in +x direction
set ext_force_sheer 0.0

#Stretch force in +z direction
set ext_force_stretch 0.0

set sequence [open Sequence1.dat]

setmd box_l $box_xy $box_xy $box_z

cellsystem layered

 set filename 1000bp_config.dat
#set filename 200bp_config.dat
# set filename 200bp_config.dat
# set filename 16bp_config.dat
# set filename 8bp_config.dat

read_configuration $n_basepairs $filename

if { $fix_lower_end == "yes" } {
    part 0 fix
    part 2 fix

    puts [part 0]
    puts [part 2]
}

set ladderlist []
while {1} {
    set line [gets $sequence]
    if {[eof $sequence]} {
        close $sequence
        break
    }
    lappend ladderlist $line
}

set int_loops [expr $total_int_steps/$steps_per_loop]

proc read_configuration { n_basepairs filename } {
    set v1 0
    set v2 0
    set v3 0
    set v4 0

    global center_xy zshift

    set n_particles [expr 4*$n_basepairs]

    set chan [open $filename r]

    # Read in configuration

    set lastid 0
    while {([gets $chan line] >= 0) && ($lastid < $n_particles)} {
	set vcnt [scan $line %s%s%s%s v1 v2 v3 v4]

	if { $v1 == 0 } {	        
	    # center particles in x and y
	    set type [expr $lastid % 4]
	    if { $type == 0 || $type == 2 } {
		set q -1
	    } else {
		set q 0
	    }
	      
	    part $lastid pos [expr $center_xy + $v2] [expr $center_xy + $v3] [expr $v4 + $zshift] type $type q $q
	    incr lastid
	}
    }
    if { $lastid != $n_particles } {
	puts "Not enough particles in file '$filename' to read [expr $n_basepairs/4] basepairs."
	exit
    }   
}


# Analysis functions

proc vecdot {v1 v2} {
    if {[llength $v1]!=[llength $v2]} {
    error "Vectors must be of equal length"
    }
   set dot [expr [lindex $v1 0]*[lindex $v2 0]+[lindex $v1 1]*[lindex $v2 1]+[lindex $v1 2]*[lindex $v2 2]]
   return $dot
}

proc analyze_pl { } {
    set Pi 3.14159
    set data1 [list]
    set data2 [list]
    
  for {set j 0} { $j < [expr ([setmd max_part] - 3)/4] } {incr j} {
    set m [expr $j*4]
    set s1 $m
    
    set b1 [expr $s1 + 1]
    set b2 [expr $s1 + 3]
    set b3 [expr $s1 + 5]
    set b4 [expr $s1 + 7]
    
    set x1 [part $b1 pr pos]
    set x2 [part $b2 pr pos]
    set x3 [part $b3 pr pos]
    set x4 [part $b4 pr pos]
    
    set m12x [expr (([lindex $x2 0] + [lindex $x1 0])/2)]
    set m12y [expr (([lindex $x2 1] + [lindex $x1 1])/2)]
    set m12z [expr (([lindex $x2 2] + [lindex $x1 2])/2)]
    set m34x [expr (([lindex $x4 0] + [lindex $x3 0])/2)]
    set m34y [expr (([lindex $x4 1] + [lindex $x3 1])/2)]
    set m34z [expr (([lindex $x4 2] + [lindex $x3 2])/2)]
    
    set v1x [expr ($m34x-$m12x)]
    set v1y [expr ($m34y-$m12y)]
    set v1z [expr ($m34z-$m12z)]
    set v1 [list $v1x $v1y $v1z]
    set lent1 [expr {sqrt($v1x*$v1x+$v1y*$v1y+$v1z*$v1z)}]
    
#     set cosang [list]
#     set actuallengths [list]
    
    # one loop in order to avoid double counting
    # loop: from base s1 to last base
    set len [list $lent1]
    set cosang [list]
    set actuallengths [list]
    
#     [expr ([setmd max_part] - $s1)/4 - 2]
    
    for { set i 0 } { $i <= [expr ([setmd max_part] - $s1)/4 - 1] } { incr i } {
      set r [expr $i*4 + 1]
      set k [expr $i*4 + 3]
      
      set b3 [expr $s1 + $r]
      set b4 [expr $s1 + $k]
      set b5 [expr $s1 + $r + 4]
      set b6 [expr $s1 + $k + 4]
      
      set x3 [part $b3 pr pos]
      set x4 [part $b4 pr pos]
      set x5 [part $b5 pr pos]
      set x6 [part $b6 pr pos]
      
      set m34x [expr (([lindex $x4 0] + [lindex $x3 0])/2)]
      set m34y [expr (([lindex $x4 1] + [lindex $x3 1])/2)]
      set m34z [expr (([lindex $x4 2] + [lindex $x3 2])/2)]
      set m56x [expr (([lindex $x6 0] + [lindex $x5 0])/2)]
      set m56y [expr (([lindex $x6 1] + [lindex $x5 1])/2)]
      set m56z [expr (([lindex $x6 2] + [lindex $x5 2])/2)]

      set v2x [expr ($m56x-$m34x)]
      set v2y [expr ($m56y-$m34y)]
      set v2z [expr ($m56z-$m34z)]
      set v2 [list $v2x $v2y $v2z]
      set lent2 [expr {sqrt($v2x*$v2x+$v2y*$v2y+$v2z*$v2z)}]
#       puts "number $i: $lent2"
      
#       lappend len $lent2
      lappend len [expr [lindex $len end] + $lent2]
      set actuallengths [lreplace $len 0 0]
      
      set dotprod [vecdot $v1 $v2]
      lappend cosang [expr $dotprod/($lent1 * $lent2)]
      
    }
#     puts "S1: $s1"
#     puts "ANGLES: $cosang"
#     puts "LENGTHS: $actuallengths"
  lappend data1 [concat $cosang]
  lappend data2 [concat $actuallengths]
#      puts "DATA1: $data1"
#      puts "DATA1: [concat $data1]"
#      puts "DATA2: $data2"

#   puts "ang $j: $actuallengths"
#   puts "data2: $data2"
    
  }
  
  set arclengths [list]
  for { set k 0 } { $k <= [expr ([setmd max_part]/4) - 2] } { incr k } {
    set value 0.0
    for { set r 0 } { $r <= [expr ([setmd max_part]/4) - 2 - $k] } { incr r } {
#       puts "r: $r"
#     puts "lindex $k: [lindex $data2 $r $k]"
    set value [expr $value + [lindex $data2 $r $k]]
    }
    set value [expr $value/((([setmd max_part]+1)/4) - 2 - $k)]
    lappend arclengths $value
#     puts "VALUE: $value"
  }

  set angles [list]
  for { set k 0 } { $k <= [expr ([setmd max_part]/4) - 2] } { incr k } {
    set value 0.0
    for { set r 0 } { $r <= [expr ([setmd max_part]/4) - 2 - $k] } { incr r } {
#     puts "lindex $k: [lindex $data1 $r $k]"
    set value [expr $value + [lindex $data1 $r $k]]
    }
    set value [expr $value/((([setmd max_part]+1)/4) - 2 - $k)]
    lappend angles $value
#     puts "VALUE: $value"
  }
  
  #number of evaluated basepairs:
  set nbp [expr ([setmd max_part]+1)/4-3]
  set arclengths [lreplace $arclengths $nbp $nbp]
  set angles [lreplace $angles $nbp $nbp]
  
#   return [concat "[lindex $data2 0] [lindex $data1 0]"]
  return [concat "$arclengths $angles"]
  }

proc analyze_bp { s1 } {
    set Pi 3.14159

    set bond [lindex [part $s1 pr bond] 0]

    set s2 [expr $s1 + 2]
    set b1 [expr $s1 + 1]
    set b2 [expr $s1 + 3]
    
    set x1 [part $s1 pr pos]
    set x2 [part $s2 pr pos]

    set dx [expr ([lindex $x2 0] - [lindex $x1 0])]
    set dy [expr ([lindex $x2 1] - [lindex $x1 1])]
    set dz [expr ([lindex $x2 2] - [lindex $x1 2])]
    set rcc [expr sqrt($dx*$dx + $dy*$dy + $dz*$dz)]

    set x2 [part $b1 pr pos]
    set dxcb [expr ([lindex $x2 0] - [lindex $x1 0])]
    set dycb [expr ([lindex $x2 1] - [lindex $x1 1])]
    set dzcb [expr ([lindex $x2 2] - [lindex $x1 2])]
    set rcb1 [expr sqrt($dxcb*$dxcb + $dycb*$dycb + $dzcb*$dzcb)]

    set psi1 [expr ($dx*$dxcb + $dy*$dycb + $dz*$dzcb)/($rcb1*$rcc)]
    set psi1 [expr 180. * acos($psi1)/$Pi]

    set x1 [part $s2 pr pos]
    set x2 [part $b2 pr pos]

    set dxcb [expr ([lindex $x2 0] - [lindex $x1 0])]
    set dycb [expr ([lindex $x2 1] - [lindex $x1 1])]
    set dzcb [expr ([lindex $x2 2] - [lindex $x1 2])]
    set rcb2 [expr sqrt($dxcb*$dxcb + $dycb*$dycb + $dzcb*$dzcb)]
    
    set psi2 [expr -($dx*$dxcb + $dy*$dycb + $dz*$dzcb)/($rcb2*$rcc)]
    set psi2 [expr 180. * acos($psi2)/$Pi]

    set x1 [part $b2 pr pos]
    set x2 [part $b1 pr pos]
    set dx [expr abs(([lindex $x2 0] - [lindex $x1 0]))]
    set dy [expr abs(([lindex $x2 1] - [lindex $x1 1]))]
    set dz [expr abs(([lindex $x2 2] - [lindex $x1 2]))]
    set rhb [expr sqrt($dx*$dx + $dy*$dy + $dz*$dz)]

    return  [list $rcc $rcb1 $rcb2 $rhb $psi1 $psi2]
}

proc analyze_bps {} {
    set rcc_avg 0.0
    set rcb1_avg 0.0
    set rcb2_avg 0.0
    set rhb_avg 0.0
    set psi1_avg 0.0
    set psi2_avg 0.0
    set n 0

    for { set i 0 } { $i <= [setmd max_part] } { incr i } {
	set type [part $i pr type]
	if { $type == 0 } {
	    set bp [analyze_bp $i]
	    set rcc_avg [expr $rcc_avg + [lindex $bp 0]]
	    set rcb1_avg [expr $rcb1_avg + [lindex $bp 1]]
	    set rcb2_avg [expr $rcb2_avg + [lindex $bp 2]]
	    set rhb_avg [expr $rhb_avg + [lindex $bp 3]]
	    set psi1_avg [expr $psi1_avg + [lindex $bp 4]]
	    set psi2_avg [expr $psi2_avg + [lindex $bp 5]]

	    incr n
	}
    }

    set rcc_avg [expr $rcc_avg / $n]
    set rcb1_avg [expr $rcb1_avg / $n]
    set rcb2_avg [expr $rcb2_avg / $n]
    set rhb_avg [expr $rhb_avg / $n]
    set psi1_avg [expr $psi1_avg / $n]
    set psi2_avg [expr $psi2_avg / $n]

    puts "<rcc> $rcc_avg, <rcb1> $rcb1_avg, <rcb2> $rcb2_avg, <rhb> $rhb_avg, <psi1> $psi1_avg, <psi2> $psi2_avg"
    return [list $rcc_avg $rcb1_avg $rcb2_avg $rhb_avg $psi1_avg $psi2_avg]
}

# Read in configuration

set lastid 0
while {[gets $chan line] >= 0} {
    set vcnt [scan $line %s%s%s%s v1 v2 v3 v4]

    if { $v1 == 0 } {	        
	# center particles in x and y
	part $lastid pos [expr $center_xy + $v2] [expr $center_xy + $v3] [expr $v4 + $zshift] type [expr $lastid % 4]
	puts [part $lastid pr id pos]
	incr lastid
    }
}

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

integrate 0

analyze_bps

for { set i 0 } { $i <= [setmd max_part] } { incr i } {
    puts [part $i]
}

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

#puts "Center of Mass is [analyze centermass 0]"

puts "Force calculation takes [time { integrate 0 }]"

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
