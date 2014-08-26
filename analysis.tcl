#
# (C) Copyright 2014, Damaris Holder, Florian Weik
#
# All rights reserved. 
#
# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#

proc vecdot {v1 v2} {
    if {[llength $v1]!=[llength $v2]} {
    error "Vectors must be of equal length"
    }
   set dot [expr [lindex $v1 0]*[lindex $v2 0]+[lindex $v1 1]*[lindex $v2 1]+[lindex $v1 2]*[lindex $v2 2]]
   return $dot
}

proc vecadd {c1 v1 c2 v2} {
    if {[llength $v1]!=[llength $v2]} {
	error "Vectors must be of equal length"
    }
    for { set i 0 } { $i < [llength $v1] } { incr i } {
	lappend ret [expr $c1*[lindex $v1 $i] + $c2*[lindex $v2 $i]]
    }
    return $ret
}

proc veclen {v1} {
    return [expr { sqrt( [vecdot $v1 $v1] ) }]
}

proc veccross {v1 v2} {
    lappend ret [expr [lindex $v1 1]*[lindex $v2 2] - [lindex $v1 2]*[lindex $v2 1]]
    lappend ret [expr [lindex $v1 2]*[lindex $v2 0] - [lindex $v1 0]*[lindex $v2 2]]
    lappend ret [expr [lindex $v1 0]*[lindex $v2 1] - [lindex $v1 1]*[lindex $v2 0]]
    
    return $ret
}

proc analyze_end_to_end_sq { } {
    set s1b 0
    set s2b 2
    
    set s1e [expr [setmd max_part] - 3]
    set s2e [expr $s1e + 2]

    set r1b [part $s1b pr pos]
    set r2b [part $s2b pr pos]
    set r1e [part $s1e pr pos]
    set r2e [part $s2e pr pos]        

    set r1 [vecadd 0.5 $r1b 0.5 $r2b]
    set r2 [vecadd 0.5 $r1e 0.5 $r2e]

    set R [vecadd 1.0 $r1 -1.0 $r2]

    return [vecdot $R $R]
}



proc analyze_contour_length {} {
    set l 0.0
    for { set i 0 } { $i < [expr [setmd max_part]/4 - 1] } { incr i } {
	set s1i [expr 4*$i]
	set s2i [expr $s1i + 2]
	set s1j [expr $s1i + 4]
	set s2j [expr $s1i + 6]

	set r1i [part $s1i pr pos]
	set r2i [part $s2i pr pos]
	set r1j [part $s1j pr pos]
	set r2j [part $s2j pr pos]        

	set r1 [vecadd 0.5 $r1i 0.5 $r2i]
	set r2 [vecadd 0.5 $r1j 0.5 $r2j]

	set R [vecadd 1.0 $r1 -1.0 $r2]

	set l [expr $l + [veclen $R]]	
    }
    return $l
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

    puts "<rcc> [format %.2f $rcc_avg], <rcb1> [format %.2f $rcb1_avg], <rcb2> [format %.2f $rcb2_avg], <rhb> [format %.2f $rhb_avg], <psi1> [format %.2f $psi1_avg], <psi2> [format %.2f $psi2_avg]"
    return [list $rcc_avg $rcb1_avg $rcb2_avg $rhb_avg $psi1_avg $psi2_avg]
}

proc analyze_stacking { s1i } {
    set Pi 3.14159

    set s2i [expr $s1i + 2]
    set b1i [expr $s1i + 1]
    set b2i [expr $s1i + 3]
    set s1j [expr $s1i + 4]
    set b1j [expr $s1i + 5]
    set s2j [expr $s1i + 6]
    set b2j [expr $s1i + 7]

    set rs2i [part $s2i pr pos]
    set rs1i [part $s1i pr pos]
    set rb2i [part $b2i pr pos]
    set rb1i [part $b1i pr pos]

    set rs2j [part $s2j pr pos]
    set rs1j [part $s1j pr pos]
    set rb2j [part $b2j pr pos]
    set rb1j [part $b1j pr pos]

    set rcci [vecadd 1.0 $rs1i -1.0 $rs2i]
    set rcci_l [veclen $rcci]

    set rccj [vecadd 1.0 $rs1j -1.0 $rs2j]
    set rccj_l [veclen $rccj]
   
    set rsi [vecadd 0.5 $rs1i 0.5 $rs2i]
    set rsj [vecadd 0.5 $rs1j 0.5 $rs2j]
    set rss [vecadd 1.0 $rsi -1.0 $rsj]
    set rss_l [veclen $rss]

    set rss1 [vecadd 1.0 $rs1i -1.0 $rs1j]
    set rss2 [vecadd 1.0 $rs2i -1.0 $rs2j]
    set rss1_l [veclen $rss1]
    set rss2_l [veclen $rss2]

    set rcb1i [vecadd 1.0 $rs1i -1.0 $rb1i]
    set rcb1i_l [veclen $rcb1i]

    set rcb2i [vecadd 1.0 $rs1i -1.0 $rb1i]
    set rcb2i_l [veclen $rcb1i]

    set rcb1j [vecadd 1.0 $rs1j -1.0 $rb1j]
    set rcb1j_l [veclen $rcb1j]

    set rcb2j [vecadd 1.0 $rs1j -1.0 $rb1j]
    set rcb2j_l [veclen $rcb1j]

    set n1i [veccross $rcci $rcb1i]
    set n2i [veccross $rcci $rcb2i]
    set n1j [veccross $rccj $rcb1j]
    set n2j [veccross $rccj $rcb2j]
    set n1i_l [veclen $n1i]

    set ani [vecadd 0.5 $n1i 0.5 $n2i]
    set anj [vecadd 0.5 $n1j 0.5 $n2j]
    set ani_l [veclen $ani]
    set anj_l [veclen $anj]

    set cos_tilt [expr [vecdot $ani $anj]/($ani_l*$anj_l)]

    set rccj_parallel [expr [vecdot $n1i $rccj]/$n1i_l]
    
    set rp [vecadd 1.0 $rccj [expr -$rccj_parallel] $n1i]
    set rp_l [veclen $rp]
    
    set cos_theta [expr [vecdot $rp $rcci] /($rp_l*$rcci_l)]

    return [list [expr 180.*acos($cos_theta)/$Pi] [expr 0.5*($rss1_l+$rss2_l)] [expr 180.*acos($cos_tilt)/$Pi ]]
}

proc analyze_stacking_all {} {
    set theta_twist_avg 0.0
    set rss_avg 0.0
    set theta_tilt_avg 0.0
    set n 0
    for { set i 0 } { $i <= [expr [setmd max_part] - 4] } { incr i } {
	set type [part $i pr type]
	if { $type == 0 } {
	    set vals [analyze_stacking $i]
	    set theta_twist_avg [expr $theta_twist_avg + [lindex $vals 0]]
	    set rss_avg [expr $rss_avg + [lindex $vals 1]]
	    set theta_tilt_avg [expr $theta_tilt_avg + [lindex $vals 2]]
	    incr n
	}    
    }
    return [list [expr $theta_twist_avg/$n] [expr $rss_avg/$n] [expr $theta_tilt_avg/$n]]
}


