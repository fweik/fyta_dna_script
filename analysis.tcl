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

proc analyze_stacking { s1i } {
    set Pi 3.14159

    set s2i [expr $s1i + 2]
    set b1i [expr $s1i + 1]
    set s1j [expr $s1i + 4]
    set s2j [expr $s1i + 6]

    set x1 [part $s2i pr pos]
    set x2 [part $s1i pr pos]

    set rccix [expr ([lindex $x2 0] - [lindex $x1 0])]
    set rcciy [expr ([lindex $x2 1] - [lindex $x1 1])]
    set rcciz [expr ([lindex $x2 2] - [lindex $x1 2])]
    set rccil [expr sqrt($rccix*$rccix + $rcciy*$rcciy + $rcciz*$rcciz)]

    set x1 [part $s2j pr pos]
    set x2 [part $s1j pr pos]

    set rccjx [expr ([lindex $x2 0] - [lindex $x1 0])]
    set rccjy [expr ([lindex $x2 1] - [lindex $x1 1])]
    set rccjz [expr ([lindex $x2 2] - [lindex $x1 2])]
    set rccjl [expr sqrt($rccjx*$rccjx + $rccjy*$rccjy + $rccjz*$rccjz)]

    set x1 [part $b1i pr pos]
    set x2 [part $s1i pr pos]
    set dxcb [expr ([lindex $x2 0] - [lindex $x1 0])]
    set dycb [expr ([lindex $x2 1] - [lindex $x1 1])]
    set dzcb [expr ([lindex $x2 2] - [lindex $x1 2])]
    set rcb1 [expr sqrt($dxcb*$dxcb + $dycb*$dycb + $dzcb*$dzcb)]
    
    # n1 = rcci x rcb1
    set nx [expr $rcciy*$dzcb - $rcciz*$dycb]
    set ny [expr $rcciz*$dxcb - $rccix*$dzcb]
    set nz [expr $rccix*$dycb - $rcciy*$dxcb]
    set nl [expr sqrt($nx*$nx + $ny*$ny + $nz*$nz)]

    set nx [expr $nx/$nl]
    set ny [expr $ny/$nl]
    set nz [expr $nz/$nl]

    set rccj_parallel [expr ($nx*$rccjx+$ny*$rccjy+$nz*$rccjz)]
    
    set rpx [expr $rccjx - $rccj_parallel*$nx]
    set rpy [expr $rccjy - $rccj_parallel*$ny]
    set rpz [expr $rccjz - $rccj_parallel*$nz]
    set rpl [expr sqrt($rpx*$rpx + $rpy*$rpy + $rpz*$rpz)]
    
    set cos_theta [expr ($rpx*$rccix + $rpy*$rcciy + $rpz*$rcciz)/($rpl*$rccil)]

    return [expr 180.*acos($cos_theta)/$Pi]
}

proc analyze_stacking_all {} {
    set theta_avg 0.0
    set n 0
    for { set i 0 } { $i <= [expr [setmd max_part] - 4] } { incr i } {
	set type [part $i pr type]
	if { $type == 0 } {
	    set theta_avg [expr $theta_avg + [analyze_stacking $i]]
	    incr n
	}    
    }
    return [expr $theta_avg/$n]
}


