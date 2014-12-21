
proc mean { l } {
  set sum 0.
  foreach x $l {
    set sum [ expr $sum + $x ]
  }
  return [ expr $sum / [ llength $l ] ]
}

proc stdev { l } {
  set sum 0.
  set themean [ mean $l ]
  foreach x $l {
    set sum [ expr $sum  + ($x-$themean)*($x-$themean) ]
  }
  return [ expr sqrt($sum / [ llength $l] ) ]
}

proc lmax { thelist } {
  set themax [ lindex $thelist 0 ]
  foreach  x  $thelist {
    set themax [ max $themax $x ]
  }
  return $themax
}

proc lmin { thelist } {
  set themin [ lindex $thelist 0 ]
  foreach { x } { $thelist } { 
    set themin [ min $themin $x ]
  }
  return $themin
}

proc maxindex { thelist } {
  set themax [ lindex $thelist 0 ]
  set maxindex 0
  for  { set i 1 } { $i < [ llength $thelist ] } { incr i } { 
    if { [ lindex $thelist $i ] > $themax } {
      set themax [ lindex  $thelist $i ]
      set maxindex $i
    }
  }
  return $maxindex
}

proc checkintegrate { steps } {
  puts "checkintegrate"
  set part_positions_old [ list ]
  set part_positions [ list ]
  set distances [ list ]
  for { set i 0 } { $i < [ setmd n_part ] } { incr i } {
    lappend part_positions_old [ part $i print pos ]
  }
  integrate $steps
  for { set i 0 } { $i < [ setmd n_part ] } { incr i } {
    lappend part_positions [ part $i print pos ]
  }
  for { set i 0 } { $i < [ setmd n_part ] } { incr i } {
    set temp1 0
    set temp2 [ expr [ lindex [ lindex $part_positions [ expr $i ] ] 0 ] - [ lindex [ lindex $part_positions_old [ expr $i ] ] 0 ] ]
    set temp1 [ expr $temp1 + $temp2*$temp2 ]
    set temp2 [ expr [ lindex [ lindex $part_positions [ expr $i ] ] 1 ] - [ lindex [ lindex $part_positions_old [ expr $i ] ] 1 ] ]
    set temp1 [ expr $temp1 + $temp2*$temp2 ]
    set temp2 [ expr [ lindex [ lindex $part_positions [ expr $i ] ] 2 ] - [ lindex [ lindex $part_positions_old [ expr $i ] ] 2 ] ]
    set temp1 [ expr $temp1 + $temp2*$temp2 ]
    lappend distances [ expr sqrt($temp1) ]
  }
  set maxi [ maxindex $distances ]
  if { [ lindex $distances $maxi ] > 0.1 } {
    puts "$maxi [ lindex $distances $maxi ] old [ lindex $part_positions_old $maxi ] new [ lindex $part_positions $maxi ]" 
  }
}

if {[ info commands galileiTransformParticles ] == "" } {
proc galileiTransformParticles {} {
  galilei_transform
}
}

if {[ info commands galilei_transform ] == "" } {
proc galilei_transform {} {
  galileiTransformParticles
}
}

proc grow_lj { ints} { 
  set timestep [ setmd time_step ]
  set n_part [ setmd n_part ]
  set n_int [ llength $ints ]
  set espresso_int [ list ]
  set sigma [ list ]
  set shifts [ list ]
  puts [ inter ]
  puts "[ analyze energy ]"
  puts "mindist [analyze mindist ]"
  for { set i 0 } { $i<$n_int } { incr i } {
    set this_int [ lindex $ints $i ]
    set this_espresso_int [ inter [ lindex $this_int 0 ]  [ lindex $this_int 1 ] ]
    lappend sigma [ lindex $this_espresso_int 4 ]
    lappend shifts [ lindex $this_espresso_int 7 ]
    lappend espresso_int $this_espresso_int
  }
  set currentsigma [ expr [ analyze mindist ] / 2. ]
  set currenttimestep [ expr $timestep * $currentsigma / [ lmax $sigma ] ]
  setmd time_step $currenttimestep
  puts "starting with timestep $currenttimestep sigma $currentsigma"
  ### setup initial interactions:
  for { set i 0 } { $i<$n_int } { incr i } {
    set int [ lindex $espresso_int $i ]
    inter [ lindex $int 0 ] [ lindex $int 1 ] [ lindex $int 2 ] [ lindex $int 3 ] [ min $currentsigma [ lindex $sigma $i ] ] [ lindex $int 5 ] [ lindex $int 6 ] [ min $currentsigma [ lindex $shifts $i ] ]
  }
  puts "[ analyze energy ]"
#  puts "sigma $sigma"
  integrate 1
  while { $currentsigma < [ lmax [ concat $sigma $shifts ] ] } {
    integrate 10
    while {  [ expr [ analyze energy kinetic ]  /[degrees_of_freedom]/$n_part ] > 0.6 } {
      puts "hallo [ expr [ analyze energy kinetic ] / [ degrees_of_freedom ] [ analyze energy ]"
      imd positions
      integrate 10
      #galilei_transform
      galileiTransformParticles
    }
    set currentsigma [ expr $currentsigma * 1.05 ]
    set currenttimestep [ expr $currenttimestep * 1.05 ]
    puts "increasing sigma to $currentsigma"
    for { set i 0 } { $i<$n_int } { incr i } {
      set int [ lindex $espresso_int $i ]
      inter [ lindex $int 0 ] [ lindex $int 1 ] [ lindex $int 2 ] [ lindex $int 3 ] [ min $currentsigma [ lindex $sigma $i ]  ] [ lindex $int 5 ] [ lindex $int 6 ] [ min $currentsigma [ lindex $shifts $i ] ] 
    }
    setmd time_step $currenttimestep
    integrate 10 
#    imd positions

  }

  setmd time_step $timestep
#  puts "total energy per degree of freedom [ expr [ lindex [ lindex [ analyze energy ] 0 ] 1 ]  /3./$n_part ]"
#  puts "kinetic energy per particle/1.5 [ expr [ lindex [ lindex [ analyze energy ] 1 ] 1 ]  /3./$n_part ]"
} 

proc violates_constraints { x y z } {
  if { $x < 1. || $x > 6 } {
    return 1
  } else {
    return 0
  }
}

proc violates_pore_constraints { id cx cy cz nx ny nz rad len pore_id } {
   set pos [ part $id print pos ] 
   set posx [ lindex $pos 0 ]
   set posy [ lindex $pos 1 ]
   set posz [ lindex $pos 2 ]
   set parttype [ part $id print type ]
   set int [ inter $parttype $pore_id ]
   set sigma 0.1
   set norm [ expr sqrt($nx*$nx + $ny*$ny + $nz*$nz) ]
   set nx [ expr $nx/$norm ]
   set ny [ expr $ny/$norm ]
   set nz [ expr $nz/$norm ]
   ## sqdistance from pore center projected to pore axis:
   set dist_n  [ expr ($posx-$cx)*$nx+($posy-$cy)*$ny+($posz-$cz)*$nz ] 
   if { $dist_n < [ expr - ($len + $sigma) ] || $dist_n > [ expr + ($len+ $sigma) ] } {
     return 0
   } else { 
    ## square distance from the center
     set sq_c_dist [ expr ($posx-$cx)*($posx-$cx)+($posy-$cy)*($posy-$cy)+($posz-$cz)*($posz-$cz) ]
  ## distance from axis by Pythagorean theorem
     set r_dist [ expr sqrt($sq_c_dist - $dist_n*$dist_n) ]
  ## if distance from pore axis < ($rad-$sigma)  
     if { $r_dist < [ expr $rad - $sigma ] } {
       return 0
     } else  {
       return 1
     }
  }
 
}

proc violates_tube_constraints { id cx cy cz nx ny nz rad } {
   set pos [ part $id print pos ] 
   set posx [ lindex $pos 0 ]
   set posy [ lindex $pos 1 ]
   set posz [ lindex $pos 2 ]
   set parttype [ part $id print type ]
   set sigma 0.1
   set norm [ expr sqrt($nx*$nx + $ny*$ny + $nz*$nz) ]
   set nx [ expr $nx/$norm ]
   set ny [ expr $ny/$norm ]
   set nz [ expr $nz/$norm ]
   ## sqdistance from pore center projected to pore axis:
   set dist_n  [ expr ($posx-$cx)*$nx+($posy-$cy)*$ny+($posz-$cz)*$nz ] 
 ## square distance from the center
   set sq_c_dist [ expr ($posx-$cx)*($posx-$cx)+($posy-$cy)*($posy-$cy)+($posz-$cz)*($posz-$cz) ]
## distance from axis by Pythagorean theorem
   set r_dist [ expr sqrt($sq_c_dist - $dist_n*$dist_n) ]
## if distance from pore axis < ($rad-$sigma)  
   if { $r_dist < [ expr $rad - $sigma ] } {
     return 0
   } else  {
     return 1
   }
}


proc init_parameters  parameterlist  {
  global globallist 
  set paramlist_flat [ list ]
  if { [ get_command_line_argument "paramfile" 0 ] == 0 } {
    for { set i 0 } { $i < [ llength $parameterlist ] } { incr i } {
      global [ lindex [ lindex $parameterlist $i ] 0 ]
      lappend paramlist_flat [ lindex [ lindex $parameterlist $i ] 0 ]
      set_global [ lindex [ lindex $parameterlist $i ] 0 ] [ get_command_line_argument [ lindex [ lindex $parameterlist $i ] 0 ] [ lindex [ lindex $parameterlist $i ] 1 ] ]
    }
    save_paramfile params $paramlist_flat
  } else { 
    set paramfilename  [ get_command_line_argument "paramfile" "" ] 
    set paramlist_flat [ list ]
    for { set i 0 } { $i < [ llength $parameterlist ] } { incr i } {
      global [ lindex [ lindex $parameterlist $i ] 0 ]
      lappend paramlist_flat [ lindex [ lindex $parameterlist $i ] 0 ]
    }
    load_paramfile $paramfilename $paramlist_flat
    puts "loaded parameterfile $paramfilename"
  }

}

#set box_l 10.
#setmd box_l $box_l $box_l $box_l
#setmd periodic 1 1 1
#set density 0.6
#
#set lj1_eps     1.0
#set lj1_sig     1.0
#set lj1_cut     1.12246
#set lj1_shift   [calc_lj_shift $lj1_sig $lj1_cut]
#
## Integration parameters
##############################################################
#
#setmd time_step 0.01
#setmd skin      0.4
#thermostat langevin 1.0 1.0
#
#inter 0 0 lennard-jones $lj1_eps $lj1_sig $lj1_cut $lj1_shift 0
#inter 0 1 lennard-jones $lj1_eps $lj1_sig $lj1_cut $lj1_shift 0
#inter 0 2 lennard-jones $lj1_eps [ expr 2 * $lj1_sig ] $lj1_cut $lj1_shift 0
#inter 1 2 lennard-jones $lj1_eps $lj1_sig $lj1_cut $lj1_shift 0
#inter 2 2 lennard-jones $lj1_eps [ expr 3 * $lj1_sig ] $lj1_cut $lj1_shift 0
#
#set volume [expr $box_l*$box_l*$box_l]
#set n_part [expr floor($volume*$density)]
#
##constraint wall normal 1. 0. 0. dist -0 type 1
##constraint wall normal -1. 0. 0. dist -7 type 1
#constraint pore center 5. 5. 5. axis 1. 0. 0. radius 2. length 1. type 1
#puts [ inter 0 0 ]
#for {set i 0} { $i < $n_part } {incr i} {
#  set violates 1
#  while { $violates } {
#    set posx [expr $box_l*[t_random]]
#    set posy [expr $box_l*[t_random]]
#    set posz [expr $box_l*[t_random]]
#    part $i pos $posx $posy $posz type [expr int(2*floor($i/($n_part/2)))] 
#    set violates [ violates_pore_constraints $i 5. 5. 5. 1. 0. 0. 2. 2. 1 ]
#  }
#}
##prepare_vmd_connection "test"
#
#puts "starting to grow w mindist [ analyze mindist ]"
#grow_lj { { 0 0 } { 0 2 } { 2 2 } } 
#puts "stopping to grow w mindist [ analyze mindist ]"
#
#writepsf "chrystal.psf"
#
proc slow_integrate steps {
  global counter blackbox
  integrate $steps
  return
  integrate 0
  for { set i 0 } { $i < $steps } { incr i } {
    for { set j 0 } { $j < [ setmd n_part ] } { incr j } {
      set f [ part $i print f ]
      for { set k 0 } { $k < 3 } { incr k } {
        if { [ expr abs([ lindex $f $k ]) ] > 1e3 } {
          puts "dangerous integration for particle $j, [ part $j ]"
        }
      }
    }
#    if { $counter == 38 && $i == 6} {
#      set ofile [ open "died_config.dat" "w" ]
#      for { set k 0 } { $k < [ setmd n_part ] } { incr k } {
#        puts $ofile [ part $k ]
#      }
#      close $ofile
#    } 
    integrate 1
#    puts $blackbox [ part 106 ]
    for { set j 0 } { $j < [ setmd n_part ] } { incr j } {
      set f [ part $i print f ]
      for { set k 0 } { $k < 3 } { incr k } {
        if { [ expr abs([ lindex $f $k ]) ] > 1e3 } {
          puts "this was a dangerous integration for particle $j, [ part $j ]"
        }
      }
    }
  }  
}

proc safe_integrate steps {
  global centerx centery centerz pore_diameter pore_length box_l lj_ion_ion_sig blackbox
  catch { integrate $steps } err
  if { $err != "" } {
    puts "found error $err"
    set errlist [ split $err { \}\}} ]
    if { [ lindex $errlist 2 ] == "\{063" } {
      set part_no [ lindex $errlist 9 ]
      puts "saving particle $part_no"
      set violates 1
      while { $violates } {
        set posx [expr $box_l*[t_random]]
        set posy [expr $box_l*[t_random]]
        set posz [expr $box_l*[t_random]]
        part $part_no pos $posx $posy $posz v 0. 0. 0.
        set violates [ violates_pore_constraints $part_no $centerx $centery $centerz 1. 0. 0. [ expr $pore_diameter/2. ] $pore_length 1 ]
        if { $violates > 0 } {
          safe_integrate 0 
          set f  [ part $part_no print f ]
          set maxf  -1
          for { set i 0 } { $i < 3 } { incr i } {
            if { [expr abs( [ lindex $f $i ]) > 1e2  ]  }  {
             set violates 1
            }
          }
        }
      }
    } else { 
    puts "could not probhibit $err" 
#    close $blackbox
    raise $err
  }
  puts "I saved  a particle violating the constraints"
  safe_integrate $steps
  }
}

