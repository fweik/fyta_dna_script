
proc setup_stiff_pe_rod { centerx centery centerz axis_x axis_y axis_z length n_monomers charge startid type } {
  if { $n_monomers == 0 } {
    return
  }
  set l_incr [ expr $length/$n_monomers ]
  set charge_per_monomer [ expr $charge / $n_monomers ]
  ## normalize axis vector 
  set l [ expr sqrt($axis_x*$axis_x + $axis_y*$axis_y + $axis_z * $axis_z) ]
  set axis_x [ expr $axis_x / $l ]
  set axis_y [ expr $axis_y / $l ]
  set axis_z [ expr $axis_z / $l ]
  for { set i 0 } { $i < $n_monomers } { incr i } {
    part [ expr $startid + $i ] pos [ expr $centerx + $i*$l_incr*$axis_x ] [ expr $centery + $i*$l_incr*$axis_y ] [ expr $centerz + $i*$l_incr*$axis_z ] q $charge_per_monomer fix 1 1 1 type $type
  }
}

proc calc_total_force { startid no } {
  set f_x 0
  set f_y 0
  set f_z 0 
  for { set i $startid } { $i < [ expr $startid + $no ] } { incr i } {
    set force [ part $i print f ]
    set f [ split $force ]
    set f_x [ expr $f_x + [ lindex $f 0 ] ]
    set f_y [ expr $f_y + [ lindex $f 1 ] ]
    set f_z [ expr $f_z + [ lindex $f 2 ] ]
  }
  return [ list $f_x $f_y $f_z ]
}

proc nix { } {
  puts "hi"
}
