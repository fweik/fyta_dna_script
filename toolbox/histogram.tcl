## construct a 1D histogram including over- and underbin
proc histogram1D { min max nbins histoname} {
  set vector [ list ]
  for { set i 0 } { $i < [ expr $nbins + 2] } { incr i } {
    lappend vector 0.
  }
  if { $min > $max } {
    set temp $max
    set max $min
    set min $temp
  }
  upvar $histoname histo 
  set histo [ list ]
  lappend histo $min
  lappend histo $max
  lappend histo $nbins
  lappend histo $vector
#  return $histo
#  return { $min $max $nbins $vector }
}

## add point to a 1D histogram
proc histogram1D_addPoint { histo_name x } {
  ## unpack histogram
  upvar $histo_name histo
  set min [ lindex $histo 0 ]
  set max [ lindex $histo 1 ]
  set nbins [ lindex $histo 2 ]
  set vector [ lindex $histo 3 ]

  if { $x > $max } {
    lset vector [ expr int($nbins + 1) ] [ expr [ lindex $vector  [ expr int($nbins + 1) ] ]+1  ]
    lset histo 3 $vector
    return 
  }
  if { $x < $min } {
    lset vector [ expr int($nbins) ] [ expr [ lindex $vector  [ expr int($nbins) ] ] +1 ]
    lset histo 3 $vector
    return 
  }
  ## else
  set bin [ expr int( floor ( ($x-$min)/($max-$min)*$nbins ) ) ]
  set new_value [ expr [ lindex $vector $bin ] + 1. ]
  lset vector $bin $new_value
  lset histo 3 $vector
  #[ expr [ lindex $vector $bin ] + 1. ]
}

proc histogram1D_printToFile { histo_name fname } {
  upvar $histo_name histo
  set ofile [ open $fname "w" ]
  set min [ lindex $histo 0 ]
  set max [ lindex $histo 1 ]
  set nbins [ lindex $histo 2 ]
  set vector [ lindex $histo 3 ]
  set binlength [ expr ($max - $min)/$nbins ]
  for { set i 0 } { $i < $nbins } { incr i } {
    puts $ofile "[ expr ($i + 0.5) * $binlength ] [ lindex $vector $i ]"
  }
  puts $ofile "under [ lindex $vector $nbins ]" 
  puts $ofile "over [ lindex $vector [ expr $nbins + 1 ] ]" 
  close $ofile
}

proc analyse_distribution_x { histo_name minpart maxpart} {
  upvar $histo_name histo
  for { set i $minpart } { $i < $maxpart } { incr i } {
    histogram1D_addPoint histo [ lindex [ part $i print pos ] 0 ] 
  }
}

proc pow2 { var } {
  return [ expr $var*$var ] 
}

proc radial_density_profile { density_profile_name resolution_x resolution_r
    {  xsize "[ lindex [ setmd box_l ] 0 ]" }  
    { rsize "[ min [ expr 0.5* [ lindex [ setmd box_l ] 1 ] ] [ expr 0.5*[ lindex [ setmd box_l ] 2 ] ] ]" }  
    { xoffset 0. }   } {
  set pi 3.1415
  upvar $density_profile_name density_profile
  set centery [ expr 0.5* [ lindex [ setmd box_l ] 1 ] ]
  set centerz [ expr 0.5* [ lindex [ setmd box_l ] 2 ] ]
  set xbins [ expr int ( floor ( $xsize / $resolution_x ) ) ]
  set rbins [ expr int ( floor ( $rsize / $resolution_r ) ) ]
  set resolution_x [ expr $xsize / $xbins ]
  set resolution_r [ expr $rsize / $rbins ]

  set data [ list ]
  set positions [ list ]
  set volumes [ list ]
  for { set i 0 } { $i < [ expr $xbins ] } { incr i } {
    for { set j 0 } { $j < $rbins } { incr j } {
      lappend data 0.
      lappend positions  [ list [ expr $xoffset + ($i + 0.5) * $resolution_x ] [ expr ($j+0.5)*$resolution_r ] ]  
      lappend volumes [ expr $resolution_x*$pi* ([ pow2 [ expr ($j+1)*$resolution_r ] ] - [ pow2 [ expr $j*$resolution_r ] ] ) ]
    }
  }
  set n_sweeps 0
  set density_profile [ list ]
  lappend density_profile $resolution_x 
  lappend density_profile $resolution_r 
  lappend density_profile $xbins
  lappend density_profile $rbins
  lappend density_profile $centery
  lappend density_profile $centerz
  lappend density_profile $n_sweeps
  lappend density_profile $data
  lappend density_profile $positions
  lappend density_profile $volumes 
  lappend density_profile $xoffset
}

proc radial_density_profile_add_particles { density_profile_name type } {
  upvar $density_profile_name density_profile
  set resolution_x [ lindex $density_profile 0 ] 
  set resolution_r [ lindex $density_profile 1 ] 
  set xbins [ lindex $density_profile 2 ]
  set rbins [ lindex $density_profile 3 ]
  set centery [ lindex $density_profile 4 ]
  set centerz [ lindex $density_profile 5 ]
  set n_sweeps [ lindex $density_profile 6 ]
  set data [ lindex $density_profile 7 ]
  set positions [ lindex $density_profile 8 ]
  set volumes [ lindex $density_profile 9 ] 
  set xoffset [ lindex $density_profile 10 ]

  for { set i 0 } { $i < [ setmd n_part ] } { incr i } {
    if { [ part $i print type ] == $type } {
      set pos [ part $i print folded_pos ]
      set x [ lindex $pos 0 ]
      set ry [ expr [ lindex $pos 1 ] - $centery ]
      set rz [ expr [ lindex $pos 2 ] - $centerz ]
      set r [ expr sqrt ( [ pow2 $ry ] + [ pow2 $rz ] ) ]
      set xbin [ expr int(floor( (($x-$xoffset)/$resolution_x )))]
      set rbin [ expr int(floor( $r/$resolution_r ))]
      if { $rbin < $rbins && $xbin < $xbins && $xbin >= 0 } {
        set bin [ expr $xbin*$rbins + $rbin ] 
        lset data $bin [ expr [ lindex $data $bin ] + 1 ]
      }
    }
  }
  incr n_sweeps 
  
  lset density_profile 0 $resolution_x
  lset density_profile 1 $resolution_r
  lset density_profile 2 $xbins
  lset density_profile 3 $rbins
  lset density_profile 4 $centery
  lset density_profile 5 $centerz
  lset density_profile 6 $n_sweeps
  lset density_profile 7 $data
  lset density_profile 8 $positions
  lset density_profile 9 $volumes
}


proc radial_density_profile_print_to_file { density_profile_name fname } {
  upvar $density_profile_name density_profile
  set resolution_x [ lindex $density_profile 0 ] 
  set resolution_r [ lindex $density_profile 1 ] 
  set xbins [ lindex $density_profile 2 ]
  set rbins [ lindex $density_profile 3 ]
  set centery [ lindex $density_profile 4 ]
  set centerz [ lindex $density_profile 5 ]
  set n_sweeps [ lindex $density_profile 6 ]
  set data [ lindex $density_profile 7 ]
  set positions [ lindex $density_profile 8 ]
  set volumes [ lindex $density_profile 9 ] 
  set xoffset [ lindex $density_profile 10 ]
  set ofile [ open $fname "w" ]
  for { set i 0 } { $i < [ llength $data ] } { incr i } {
    puts $ofile "[ lindex [ lindex $positions $i ] 0 ] [ lindex [ lindex $positions $i ] 1 ] [ expr [ lindex $data $i ] / [ lindex $volumes $i ] / $n_sweeps / 6.e-4]"
  }
  close $ofile
}

proc radial_j_density_profile { density_profile_name resolution_x resolution_r
    {  xsize "[ lindex [ setmd box_l ] 0 ]" }  
    { rsize "[ min [ expr 0.5* [ lindex [ setmd box_l ] 1 ] ] [ expr 0.5*[ lindex [ setmd box_l ] 2 ] ] ]" }  
    { xoffset 0. }   } {
  set pi 3.1415
  upvar $density_profile_name density_profile
  set centery [ expr 0.5* [ lindex [ setmd box_l ] 1 ] ]
  set centerz [ expr 0.5* [ lindex [ setmd box_l ] 2 ] ]
  set xbins [ expr int ( floor ( $xsize / $resolution_x ) ) ]
  set rbins [ expr int ( floor ( $rsize / $resolution_r ) ) ]
  set resolution_x [ expr $xsize / $xbins ]
  set resolution_r [ expr $rsize / $rbins ]

  set data [ list ]
  set positions [ list ]
  set volumes [ list ]
  for { set i 0 } { $i < [ expr $xbins ] } { incr i } {
    for { set j 0 } { $j < $rbins } { incr j } {
      lappend data 0.
      lappend positions  [ list [ expr $xoffset + ($i + 0.5) * $resolution_x ] [ expr ($j+0.5)*$resolution_r ] ]  
      lappend volumes [ expr $resolution_x*$pi* ([ pow2 [ expr ($j+1)*$resolution_r ] ] - [ pow2 [ expr $j*$resolution_r ] ] ) ]
    }
  }
  set n_sweeps 0
  set density_profile [ list ]
  lappend density_profile $resolution_x 
  lappend density_profile $resolution_r 
  lappend density_profile $xbins
  lappend density_profile $rbins
  lappend density_profile $centery
  lappend density_profile $centerz
  lappend density_profile $n_sweeps
  lappend density_profile $data
  lappend density_profile $positions
  lappend density_profile $volumes 
  lappend density_profile $xoffset
}

proc radial_j_density_profile_add_particles { density_profile_name type } {
  upvar $density_profile_name density_profile
  set resolution_x [ lindex $density_profile 0 ] 
  set resolution_r [ lindex $density_profile 1 ] 
  set xbins [ lindex $density_profile 2 ]
  set rbins [ lindex $density_profile 3 ]
  set centery [ lindex $density_profile 4 ]
  set centerz [ lindex $density_profile 5 ]
  set n_sweeps [ lindex $density_profile 6 ]
  set data [ lindex $density_profile 7 ]
  set positions [ lindex $density_profile 8 ]
  set volumes [ lindex $density_profile 9 ] 
  set xoffset [ lindex $density_profile 10 ]

  for { set i 0 } { $i < [ setmd n_part ] } { incr i } {
    if { [ part $i print type ] == $type } {
      set pos [ part $i print folded_pos ]
      set x [ lindex $pos 0 ]
      set ry [ expr [ lindex $pos 1 ] - $centery ]
      set rz [ expr [ lindex $pos 2 ] - $centerz ]
      set r [ expr sqrt ( [ pow2 $ry ] + [ pow2 $rz ] ) ]
      set xbin [ expr int(floor( (($x-$xoffset)/$resolution_x )))]
      set rbin [ expr int(floor( $r/$resolution_r ))]
      if { $rbin < $rbins && $xbin < $xbins && $xbin >= 0 } {
        set bin [ expr $xbin*$rbins + $rbin ] 
        lset data $bin [ expr [ lindex $data $bin ] + [ lindex [ part $i print v ] 0 ] ]
      }
    }
  }
  incr n_sweeps 
  
  lset density_profile 0 $resolution_x
  lset density_profile 1 $resolution_r
  lset density_profile 2 $xbins
  lset density_profile 3 $rbins
  lset density_profile 4 $centery
  lset density_profile 5 $centerz
  lset density_profile 6 $n_sweeps
  lset density_profile 7 $data
  lset density_profile 8 $positions
  lset density_profile 9 $volumes
}


proc radial_j_density_profile_print_to_file { density_profile_name fname } {
  upvar $density_profile_name density_profile
  set resolution_x [ lindex $density_profile 0 ] 
  set resolution_r [ lindex $density_profile 1 ] 
  set xbins [ lindex $density_profile 2 ]
  set rbins [ lindex $density_profile 3 ]
  set centery [ lindex $density_profile 4 ]
  set centerz [ lindex $density_profile 5 ]
  set n_sweeps [ lindex $density_profile 6 ]
  set data [ lindex $density_profile 7 ]
  set positions [ lindex $density_profile 8 ]
  set volumes [ lindex $density_profile 9 ] 
  set xoffset [ lindex $density_profile 10 ]
  set ofile [ open $fname "w" ]
  for { set i 0 } { $i < [ llength $data ] } { incr i } {
    puts $ofile "[ lindex [ lindex $positions $i ] 0 ] [ lindex [ lindex $positions $i ] 1 ] [ expr [ lindex $data $i ] / [ lindex $volumes $i ] / $n_sweeps ]"
  }
  close $ofile
}

