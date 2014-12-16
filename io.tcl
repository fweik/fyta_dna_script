#
# (C) Copyright 2014, Damaris Holder, Florian Weik
#
# All rights reserved. 
#
# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#

proc write_configuration { fd } {
    for { set i 0 } { $i <= [setmd max_part] } { incr i } {
	puts $fd [part $i pr pos]
    }
}

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
	puts "Not enough particles in file '$filename' to read [expr $n_particles] basepairs."
	exit
    }   
}

proc read_sequence { filename } {
    set sequence [open $filename "r"]

    set ladderlist []
    while {1} {
	set line [gets $sequence]
	if {[eof $sequence]} {
	    close $sequence
	    break
	}
	lappend ladderlist $line
    }
    return $ladderlist
}
