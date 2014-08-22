
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
source interactions.tcl

# General MD parameters
set time_step 0.1
set total_int_steps 2000000
set steps_per_loop 1000
set equilibration_loops 1000
set skin 1.0

# Langevin parameters
# kT = 0.025eV ~ 300K (k_B = 8.617e-5 eV/K)
set kT 0.026
set gamma 1.0

# Output options
set vmd "no"
set vtf "no"
set vtf_filename "/work/fweik/dna.vtf"

# Molecule
set n_basepairs 200
set configuration_filename "configurations/1000bp_conf.dat"
set sequence_filename "sequences/Sequence_polyAT.dat"
# Fix one end of molecule?
set fix_lower_end "no"
# External forces on molecule
#Sheer force in +x direction
set ext_force_sheer 0.0
#Stretch force in +z direction
set ext_force_stretch 0.0

# Analysis
set analyse_persistence_length "no"
set persistence_length_file "peristence_length.dat"
set analyse_chain_parameters "no"
set chain_parameter_file "chain_parameters.dat"
set analyse_end_to_end_dist "yes"
set end_to_end_file "ete.dat"
set analyse_avg_end_to_end_dist "yes"
set avg_end_to_end_file "avg_ete.dat" 

# Electrostatics
set lB 561
set lambdaDB 9.6
set alpha [expr -14.23]

# Check for command lineparameters

if { $argc == 3 } {
    set n_basepairs [lindex $argv 0]
    set lambdaDB [lindex $argv 1]
    set avg_end_to_end_file [lindex $argv 2]
    if { $analyse_avg_end_to_end_dist == "yes" } {
	set end_to_end_file "ete_[expr $n_basepairs]_[expr $lambdaDB].dat"
    }
}

puts "n_basepairs $n_basepairs, end_to_end_file $end_to_end_file avg_end_to_end_file $avg_end_to_end_file"
puts "bejrrum length $lB debye length $lambdaDB"

# Box geometry
# Length along the molecule
set box_z [expr 4*$n_basepairs + 500.]
# Other directions
set box_xy 500.
# Shift along molecule axis
set zshift 250.
# Shift in other directions
set center_xy [expr 0.5*$box_xy]

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

if { $n_basepairs > [expr [llength $ladderlist]] } {
    puts "Sequence '$sequence_filename' is too short for molecule of length $n_basepairs."
    exit
}

# Configure molecule
set_charges

set_masses $ladderlist

# electrostatic interactions
setup_electrostatics $lB $lambdaDB [expr 5*$lambdaDB] $alpha

setup_bonded_interactions $ladderlist

if { $fix_lower_end == "yes" } {
    part 0 fix
    part 2 fix
}

#for { set i 0 } { $i <= [setmd max_part] } { incr i } {
#    puts [part $i]
#}

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

if { $analyse_persistence_length == "yes" } {
    set pers [open $persistence_length_file "w"]
}
if { $analyse_chain_parameters == "yes" } {
    set fo [open $chain_parameter_file "w"]
}

if { $analyse_end_to_end_dist == "yes" } {
    set e2e [open $end_to_end_file "w"]    
}

if { $analyse_avg_end_to_end_dist == "yes" } {
    set samples 0
    set R_avg 0.0
    set R2_avg 0.0
    set l_avg 0.0
    set l2_avg 0.0
}

set largest 0

for { set i 0 } { $i < $equilibration_loops } { incr i } {
    puts "Equilibration Loop [expr $i+1] of $equilibration_loops"
#    analyze_bps
#    puts "<theta_twist> = [analyze_stacking_all]"
    integrate $steps_per_loop    
}

for { set i 0 } { $i <= $int_loops } { incr i } {
    puts "Loop $i of $int_loops, (time [format %.2g [expr $i*$time_step*$steps_per_loop]] of [format %.2g [expr $total_int_steps*$time_step]])."
    puts [time { integrate $steps_per_loop }]

    if { $analyse_chain_parameters == "yes" } {
	puts $fo [analyze_bps]
    }

    if { $analyse_persistence_length == "yes" } {
	puts $pers [analyze_pl]
    }

    if { $analyse_end_to_end_dist == "yes" } {
	puts $e2e [analyze_end_to_end_sq]
	flush $e2e
	puts [analyze_end_to_end_sq]	
    }

    if { $analyse_avg_end_to_end_dist == "yes" } {
	set R [analyze_end_to_end_sq]
	set R_avg [expr $R_avg + $R]
	set R2_avg [expr $R2_avg + $R*$R]
	set l [analyze_contour_length]
	set l_avg [expr $l_avg + $l]
	set l2_avg [expr $l2_avg + $l*$l]
	incr samples
    }

    if { $vmd == "yes" } { 
	imd positions
    }

    if { $vtf == "yes" } {
	writevcf $vtfchan
    }
}

if { $analyse_avg_end_to_end_dist == "yes" } {
    set avg_e2e_fd [open $avg_end_to_end_file "a"]
    set R_mean [expr $R_avg/$samples]
    set R_sd [expr { sqrt( $R2_avg / $samples - $R_mean*$R_mean ) }]
    set l_mean [expr $l_avg/$samples]
    set l_sd [expr { sqrt( $l2_avg / $samples - $l_mean*$l_mean ) }]
    puts $avg_e2e_fd "$n_basepairs $R_mean $R_sd $l_mean $l_sd"
    close $avg_e2e_fd
}

if { $analyse_persistence_length == "yes" } {
    close $pers
}
if { $analyse_chain_parameters == "yes" } {
    close $fo
}
if { $analyse_end_to_end_dist == "yes" } {
    close $e2e     
}

