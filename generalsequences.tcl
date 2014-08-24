#
# (C) Copyright 2014, Damaris Holder, Florian Weik
#
# All rights reserved. 
#
# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#

source ./io.tcl
source ./analysis.tcl
source ./interactions.tcl

# General MD parameters
set time_step 0.1
set total_int_steps 100000000
set steps_per_loop 10000
set skin 1.0

# Langevin parameters
# kT = 0.025eV ~ 300K (k_B = 8.617e-5 eV/K)
set kT 0.025
#set kT 0.05
set gamma 0.01

# Output options
set vmd "no"
set vtf "yes"
set vtf_filename "/work/fweik/dna.vtf"

# Molecule
set n_basepairs 200
set configuration_filename "configurations/config_1000bp_31.4deg_raise4.dat"
set sequence_filename "sequences/Sequence_polyAT.dat"
# Fix one end of molecule?
set fix_lower_end "no"
# External forces on molecule
#Sheer force in +x direction
set ext_force_sheer 0.0
#Stretch force in +z direction
set ext_force_stretch 0.0

# Box geometry
# Length along the molecule
set box_z [expr 5*$n_basepairs + 500.]
# Other directions
set box_xy 250.
# Shift along molecule axis
set zshift 250.
# Shift in other directions
set center_xy [expr 0.5*$box_xy]

# Analysis
set analyse_persistence_length "yes"
set persistence_length_file "/work/fweik/peristence_length.dat"
set analyse_chain_parameters "yes"
set chain_parameter_file "/work/fweik/chain_parameters.dat"
set analyse_energy "yes"
set energy_file "/work/fweik/energy.dat"

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
set lB 561
set lambdaDB 9.6
set alpha [expr -14.23]

setup_electrostatics $lB $lambdaDB [expr 5*$lambdaDB] $kT

setup_bonded_interactions $ladderlist

if { $fix_lower_end == "yes" } {
    part 0 fix
    part 2 fix
}

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

if { $analyse_persistence_length == "yes" } {
    set pers [open $persistence_length_file "w"]
}
if { $analyse_chain_parameters == "yes" } {
    set fo [open $chain_parameter_file "w"]
}

set largest 0

puts "cell_grid [ setmd cell_grid ] cell_size [ setmd cell_size ] max_cut [ setmd max_cut ] max_range [ setmd max_range ] box_l [ setmd box_l ]"

if { $analyse_energy == "yes"} {
    set energy_fd [open $energy_file "w"]
    puts $energy_fd "#total coulomb angular basepair stacking kinetic backbone"
}

for { set i 0 } { $i <= $int_loops } { incr i } {
    puts "Loop $i of $int_loops, (time [format %.2g [expr $i*$time_step*$steps_per_loop]] of [format %.2g [expr $total_int_steps*$time_step]])."
    puts [time { integrate $steps_per_loop }]

    if { $analyse_energy == "yes" } {
	set coulomb [expr [analyze energy coulomb] + [analyze energy bonded 1]]
	set backbone [analyze energy bonded 0]
	set angular [expr [analyze energy bonded 111] + [analyze energy bonded 112] + [analyze energy bonded 121] + [analyze energy bonded 122] + [analyze energy bonded 131] + [analyze energy bonded 132] + [analyze energy bonded 141] + [analyze energy bonded 142]]

	set basepair [expr [analyze energy bonded 211] + [analyze energy bonded 212] + [analyze energy bonded 221] + [analyze energy bonded 222]]
	set total [analyze energy total]
	set stacking [expr [analyze energy bonded 312]]
	set kinetic [analyze energy kinetic]

	puts "total $total coulomb $coulomb angular $angular basepair $basepair stacking $stacking kinetic $kinetic backbone $backbone"
	puts $energy_fd "$total $coulomb $angular $basepair $stacking $kinetic $backbone"
	flush $energy_fd
    }

    if { $analyse_chain_parameters == "yes" } {
	set stacking [analyze_stacking_all]
	puts $fo "[analyze_bps] $stacking"
	flush $fo
	puts "<theta_tw> [format %.2f [lindex $stacking 0]] <rss> [format %.2f [lindex $stacking 1]]"
    }

    if { $analyse_persistence_length == "yes" } {
	puts $pers [analyze_pl]
	flush $pers
    }

    if { $vmd == "yes" } { 
	imd positions
    }

    if { $vtf == "yes" } {
	writevcf $vtfchan
    }
}

if { $analyse_energy } {
    close $energy_fd
}

if { $analyse_chain_parameters == "yes" } {
	close $fo	
}

if { $analyse_persistence_length == "yes" } {
    close $pers
}
