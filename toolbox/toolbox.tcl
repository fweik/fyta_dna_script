
proc init_toolbox dir {
  source "$dir/histogram.tcl"
  source "$dir/command_line_parser.tcl"
  source "$dir/md_tools.tcl"
  source "$dir/polyelectrolytes.tcl"
  global globallist
  set globallist [ list ]
}

proc set_global { name val } {
  global globallist
  global $name
  set $name $val
  if { [ lsearch $globallist $name ] == -1 } {
    lappend globallist $name
  } 
}

