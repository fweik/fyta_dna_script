
global novalue invalidvalue
#global novalue invalidvalue
set novalue "_no _value"
set invalidvalue "_invalid _value"

proc parse_string_for_variable { the_string varname } {
  global novalue 
## first check if there is a string in double quotes (") assigned to the variable
  set match ""
  set re "$varname\[\[:space:\]\]*=\[\[:space:\]\]*\"(.*)\""
  regexp -expanded $re $the_string match sub
  if { $match != ""} {
    return $sub
  } 
## then check if there is a string in single quotes (') assigned to the variable
  set re "$varname \[\[:space:\]\]* = \[\[:space:\]\]*'(.*)'"
  regexp -expanded $re $the_string match sub
  if { $match != ""} {
    return $sub
  } 
## If both is not the case, read all characters till the next whitespace
  set re "$varname\[\[:space:\]\]*=\[\[:space:\]\]*(\[^\[:space:\]\]*)"
  regexp -expanded $re $the_string match sub
  if { $match != ""} {
    return $sub
  } 
  return $novalue 
}


## Now include typechecking!
proc parse_string_for_int_variable { the_string varname } {
  global novalue invalidvalue
  set varvalue [ parse_string_for_variable $the_string $varname ]
  if { $varvalue == $novalue } {
    return $novalue
  }
  set match ""
  regexp "\[+-\]?\[0-9\]+" $varvalue match 
  if { $match == $varvalue } {
    return $match
  }
  return $invalidvalue
}

proc parse_string_for_float_variable { the_string varname } {
  global novalue
  global invalidvalue 
  set varvalue [ parse_string_for_variable $the_string $varname ]
  if { $varvalue == $novalue } {
    return $novalue
  }
  set match ""
  regexp "\[+-\]?\[0-9\]+\.?\[0-9\]*(\[eE\]\[+-\]?\[0-9\]+)*" $varvalue match 
  if { $match == $varvalue } {
    return $match
  }
  return $invalidvalue
}


### Now a function that can be conviently used to read input from the
### command line, and assigns the given default value to the matching 
### variable if it is not given
proc get_command_line_argument {argname default_value } {
  global argc
  global argv 
  global novalue 
  global invalidvalue 
  set argvalue [ parse_string_for_variable $argv $argname ]
  if { $argvalue == $invalidvalue } {
    puts "Invalid Argument $argname"
    return $invalidvalue
  }
  if { $argvalue == $novalue } {
    return $default_value
  }
  return $argvalue
}


## The same, but with typechecking
proc get_command_line_int_argument { argname { default_value $novalue } } {
  global argc 
  global argv
  global novalue 
  global invalidvalue 
  set argvalue [ parse_string_for_int_variable $argv $argname ]
  puts $argvalue
  if { $argvalue == $invalidvalue } {
    puts "Invalid Argument $argname"
    return $invalidvalue
  }
  if { $argvalue == $novalue } {
    return $default_value
  }
  return $argvalue
}

proc get_command_line_float_argument { argname default_value } {
  global argc
  global argv 
  global novalue 
  global invalidvalue 
  set argvalue [ parse_string_for_float_variable $argv $argname ]
  if { $argvalue == $invalidvalue } {
    puts "Invalid Argument $argname"
    return $invalidvalue
  }
  if { $argvalue == $novalue } {
    return $default_value
  }
  return $argvalue
}

proc save_paramfile { fname listofparamnames } {
  set listofparamvalues [ list ]
  foreach param $listofparamnames {
    global $param 
    lappend listofparamvalues [ set $param ]
  } 
  set ofile [ open $fname "w" ]
  for { set i 0 } {$i<[llength $listofparamnames]} { incr i } {
    puts $ofile "[ lindex $listofparamnames $i ] = [ lindex $listofparamvalues $i ]"
  }
  close $ofile
}


proc load_paramfile { fname listofparamnames } {
  set ifile [ open $fname "r" ]
  set data [ read $ifile ]

  foreach param $listofparamnames {
    global $param 
    set_global $param [ parse_string_for_variable $data $param ]
    puts "set_global $param [ lindex [ parse_string_for_variable $data $param ] 0 ]"
  } 
  close $ifile
}

