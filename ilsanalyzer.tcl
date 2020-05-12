# VERSION 2.0 - 05/2020

puts "\nWelcome to";
puts "-----------"
puts "................................................................................"
puts "................................................................................"
puts ".............OOO,...OOO........ZOOOOOOO........................................."
puts ".............OOO...,OO8.......=OOO.............................................."
puts ".............OOO...?OOZ........OOZ.............................................."
puts ".............OOO...OOOZ........:OOO8............................................"
puts ".............OOO...OOO=..........8OOO:.........................................."
puts ".............OOO...OOO.............OOOO........................................."
puts ".............OOO...OOO..............OOO........................................."
puts ".............OOO...OOO8OOOOOZ.O8Z,0OOOO........................................."
puts ".............OOO...OOOOOOOOOO.8OOOOOO..........................................."
puts "................................................................................"
puts "................................................................................"
puts ".......................................,........................................"
puts ".....7OOOO.............................OO,......................................"
puts ".....OO.OOO..........,,.........~......OO..........................=..........:."
puts "....OO...OO ....OO+OOOOOO...OO7,.OOO...OO..OOO....OO..OOOOOOO..,OO8,8OO...OOOOO?"
puts "...OOO...OOO....OOO....OOI.......ZOO...OO..,OO...OOI.....OOO...OO....OOO..OOO..."
puts "..OOOOOOOOOO....OO.....OO,..,OOOOOOO...OO...OO0.ZOO.....OOZ...OOOOOOOOOO.,OO...."
puts ".OOO .....OOO...OO.....OO..OOO...,OO..7OO....OO.OO....ZOO.....?OO .......7OO...."
puts "~OO:.......OO=.ZOO.....OO..OOO...ZOO..OOO....OOOO....OOOI+??0..OOO?...O..OOO...."
puts "OOO........OOO.OOO.....O8...OOOO.OOO..OOO.....OOO....OOOOOOOZ....OOOOOO..OOO...."
puts ".............................................OOO................................"
puts "..........................................OOOOO................................."
puts "................................................................................"
puts "................................................................................"
puts "....................................by...Diego J. Alonso de Armiño.......(2020)."
puts "................................................................................"
puts "...........and..Juan.P.Bustamante..Victoria.G.Dumas..&..Marcelo.A.Marti..(2012)."
puts "................................................................................"

puts "\n\n--------- a programm that analyzes ligand migration thought ILS grids ----------\n\n"

##########################################################################
#            PROCEDURE DEFINITIONS
##########################################################################
 
 
proc load_pdb_to_vmd {a_pdb_file} {
   # LOAD A PDB FILE WITH PROTEIN STRUCTURE INTO VMD
   if {[file exists $a_pdb_file] == 1} {
      axes location Off
      mol new $a_pdb_file type pdb first 0 last -1 step 1 waitfor 1
      mol modstyle 0 0 NewRibbons 0.300000 10.000000 3.000000 0
      display projection Orthographic
      #
      mol selection "resname HEM"
      mol representation Licorice
      mol color Name
      mol addrep top
      #
      mol selection "noh and resid 75 37 88 48 21"
      mol representation Licorice 0.200000 12.000000 12.000000
      mol modmaterial 2 0 AOChalky
      mol color Name
      mol addrep top
      #
      set pdb_ok 1;
   } else {
      puts "File not found. Exiting..."
      quit
   }
}

proc load_dx_to_vmd {a_dx_file} {
   # LOAD DX FILE WITH ILS RESULTS INTO VMD.
   if {[file exists $a_dx_file] == 1} {
      set dx_file [open $a_dx_file r]
      
      mol new $a_dx_file type dx first 0 last -1 step 1 waitfor 1 volsets 0
      mol modstyle 0 1 Isosurface 0.723546 0 0 0 1 1
      mol modmaterial 0 1 Transparent
      display resetview
      
      set dx_ok 1;
      set seguir 1;
   } else {
      puts "DX file not found. Exiting..."
      quit
   }
}


proc read_DX {a_dx_file} {
   # Reads ILS .dx grid file. 
   # Ouputs 
   #   * grilla(i,j,k) --> Array variable. Contains energy values as a 
   #                       function of grid coordinates.
   #   * varlist(i)    --> Array varible. Contains energy values as a 
   #                       one-dimensional list.
   #

   global grilla
   global valList
   global X
   global Y
   global Z
   global xdelta
   global ydelta
   global zdelta
   global xOrigen
   global yOrigen
   global zOrigen
   
   set in [open $a_dx_file r]

   # Reading Title
   set InputLine [gets $in]
   set InputLine [gets $in]

   # Reading maximum grid indexes.
   set InputLine [gets $in]
   scan $InputLine "object 1 class gridpositions counts %i %i %i" X Y Z

   # Reading origin of coordinates.
   set InputLine [gets $in]
   scan $InputLine "origin %e %e %e" xOrigen yOrigen zOrigen
   set origin [list $xOrigen $yOrigen $zOrigen ]

   # Reading deltas for each coordinate.
   set InputLine [gets $in]
   scan $InputLine "delta %e %e %e" xdelta dum2 dum1
   set InputLine [gets $in]
   scan $InputLine "delta %e %e %e" dum1 ydelta dum2
   set InputLine [gets $in]
   scan $InputLine "delta %e %e %e" dum2 dum1 zdelta
   set xVec [list [expr $xdelta * $X] 0 0 ]
   set yVec [list 0 [expr $ydelta * $Y] 0 ]
   set zVec [list 0 0 [expr $zdelta * $Z] ]

   # Reading inconsequential stuff.
   set InputLine [gets $in]
   set InputLine [gets $in]
   set total [expr int($X* $Y * $Z / 3)]

   # Printing summary of .DX file prologue.
   puts "Dimensions = x:$X y:$Y z:$Z"
   puts "Origin: $xOrigen $yOrigen $zOrigen" 
   puts "Resolution = x:$xdelta y:$ydelta z:$zdelta"
   puts "Reading values..."

   # Now reading energy values.
   set a 0;
   set b 1;
   set c 2;

   set count 0
   while {$count < $total} {
      set InputLine [gets $in]
      scan $InputLine "%e %e %e" v1 v2 v3
      set valList($a) $v1
      set valList($b) $v2
      set valList($c) $v3
      set a [expr $a+3]
      set b [expr $b+3]
      set c [expr $c+3]
      incr count
   }
   if {[expr $X * $Y * $Z - 3 * $total] == 2} {
      scan $InputLine "%e %e" v1 v2
      set valList($a) $v1
      set valList($b) $v2
      set a [expr $a+3]
      set b [expr $b+3]
   }
   if {[expr $X * $Y * $Z - 3 * $total] == 1} {
      scan $InputLine "%e"  v1
      set valList($a) $v1
      set a [expr $a+3]
   }

   puts "Loaded all energy values."
   
   set v 0
   for {set i 0} {$i < $X} {incr i} {
   for {set j 0} {$j < $Y} {incr j} {
   for {set k 0} {$k < $Z} {incr k} {
      # Compiling energies as a 3D array.
      set grilla($i,$j,$k) $valList($v);
      incr v;
   } } }

   puts "3D grid generated..."
   puts "Searching for energy minima..."
}

# -----------------------------------------------------------------------------
# ---------------------------- BUSQUEDA DEL MINIMO ----------------------------
# -----------------------------------------------------------------------------
 

proc optimizer {ref_coord_x ref_coord_y ref_coord_z xOrigen yOrigen zOrigen xdelta ydelta zdelta} {

   # MAIN OPTIMIZER
   # Given an initial position, finds the nearest representation of it in the ILS grid
   # and optimizes until a real minimum is found. The criteria for a minimum is that all
   # neighbouring cells have higher energy values than the cell we a are on.

   # Find the nearest grid point to a single reference coordinate.
   # This is done by Ref = (x0 + DeltaX * gridx, y0 + DeltaY * gridy, z0 + DeltaZ * gridz)
   # and finding the values of gridx, gridy and gridz by solving the linear equation.
   # gridx, gridy and gridz are real values, but  we need intger values so we round off.
   
   # Due to limitations intrisic to TCL, it is impossible to pass array variables to a
   # procedure directly (Yeah! I know, right?). There are rather cumbersome ways to do it
   # but I'll risk just making these following variables global.
   global grilla

   #set refx [ expr {int( ($ref_coord_x - $xOrigen)/$xdelta )} ]
   #set refy [ expr {int( ($ref_coord_y - $yOrigen)/$ydelta )} ]
   #set refz [ expr {int( ($ref_coord_z - $zOrigen)/$zdelta )} ]
   set refx $ref_coord_x
   set refy $ref_coord_y
   set refz $ref_coord_z
   
   set minfound 0
   if { [ info exists $grilla($refx,$refy,$refz) ] } { 
      return 
   }
   set Gmin $grilla($refx,$refy,$refz)
   set minx $refx
   set miny $refy
   set minz $refz
   set cnt 1
   while { $minfound == 0 } {
      set moved 0
      puts -nonewline "Starting step no "
      puts $cnt
      foreach zz [list -1 0 1] {
         foreach yy [list -1 0 1] {
            foreach xx [list -1 0 1] {
               set testx [expr $refx + $xx]
               set testy [expr $refy + $yy]
               set testz [expr $refz + $zz]
               if { [ info exists $grilla($refx,$refy,$refz) ] } { continue }
               set Gtest $grilla($testx,$testy,$testz)
               if { $Gtest < $Gmin } {
                  set Gmin $Gtest
                  set minx $testx
                  set miny $testy
                  set minz $testz
                  set moved 1
               }
            }
         }
      }
      if {$moved == 1} {
         set refx $minx
         set refy $miny
         set refz $minz
      } else {
         set minfound 1
         puts "Minimum found!"
      }
      set cnt [expr $cnt + 1]
   }
   return [list $refx $refy $refz $Gmin]
}    

proc global_search_min {ref_coord_list xOrigen yOrigen zOrigen xdelta ydelta zdelta radius outfile} {

   # GLOBAL OPTIMIZATION
   # Given an initial position, creates a square grid of given radius around said position. 
   # Each point in the grid is treated as a starting position for optimization using the
   # single point optimizer procedure.
   # All optmized positions are printed into a PDB file.

   
   # Due to limitations intrisic to TCL, it is impossible to pass array variables to a
   # procedure directly (Yeah! I know, right?). There are rather cumbersome ways to do it
   # but I'll risk just making these following variables global.
   global grilla

   # Find the nearest grid point to a single reference coordinate.
   # This is done by Ref = (x0 + DeltaX * gridx, y0 + DeltaY * gridy, z0 + DeltaZ * gridz)
   # and finding the values of gridx, gridy and gridz by solving the linear equation.
   # gridx, gridy and gridz are real values, but  we need intger values so we round off.
   set opt_list [list ]
   set out_fh [open $outfile w]
   foreach a $ref_coord_list {
      set refx [ expr {round( ([lindex $a 0] - $xOrigen)/$xdelta )} ]
      set refy [ expr {round( ([lindex $a 1] - $yOrigen)/$ydelta )} ]
      set refz [ expr {round( ([lindex $a 2] - $zOrigen)/$zdelta )} ]
      
      set minfound 0
      set Gmin $grilla($refx,$refy,$refz)
      set minx $refx
      set miny $refy
      set minz $refz
      set total [ expr 8*($radius)**3 ]
      
      set cnt 1
      for {set x [expr -1 * $radius]} {$x <= $radius} {incr x} {
      for {set y [expr -1 * $radius]} {$y <= $radius} {incr y} {
      for {set z [expr -1 * $radius]} {$z <= $radius} {incr z} {
         puts ""
         puts "POINT No $cnt / $total"
         set inix [ expr $refx + $x ]
         set iniy [ expr $refy + $y ]
         set iniz [ expr $refz + $z ]
   
         puts " $refx  $refy  $refz "
         puts " $inix  $iniy  $iniz "
   
         set opt [optimizer $inix $iniy $iniz $xOrigen $yOrigen $zOrigen $xdelta $ydelta $zdelta]
         # Now I have to separate the previous output into different variables because TCL is
         # that full of shit.
         set Gmin [lindex $opt 3]
         if { $Gmin < 20.00 && [ lsearch $opt_list $opt ] == -1 } {
            lappend opt_list $opt
         }
         incr cnt
      }}}
   }
   puts $opt_list
   # Print PDB file with optimized geometries.
   set leaveopen 1
   set outfile "min.pdb"
   foreach a $opt_list {
      set minx [lindex $a 0]
      set miny [lindex $a 1]
      set minz [lindex $a 2]
      set Gmin [lindex $a 3]
      puts "$minx $miny $minz $Gmin"
      print_pdb $out_fh $xOrigen $yOrigen $zOrigen $minx $miny $minz $xdelta $ydelta $zdelta $Gmin $cnt $leaveopen
   }
   close $out_fh
}    



proc print_pdb {out_fh xOrigen yOrigen zOrigen refx refy refz xdelta ydelta zdelta Gmin g leaveopen} {
   # Prints a PDB file with the optimized coordinates.
   set fmt "ATOM  %5d  N   MIN %5d    %8.3f%8.3f%8.3f %5.2f  0.00"
   set pos_x [expr $xOrigen + $refx * $xdelta];
   set pos_y [expr $yOrigen + $refy * $ydelta];
   set pos_z [expr $zOrigen + $refz * $zdelta];

   puts $out_fh [format $fmt $g $g $pos_x $pos_y $pos_z $Gmin ]
}



######################################################################################
#
#                  MINIMUM ENERGY TRANSITION PATHS SEARCH
#
######################################################################################

proc distance { x1 y1 z1 x2 y2 z2 } {
   set dist [ expr { sqrt ( ($x2 - $x1 )**2 + ($y2 - $y1)**2 + ($z2 - $z1)**2 ) } ]
   return $dist
}
 
 
proc get_position_from_index { indx delta origin } {
   set pos [ expr { $origin + $delta * double($indx) } ]
   return $pos
}

proc is_inside {border1_x border1_y border1_z border2_x border2_y border2_z point_x point_y point_z  line_x line_y line_z} {
   # create vectors
   set v1p_x [ expr { $point_x - $border1_x } ]
   set v1p_y [ expr { $point_y - $border1_y } ]
   set v1p_z [ expr { $point_z - $border1_z } ]

   set v2p_x [ expr { $point_x - $border2_x } ]
   set v2p_y [ expr { $point_y - $border2_y } ]
   set v2p_z [ expr { $point_z - $border2_z } ]

   set proj1 [ expr { $v1p_x * $line_x  + $v1p_y * $line_y + $v1p_z * $line_z } ]
   set proj2 [ expr { $v2p_x * $line_x  + $v2p_y * $line_y + $v2p_z * $line_z } ]

   if { $proj1 > 0.0 && $proj2 < 0.0 } {
      set result 1
   } else {
      set result 0
   }

   return $result

}
 

proc NEB_optimizer {ref_indx_x ref_indx_y ref_indx_z \
                    xOrigen yOrigen zOrigen \
          xdelta ydelta zdelta \
          pnt \
          border_list_x border_list_y border_list_z \
          line_x line_y line_z} {

   # NEB OPTIMIZER
   # Optimization of a single point between two planes defined by two points and a vector.
   # This procedure is part of the grid-NEB subrutine which finds the optimal path between 
   # two points in a grid.

   # Due to limitations intrisic to TCL, it is impossible to pass array variables to a
   # procedure directly (Yeah! I know, right?). There are rather cumbersome ways to do it
   # but I'll risk just making this following variable global.
   global grilla

   set refx $ref_indx_x
   set refy $ref_indx_y
   set refz $ref_indx_z
   
   set minfound 0
   set Gmin $grilla($refx,$refy,$refz)
   set minx $refx
   set miny $refy
   set minz $refz

   set border1_x [lindex $border_list_x $pnt]
   set border1_y [lindex $border_list_y $pnt]
   set border1_z [lindex $border_list_z $pnt]

   set border2_x [lindex $border_list_x [expr {$pnt + 1}]]
   set border2_y [lindex $border_list_y [expr {$pnt + 1}]]
   set border2_z [lindex $border_list_z [expr {$pnt + 1}]]

   set cnt 1
   while { $minfound == 0 } {
      set moved 0
      puts "Starting step no $cnt"
      foreach zz [list -1 0 1] {
         foreach yy [list -1 0 1] {
            foreach xx [list -1 0 1] {
               set testx [expr $refx + $xx]
               set testy [expr $refy + $yy]
               set testz [expr $refz + $zz]
               set test_pos_x [get_position_from_index  $testx $xdelta $xOrigen]
               set test_pos_y [get_position_from_index  $testy $ydelta $yOrigen]
               set test_pos_z [get_position_from_index  $testz $zdelta $zOrigen]
               set is_allowed [is_inside  $border1_x $border1_y $border1_z \
                                          $border2_x $border2_y $border2_z \
                                          $test_pos_x $test_pos_y $test_pos_z  \
                                          $line_x $line_y $line_z ]
               if {$is_allowed == 1} {
                  set Gtest $grilla($testx,$testy,$testz)
                  if { $Gtest < $Gmin } {
                     set Gmin $Gtest
                     set minx $testx
                     set miny $testy
                     set minz $testz
                     set moved 1
                  }
               }
            }
         }
      }
      if {$moved == 1} {
         set refx $minx
         set refy $miny
         set refz $minz
      } else {
         set minfound 1
         puts "Minimum found!"
      }
      set cnt [expr $cnt + 1]
   }
   return [list $refx $refy $refz $Gmin]
}    

proc grid_NEB {pos1x pos1y pos1z pos2x pos2y pos2z xOrigen yOrigen zOrigen xdelta ydelta zdelta out_fh3 out_fh2 out_fh1 ntot} {
   
   global grilla
   global connect_list
   # Find middlepoint
   # First compute difference vector between the two positions
   set posRx [expr {($pos2x - $pos1x)}] 
   set posRy [expr {($pos2y - $pos1y)}] 
   set posRz [expr {($pos2z - $pos1z)}] 


   # Find the nearest grid point to a single reference coordinate.
   # This is done by Ref = (x0 + DeltaX * gridx, y0 + DeltaY * gridy, z0 + DeltaZ * gridz)
   # and finding the values of gridx, gridy and gridz by solving the linear equation.
   # gridx, gridy and gridz are real values, but  we need intger values so we round off.
   set indx_ini_x [ expr {round( ($pos1x - $xOrigen)/$xdelta )} ]
   set indx_ini_y [ expr {round( ($pos1y - $yOrigen)/$ydelta )} ]
   set indx_ini_z [ expr {round( ($pos1z - $zOrigen)/$zdelta )} ]

   set indx_end_x [ expr {round( ($pos2x - $xOrigen)/$xdelta )} ]
   set indx_end_y [ expr {round( ($pos2y - $yOrigen)/$ydelta )} ]
   set indx_end_z [ expr {round( ($pos2z - $zOrigen)/$zdelta )} ]

   set indx_x [ list [ expr {round( ($posRx -$xOrigen)/$xdelta )} ] ]
   set indx_y [ list [ expr {round( ($posRy -$xOrigen)/$ydelta )} ] ]
   set indx_z [ list [ expr {round( ($posRz -$xOrigen)/$zdelta )} ] ]

############################################################################
#  GENERATE A LINEAR STRING OF POINTS CONNECTING INITIAL AND FINAL POSITIONS
############################################################################

   set cnt 1
   set endpoint 0
   set refx $indx_ini_x
   set refy $indx_ini_y
   set refz $indx_ini_z
   set ref_dist [ distance $pos2x $pos2y $pos2z $pos1x $pos1y $pos1z ]
   while {$endpoint == 0} {
      set elements_left 1
      set nno 1
      foreach zz [list -1 0 1] {
      foreach yy [list -1 0 1] {
      foreach xx [list -1 0 1] {
         set indx_test_x [expr $refx + $xx]
         set indx_test_y [expr $refy + $yy]
         set indx_test_z [expr $refz + $zz]
         set pos_test_x [ expr { $xOrigen + $xdelta * double($indx_test_x) } ]
         set pos_test_y [ expr { $yOrigen + $ydelta * double($indx_test_y) } ]
         set pos_test_z [ expr { $zOrigen + $zdelta * double($indx_test_z) } ]
         set dist [ distance $pos2x $pos2y $pos2z $pos_test_x $pos_test_y $pos_test_z ]
         if { $dist < $ref_dist } {
            set selected_x $indx_test_x
            set selected_y $indx_test_y
            set selected_z $indx_test_z
            set ref_dist $dist
         }
         if { $indx_test_x == $indx_end_x && \
              $indx_test_y == $indx_end_y && \
              $indx_test_z == $indx_end_z } {
            set endpoint 1
         }
         set nno [expr { $nno + 1 } ]
      }}}
      if {$endpoint == 0} {
         puts "Found point No $cnt in initial pathway."
         lappend ini_path_x $selected_x
         lappend ini_path_y $selected_y
         lappend ini_path_z $selected_z
         set refx $selected_x
         set refy $selected_y
         set refz $selected_z
      } else {
         puts "SEARCH DONE"
         puts "FOUND A TOTAL OF [expr $cnt - 1] INTERMEDIATE POINTS."
      }
      set cnt [expr { $cnt + 1 } ]
   }

############################################################################
#  RESTRAINED OPTIMIZATION
############################################################################

   set full_band_points_x [ linsert $ini_path_x 0 $indx_ini_x ]
   set full_band_points_y [ linsert $ini_path_y 0 $indx_ini_y ]
   set full_band_points_z [ linsert $ini_path_z 0 $indx_ini_z ]

   lappend full_band_points_x $indx_end_x 
   lappend full_band_points_y $indx_end_y 
   lappend full_band_points_z $indx_end_z 


   set band_length [ llength $ini_path_x ]
   set full_band_length [ llength $full_band_points_x ]

   # GENERATE POINTS MARKING THE LIMITS OF THE OPTIMIZATION SPACE FOR EACH POINT.
   set wx [ get_position_from_index [lindex $full_band_points_x 0] $xdelta $xOrigen ]
   set wy [ get_position_from_index [lindex $full_band_points_y 0] $ydelta $yOrigen ]
   set wz [ get_position_from_index [lindex $full_band_points_z 0] $zdelta $zOrigen ]

   set normR [ expr { sqrt( $posRx**2 + $posRy**2 + $posRz**2) } ]
   set normRx [ expr { $posRx / $normR } ]
   set normRy [ expr { $posRy / $normR } ]
   set normRz [ expr { $posRz / $normR } ]
  puts "STARTING BORDERS SETUP"
   for {set pnt1 0} {$pnt1 < [expr $full_band_length - 1]} {incr pnt1} {
      puts "Step $pnt1" 
      set pnt2 [ expr $pnt1 + 1 ]

      set pnt2_x [ get_position_from_index [lindex $full_band_points_x $pnt2] $xdelta $xOrigen ]
      set pnt2_y [ get_position_from_index [lindex $full_band_points_y $pnt2] $ydelta $yOrigen ]
      set pnt2_z [ get_position_from_index [lindex $full_band_points_z $pnt2] $zdelta $zOrigen ]

      set vx [ expr { $pnt2_x - $wx } ]
      set vy [ expr { $pnt2_y - $wy } ]
      set vz [ expr { $pnt2_z - $wz } ]

      #set dotprod [ expr {$vx*$posRx + $vy*$posRy + $vz*$posRz} ]
      set dotprod [ expr {$vx*$normRx + $vy*$normRy + $vz*$normRz} ]
      set projx [ expr { $dotprod*$normRx } ]
      set projy [ expr { $dotprod*$normRy } ]
      set projz [ expr { $dotprod*$normRz } ]

      lappend border_list_x [ expr { $wx + 0.5 * $projx } ]
      lappend border_list_y [ expr { $wy + 0.5 * $projy } ]
      lappend border_list_z [ expr { $wz + 0.5 * $projz } ]

      set wx [ expr { $wx + $projx } ]
      set wy [ expr { $wy + $projy } ]
      set wz [ expr { $wz + $projz } ]

   }
   # OPTIMIZE EACH POINT OF THE PATH WITHIN ITS RESTRICTIONS.
   
   lappend min_list_x $indx_ini_x
   lappend min_list_y $indx_ini_y
   lappend min_list_z $indx_ini_z
   lappend Gmin_list  1.0
   for { set pnt1 0 } { $pnt1 < $band_length } { incr pnt1 } {
      set ref_indx_x [lindex $ini_path_x $pnt1]
      set ref_indx_y [lindex $ini_path_y $pnt1]
      set ref_indx_z [lindex $ini_path_z $pnt1]
      set opt [ NEB_optimizer $ref_indx_x $ref_indx_y $ref_indx_z $xOrigen $yOrigen $zOrigen $xdelta $ydelta $zdelta $pnt1 $border_list_x $border_list_y $border_list_z $posRx $posRy $posRz]
      lappend min_list_x [lindex $opt 0]
      lappend min_list_y [lindex $opt 1]
      lappend min_list_z [lindex $opt 2]
      lappend Gmin_list  [lindex $opt 3]
   }
   lappend min_list_x $indx_end_x
   lappend min_list_y $indx_end_y
   lappend min_list_z $indx_end_z
   lappend Gmin_list  1.0

############################################################################
#  PRINT RESULTS
############################################################################
#   set out_fh1 [open "borders.pdb" w]
   set fmt "ATOM  %5d  N   MIN %5d    %8.3f%8.3f%8.3f %5.2f  0.00"
   for {set i 0} {$i < [llength $border_list_x]} {incr i} {
      set x [ lindex $border_list_x $i ]
      set y [ lindex $border_list_y $i ]
      set z [ lindex $border_list_z $i ]
      puts $out_fh1 [format $fmt 1 1 $x $y $z 1.0 ]
   }
#   close $out_fh1
#   set out_fh2 [open "initial_path.pdb" w]
   set fmt "ATOM  %5d  N   MIN %5d    %8.3f%8.3f%8.3f %5.2f  0.00"
   for {set i 0} {$i < [llength $full_band_points_x]} {incr i} {
      set x [ lindex $full_band_points_x $i ]
      set y [ lindex $full_band_points_y $i ]
      set z [ lindex $full_band_points_z $i ]
      print_pdb $out_fh2 $xOrigen $yOrigen $zOrigen $x $y $z $xdelta $ydelta $zdelta 1.0 $i 1
   }
#   close $out_fh2
#   set out_fh3 [open "optimized_path.pdb" w]
   set fmt "ATOM  %5d  N   MIN %5d    %8.3f%8.3f%8.3f %5.2f  0.00"
   set fmtcon "CONECT  %4d %4d %4d %4d "
   for {set i 0} {$i < [llength $min_list_x]} {incr i} {
      set Gmin [ lindex $Gmin_list $i ]
      set x [ lindex $min_list_x $i ]
      set y [ lindex $min_list_y $i ]
      set z [ lindex $min_list_z $i ]
      set index [expr $ntot + $i ]
      print_pdb $out_fh3 $xOrigen $yOrigen $zOrigen $x $y $z $xdelta $ydelta $zdelta $Gmin $index 1
      if { $i == 0 } { 
         #lappend connect_list [list $index [expr $index + 1]] 
         lappend connect_list [list [expr $index + 1]] 
      } elseif { $i == [expr [llength $min_list_x] - 1] } {
         #lappend connect_list [list $index [expr $index - 1]]
         lappend connect_list [list [expr $index - 1]]
      } else {
         #lappend connect_list [list $index [expr $index - 1] [expr $index + 1]]
         lappend connect_list [list [expr $index - 1] [expr $index + 1]]
      }
   }
   return [llength $min_list_x]
#   close $out_fh3


#   puts $out_fh [format $fmt 2 2 $pos2x $pos2y $pos2z 1.0 ]
#   puts $out_fh [format $fmt 3 3 $posMx $posMy $posMz 1.0 ]

   
}
