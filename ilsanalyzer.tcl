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
      mol selection "noh and same residue as within 5 of name FE"
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
   global endgridX
   global endgridY
   global endgridZ
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
   scan $InputLine "object 1 class gridpositions counts %i %i %i" endgridX endgridY endgridZ

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
   set xVec [list [expr $xdelta * $endgridX] 0 0 ]
   set yVec [list 0 [expr $ydelta * $endgridY] 0 ]
   set zVec [list 0 0 [expr $zdelta * $endgridZ] ]

   # Reading inconsequential stuff.
   set InputLine [gets $in]
   set InputLine [gets $in]
   set total [expr int($endgridX* $endgridY * $endgridZ / 3)]

   # Printing summary of .DX file prologue.
   puts "Dimensions = x:$endgridX y:$endgridY z:$endgridZ"
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
   if {[expr $endgridX * $endgridY * $endgridZ - 3 * $total] == 2} {
      scan $InputLine "%e %e" v1 v2
      set valList($a) $v1
      set valList($b) $v2
      set a [expr $a+3]
      set b [expr $b+3]
   }
   if {[expr $endgridX * $endgridY * $endgridZ - 3 * $total] == 1} {
      scan $InputLine "%e"  v1
      set valList($a) $v1
      set a [expr $a+3]
   }

   puts "Loaded all energy values."

   set v 0
   for {set i 0} {$i < $endgridX} {incr i} {
   for {set j 0} {$j < $endgridY} {incr j} {
   for {set k 0} {$k < $endgridZ} {incr k} {
      # Compiling energies as a 3D array.
      set grilla($i,$j,$k) $valList($v);
      incr v;
   } } }

   puts "3D grid generated..."
   puts "Searching for energy minima..."
}

## lremove - remove items from a list
# OPTS:
# -all remove all instances of each item
# -glob remove all instances matching glob pattern
# -regexp remove all instances matching regexp pattern
# ARGS: l a list to remove items from
# args items to remove (these are 'join'ed together)
##
proc lremove {args} {
    array set opts {-all 0 pattern -exact}
    while {[string match -* [lindex $args 0]]} {
        switch -glob -- [lindex $args 0] {
            -a* { set opts(-all) 1 }
            -g* { set opts(pattern) -glob }
            -r* { set opts(pattern) -regexp }
            -- { set args [lreplace $args 0 0]; break }
            default {return -code error "unknown option \"[lindex $args 0]\""}
        }
        set args [lreplace $args 0 0]
    }
    set l [lindex $args 0]
    foreach i [join [lreplace $args 0 0]] {
        if {[set ix [lsearch $opts(pattern) $l $i]] == -1} continue
        set l [lreplace $l $ix $ix]
        if {$opts(-all)} {
            while {[set ix [lsearch $opts(pattern) $l $i]] != -1} {
                set l [lreplace $l $ix $ix]
            }
        }
    }
    return $l
}



proc write_pruned_DX {p_dx_file path_pdb cutoff} {
   # Writes a pruned ILS .dx grid file for better visualization. The pruning
   # is performed with a proximity criteria with the previously computed reaction
   # paths. All grid points within a cutoff radius of the reaction path points
   # are asigned their original energy value. The rest get a very high energy value
   # (100.00 by default) so they are not visualized in VMD. This is intended to help
   # with the visualization of complex pathway systems inside proteins.

   global grilla
   global valList
   global endgridX
   global endgridY
   global endgridZ
   global xdelta
   global ydelta
   global zdelta
   global xOrigen
   global yOrigen
   global zOrigen

   set outf [open $p_dx_file w]

   # Reading Title
   puts $outf "#"
   puts $outf "#"

   # Reading maximum grid indexes.
   puts $outf  "object 1 class gridpositions counts $endgridX $endgridY $endgridZ"

   # Reading origin of coordinates.
   puts $outf  "origin $xOrigen $yOrigen $zOrigen"

   # Reading deltas for each coordinate.
   puts $outf  "delta $xdelta 0 0"
   puts $outf  "delta 0 $ydelta 0"
   puts $outf  "delta 0 0 $zdelta"
   set xVec [list [expr $xdelta * $endgridX] 0 0 ]
   set yVec [list 0 [expr $ydelta * $endgridY] 0 ]
   set zVec [list 0 0 [expr $zdelta * $endgridZ] ]

   # Reading inconsequential stuff.
   puts $outf  "object 2 class gridconnections counts $endgridX $endgridY $endgridZ"
   puts $outf  "object 3 class array double rank 0 items [expr $endgridX*$endgridY*$endgridZ]"

   puts "CREATING PRUNED .DX FILE"

   # Read PDB with reaction paths.
   set pathf [open $path_pdb r]
   set path_data [read $pathf]
   close $pathf
   set path_lines [ split $path_data "\n" ]
   set path_lines [ lrange $path_lines 0 end-1 ]
   set npoints [ llength $path_lines ]
   foreach l $path_lines {
      lappend path_x [ lindex $l 5 ]
      lappend path_y [ lindex $l 6 ]
      lappend path_z [ lindex $l 7 ]
   }

   # Now writing energy values.
   set dbcnt 1
   set cnt 1
   set fmtdx "%11.5f%11.5f%11.5f"
   set fmtdx2 "%11.5f%11.5f"
   set fmtdx1 "%11.5f"
   set gridline [ list ]
   for {set i 0} {$i < $endgridX} {incr i} {
   for {set j 0} {$j < $endgridY} {incr j} {
   for {set k 0} {$k < $endgridZ} {incr k} {
      # Compiling energies as a 3D array.
      set grid_pos_x [get_position_from_index  $i $xdelta $xOrigen]
      set grid_pos_y [get_position_from_index  $j $ydelta $yOrigen]
      set grid_pos_z [get_position_from_index  $k $zdelta $zOrigen]
      set is_near_min 0
      for { set m 0 } { $m < $npoints } { incr m } {
         set px [ lindex $path_x $m ]
         set py [ lindex $path_y $m ]
         set pz [ lindex $path_z $m ]
         set dist [ distance $px $py $pz $grid_pos_x $grid_pos_y $grid_pos_z ]
         if { $dist < $cutoff } {
            lappend gridline $grilla($i,$j,$k)
            lappend full_list $grilla($i,$j,$k)
            set dbcnt [expr $dbcnt + 1]
            set is_near_min 1
            break
         }
      }
      if { $is_near_min == 0 } {
         lappend gridline 100.00
         lappend full_list 60.00
         set dbcnt [expr $dbcnt + 1]
      }
#      puts "is_near_min = $is_near_min"
#      if { [ expr { $cnt % 3 } ] == 0 } {
#         set g0 [ lindex $gridline 0 ]
#         set g1 [ lindex $gridline 1 ]
#         set g2 [ lindex $gridline 2 ]
#         puts $outf [ format $fmtdx $g0 $g1 $g2 ]
#         set gridline [list ]
#      }
#      puts "Element Nº: $cnt"
      set cnt [expr $cnt + 1]
   }}}

   # Now reading energy values.
   set a 0;
   set b 1;
   set c 2;

   set total [expr int($endgridX* $endgridY * $endgridZ / 3)]

   set count 0
   while {$count < $total} {
      set g0 [ lindex $full_list $a ]
      set g1 [ lindex $full_list $b ]
      set g2 [ lindex $full_list $c ]
      puts $outf [ format $fmtdx $g0 $g1 $g2 ]

      set a [expr $a+3]
      set b [expr $b+3]
      set c [expr $c+3]
#      puts "Line No: $count"
      incr count
   }
   if {[expr $endgridX * $endgridY * $endgridZ - 3 * $total] == 2} {
      set g0 [ lindex $full_list $a ]
      set g1 [ lindex $full_list $b ]
      puts $outf [ format $fmtdx2 $g0 $g1 ]
   }
   if {[expr $endgridX * $endgridY * $endgridZ - 3 * $total] == 1} {
      set g0 [ lindex $full_list $a ]
      puts $outf [ format $fmtdx1 $g0 $g1 ]
   }
   close $outf

#   if { [ llength $gridline ] == 2 } {
#      set g0 [ lindex $gridline 0 ]
#      set g1 [ lindex $gridline 1 ]
#      puts $outf [ format $fmtdx2 $g0 $g1 ]
#   } elseif { [ llength $gridline ] == 1 } {
#      set g0 [ lindex $gridline 0 ]
#      puts $outf [ format $fmtdx1 $g0 ]
#   }

   puts "Pruned 3D grid generated."
   puts "dbcnt Nº: $dbcnt"
   puts "Total elements in grid: [llength $full_list]"

}



proc read_occupation_volmap_DX {a_dx_file} {
   # Reads ILS .dx grid file.
   # Ouputs
   #   * grilla(i,j,k) --> Array variable. Contains energy values as a
   #                       function of grid coordinates.
   #   * varlist(i)    --> Array varible. Contains energy values as a
   #                       one-dimensional list.
   #

   global grilla
   global valList
   global endgridX
   global endgridY
   global endgridZ
   global xdelta
   global ydelta
   global zdelta
   global xOrigen
   global yOrigen
   global zOrigen


   set in [open $a_dx_file r]

   # Reading Title
   set InputLine [gets $in]
#   set InputLine [gets $in]

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
   set xVec [list [expr $xdelta * $endgridX] 0 0 ]
   set yVec [list 0 [expr $ydelta * $endgridY] 0 ]
   set zVec [list 0 0 [expr $zdelta * $endgridZ] ]

   # Reading inconsequential stuff.
   set InputLine [gets $in]
   set InputLine [gets $in]
   set total [expr int($endgridX* $endgridY * $endgridZ / 3)]

   # Printing summary of .DX file prologue.
   puts "Dimensions = x:$endgridX y:$endgridY z:$endgridZ"
   puts "Origin: $xOrigen $yOrigen $zOrigen"
   puts "Resolution = x:$xdelta y:$ydelta z:$zdelta"
   puts "Reading values..."

   # Now reading energy values.
   set a 0;
   set b 1;
   set c 2;

   set count 0
   set sumval 0.0
   while {$count < $total} {
      set InputLine [gets $in]
      scan $InputLine "%e %e %e" v1 v2 v3
      set valList($a) $v1
      set valList($b) $v2
      set valList($c) $v3
      set a [expr $a+3]
      set b [expr $b+3]
      set c [expr $c+3]
      set sumval [ expr { $sumval + $v1 + $v2 + $v3 } ]
      incr count
   }
   if {[expr $endgridX * $endgridY * $endgridZ - 3 * $total] == 2} {
      scan $InputLine "%e %e" v1 v2
      set valList($a) $v1
      set valList($b) $v2
      set a [expr $a+3]
      set b [expr $b+3]
      set sumval [ expr { $sumval + $v1 + $v2 } ]
   }
   if {[expr $endgridX * $endgridY * $endgridZ - 3 * $total] == 1} {
      scan $InputLine "%e"  v1
      set valList($a) $v1
      set sumval [ expr { $sumval + $v1 } ]
      set a [expr $a+3]
   }

   puts "Loaded all energy values."

   set v 0
   for {set i 0} {$i < $endgridX} {incr i} {
   for {set j 0} {$j < $endgridY} {incr j} {
   for {set k 0} {$k < $endgridZ} {incr k} {
      # Compiling energies as a 3D array.
      set deltaG [ expr {-log($valList($v)/$sumval) } ]
      set grilla($i,$j,$k) $deltaG;
      incr v;
   } } }

   puts "3D grid generated..."
   puts "Searching for energy minima..."
}


proc write_DX {out_dx_file} {
   # Writes a pruned ILS .dx grid file for better visualization. The pruning
   # is performed with a proximity criteria with the previously computed reaction
   # paths. All grid points within a cutoff radius of the reaction path points
   # are asigned their original energy value. The rest get a very high energy value
   # (100.00 by default) so they are not visualized in VMD. This is intended to help
   # with the visualization of complex pathway systems inside proteins.

   global grid
   global nx
   global ny
   global nz
   global deltax
   global deltay
   global deltaz
   global originx
   global originy
   global originz

   set outf [open $out_dx_file w]

   # Reading Title
   puts $outf "#"
   puts $outf "#"

   # Reading maximum grid indexes.
   puts $outf  "object 1 class gridpositions counts $nx $ny $nz"

   # Reading origin of coordinates.
   puts $outf  "origin $originx $originy $originz"

   # Reading deltas for each coordinate.
   puts $outf  "delta $deltax 0 0"
   puts $outf  "delta 0 $deltay 0"
   puts $outf  "delta 0 0 $deltaz"
   set xVec [list [expr $deltax * $nx] 0 0 ]
   set yVec [list 0 [expr $deltay * $ny] 0 ]
   set zVec [list 0 0 [expr $deltaz * $nz] ]

   # Reading inconsequential stuff.
   puts $outf  "object 2 class gridconnections counts $nx $ny $nz"
   puts $outf  "object 3 class array double rank 0 items [expr $nx*$ny*$nz]"

   puts "CREATING PRUNED .DX FILE"

   # Now writing energy values.
   set dbcnt 1
   set cnt 1
   set fmtdx "%11.5f%11.5f%11.5f"
   set fmtdx2 "%11.5f%11.5f"
   set fmtdx1 "%11.5f"
   set gridline [ list ]
   for {set i 0} {$i < $nx} {incr i} {
   for {set j 0} {$j < $ny} {incr j} {
   for {set k 0} {$k < $nz} {incr k} {
      # Compiling energies as a list.
      lappend full_list $grid($i,$j,$k)
      set cnt [expr $cnt + 1]
   }}}

   # Now reading energy values.
   set a 0;
   set b 1;
   set c 2;

   set total [expr int($nx* $ny * $nz / 3)]

   set count 0
   while {$count < $total} {
      set g0 [ lindex $full_list $a ]
      set g1 [ lindex $full_list $b ]
      set g2 [ lindex $full_list $c ]
      puts $outf [ format $fmtdx $g0 $g1 $g2 ]

      set a [expr $a+3]
      set b [expr $b+3]
      set c [expr $c+3]
      incr count
   }
   if {[expr $nx * $ny * $nz - 3 * $total] == 2} {
      set g0 [ lindex $full_list $a ]
      set g1 [ lindex $full_list $b ]
      puts $outf [ format $fmtdx2 $g0 $g1 ]
   }
   if {[expr $nx * $ny * $nz - 3 * $total] == 1} {
      set g0 [ lindex $full_list $a ]
      puts $outf [ format $fmtdx1 $g0 $g1 ]
   }
   close $outf

   puts "DX 3D grid generated."
   puts "dbcnt Nº: $dbcnt"
   puts "Total elements in grid: [llength $full_list]"

}

##########################################################################################
##                  GAUSSIAN KERNEL DENSITY ESTIMATOR ANALYSIS                             ##
##########################################################################################
proc KDEanalysis { pdbin dxout minmax delta sigma } {
   global grid
   global nx
   global ny
   global nz
   global deltax
   global deltay
   global deltaz
   global originx
   global originy
   global originz

   # PREPARING GRID
   puts "DEFINING GRID PARAMETERS"
   set deltax [lindex $delta 0]
   set deltay [lindex $delta 1]
   set deltaz [lindex $delta 2]

   set minx [expr [ lindex $minmax 0 0 ] - $deltax ]
   set miny [expr [ lindex $minmax 0 1 ] - $deltay ]
   set minz [expr [ lindex $minmax 0 2 ] - $deltaz ]

   set origin [ list $minx $miny $minz ]
   set originx $minx
   set originy $miny
   set originz $minz

   set maxx [ expr [lindex $minmax 1 0] + $deltax ]
   set maxy [ expr [lindex $minmax 1 1] + $deltay ]
   set maxz [ expr [lindex $minmax 1 2] + $deltaz ]

   set nx [ expr int((ceil($maxx) - floor($minx))/$deltax) ]
   set ny [ expr int((ceil($maxy) - floor($miny))/$deltay) ]
   set nz [ expr int((ceil($maxz) - floor($minz))/$deltaz) ]

   set totgp [ expr $nx*$ny*$nz ]

   # INITIALIZE GRID
   puts "INITIALIZING GRID"
   for {set i 0} {$i < $nx} {incr i} {
   for {set j 0} {$j < $ny} {incr j} {
   for {set k 0} {$k < $nz} {incr k} {
      set grid($i,$j,$k) 0.0
   }}}

   # OPEN COORDINATE FILE
   set in [open $pdbin r]

   set PI 3.1415926535897931
   set PI2 [ expr { 2.0*$PI } ]
   set sqrtPI2 [ expr { sqrt($PI2) } ]
   set Gnorm [ expr { 1.0/($sqrtPI2**3) } ]
   set Gnorm [ expr { $Gnorm/($sigma**3) } ]

   puts "STARTING MAIN LOOP"
   set counter 0
   set n_gaus 0
   while { [gets $in InputLine] >= 0 } {

      # READ LINE OFF COORDINATE FILE.
      scan $InputLine "%e %e %e" XX YY ZZ
      lappend lxx $XX
      lappend lyy $YY
      lappend lzz $ZZ
   }

   for {set i 0} {$i < [llength $lxx]} {incr i} {
      set XX [ lindex $lxx $i ]
      set YY [ lindex $lyy $i ]
      set ZZ [ lindex $lzz $i ]

      # DETERMINE THE CLOSEST GRID POINT TO THIS POINT.
      set gridx [ expr {round( ($XX - [lindex $origin 0])/$deltax )} ]
      set gridy [ expr {round( ($YY - [lindex $origin 1])/$deltay )} ]
      set gridz [ expr {round( ($ZZ - [lindex $origin 2])/$deltaz )} ]

      # DETERMINE GRID POINTS TO ACCUMULATE (CUTOFF USED).
      set strtx [ expr $gridx - 5 ]
      set strty [ expr $gridy - 5 ]
      set strtz [ expr $gridz - 5 ]
      set endx [ expr $gridx + 6 ]
      set endy [ expr $gridy + 6 ]
      set endz [ expr $gridz + 6 ]

      # COMPUTE GAUSSIAN VALUES ON AFFECTED GRID POINTS.
      for {set i $strtx} {$i < $endx} {incr i} {
      for {set j $strty} {$j < $endy} {incr j} {
      for {set k $strtz} {$k < $endz} {incr k} {
         # Check we are inside grid.
         if { $i < 0 || $i >= $nx } continue
         if { $j < 0 || $j >= $ny } continue
         if { $k < 0 || $k >= $nz } continue
         # Get cartesians corresponging to this grid point.
         set rx [ get_position_from_index $i $deltax $originx  ]
         set ry [ get_position_from_index $j $deltay $originy  ]
         set rz [ get_position_from_index $k $deltaz $originz  ]
         # Compute gaussian kernel.
         set gausx [ expr { exp(-0.5 * (($rx-$XX)/$sigma)**2) } ]
         set gausy [ expr { exp(-0.5 * (($ry-$YY)/$sigma)**2) } ]
         set gausz [ expr { exp(-0.5 * (($rz-$ZZ)/$sigma)**2) } ]
         # Acumulate on grid point.
         set grid($i,$j,$k) [ expr { $grid($i,$j,$k) + $gausx*$gausy*$gausz*$Gnorm } ]
      }}}
      set n_gaus [ expr { $n_gaus + 1 } ]
      set counter [ expr $counter + 1 ]
      puts "STEP No $counter"
   }

   # NORMALIZE GRID AND COMPUTE FREE ENERGY.
  # puts "NORMALIZING GRID AND COMPUTING FREE ENERGIES"
  # for {set i 0} {$i < $nx} {incr i} {
  # for {set j 0} {$j < $nx} {incr j} {
  # for {set k 0} {$k < $nx} {incr k} {
  #    if { $grid($i,$j,$k) > 1e-10 } {
  #       set grid($i,$j,$j) [ expr { -log($grid($i,$j,$k)/$n_gaus) } ]
  #    } else {
  #       set grid($i,$j,$k) 1000.0
  #    }
  # }}}

   puts $origin
   puts "WRITING RESULTS TO .DX FILE"
   write_DX $dxout


}


##########################################################################################
##                             MINIMA SEARCH                                            ##
##########################################################################################

proc optimizer {ref_coord_x ref_coord_y ref_coord_z xOrigen yOrigen zOrigen xdelta ydelta zdelta endgridX endgridY endgridZ } {

   # MAIN OPTIMIZER
   # Given an initial position, finds the nearest representation of it in the ILS grid
   # and optimizes until a real minimum is found. The criteria for a minimum is that all
   # neighbouring cells have higher energy values than the cell we a are on.

   # Find the nearest grid point to a single reference coordinate.
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
#   puts "$refx $refy $refz exists? ->> [ info exists $grilla($refx,$refy,$refz) ]"
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
      puts "Starting step no $cnt"
      foreach zz [list -1 0 1] {
         foreach yy [list -1 0 1] {
            foreach xx [list -1 0 1] {

               set testx [expr $refx + $xx]
               set testy [expr $refy + $yy]
               set testz [expr $refz + $zz]

               if {$testx < 0 || ($testx >= $endgridX)} {puts "Point out of grid in X. Ignoring.";continue}
               if {$testy < 0 || ($testy >= $endgridY)} {puts "Point out of grid in Y. Ignoring.";continue}
               if {$testz < 0 || ($testz >= $endgridZ)} {puts "Point out of grid in Z. Ignoring.";continue}

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

proc global_search_min {ref_coord_list xOrigen yOrigen zOrigen xdelta ydelta zdelta endgridX endgridY endgridZ radius outfile} {

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
      set total [ expr int((2*$radius+1)**3) ]

      set cnt 1
      for {set x [expr -1 * round($radius)]} {$x <= round($radius)} {incr x} {
      for {set y [expr -1 * round($radius)]} {$y <= round($radius)} {incr y} {
      for {set z [expr -1 * round($radius)]} {$z <= round($radius)} {incr z} {
         puts ""
         puts "POINT No $cnt / $total"
         set inix [ expr $refx + round($x) ]
         set iniy [ expr $refy + round($y) ]
         set iniz [ expr $refz + round($z) ]

#         puts " $refx  $refy  $refz "
#         puts " $inix  $iniy  $iniz "

         if { $inix < 0 || $inix > $endgridX } { continue }
         if { $iniy < 0 || $iniy > $endgridY } { continue }
         if { $iniz < 0 || $iniz > $endgridZ } { continue }

         set opt [optimizer $inix $iniy $iniz $xOrigen $yOrigen $zOrigen $xdelta $ydelta $zdelta $endgridX $endgridY $endgridZ]
         set Gmin [lindex $opt 3]
         if { $Gmin < 20.00 && [ lsearch $opt_list $opt ] == -1 } {
            lappend opt_list $opt
         }
	 set cnt [expr $cnt + 1]
      }}}
   }
#   puts $opt_list
   # Print PDB file with optimized geometries.
   set leaveopen 1
   set outfile "min.pdb"
   foreach a $opt_list {
      set minx [lindex $a 0]
      set miny [lindex $a 1]
      set minz [lindex $a 2]
      set Gmin [lindex $a 3]
#      puts "$minx $miny $minz $Gmin"
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

proc grid_NEB {pos1x pos1y pos1z pos2x pos2y pos2z xOrigen yOrigen zOrigen xdelta ydelta zdelta out_fh3 out_fh2 out_fh1 out_netw out_kmodel ntot ini end ene_ini ene_end} {

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
   #lappend Gmin_list  9.0
   lappend Gmin_list  $ene_ini
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
   lappend Gmin_list  $ene_end

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
   set enefh [open "energy_profile_${ini}_${end}.dat" w]
   set fmt "ATOM  %5d  N   MIN %5d    %8.3f%8.3f%8.3f %5.2f  0.00"
   set fmtene "%5d %5d  %8.3f"
   set fmtcon "CONECT  %4d %4d %4d %4d "
   for {set i 0} {$i < [llength $min_list_x]} {incr i} {
      set Gmin [ lindex $Gmin_list $i ]
      set x [ lindex $min_list_x $i ]
      set y [ lindex $min_list_y $i ]
      set z [ lindex $min_list_z $i ]
      set index [expr $ntot + $i ]
      print_pdb $out_fh3 $xOrigen $yOrigen $zOrigen $x $y $z $xdelta $ydelta $zdelta $Gmin $index 1
      # Print energy profiles
      if { $i == 0 } {
         puts $enefh [ format $fmtene $ini $end $ene_ini ]
      } elseif { $i == [expr [llength $min_list_x]-1] } {
         puts $enefh [ format $fmtene $ini $end $ene_end ]
      } else {
         puts $enefh [ format $fmtene $ini $end $Gmin ]
      }

      # Add points to conectivity list.
      if { $i == 0 } {
         lappend connect_list [list [expr $index + 1]]
      } elseif { $i == [expr [llength $min_list_x] - 1] } {
         lappend connect_list [list [expr $index - 1]]
      } else {
         lappend connect_list [list [expr $index - 1] [expr $index + 1]]
      }
   }
   close $enefh

   # Computing transition state energies, activation energies, 
   # kinetic constants and equilibrium constants.

   # Transition state energy
   puts "LIST OF PATHWAY ENERGIES"
   puts $Gmin_list
   set Gts [tcl::mathfunc::max {*}$Gmin_list] 
   puts "TRANSITION STATE ENERGY"
   puts $Gts

   # Activation energies
   set Gact1 [ expr {$Gts - $ene_ini} ]
   set Gact2 [ expr {$Gts - $ene_end} ]

   # Equilibrium constants
   set Keq1 [ expr {exp(-( $ene_end - $ene_ini )) } ]
   set Keq2 [ expr {exp(-( $ene_ini - $ene_end )) } ]

   # Kinetic constants at T = 283K, 298K and 323K
   # DEBUG
#   set k1_283 [ expr { 2.087e10 * 283.0 * exp(-$Gact1) } ]
   set k1_283 [ expr { 2.087e10 * 283.0 * exp(-$Gact1/0.562102) } ]
   set k1_298 [ expr { 2.087e10 * 298.0 * exp(-$Gact1) } ]
   set k1_323 [ expr { 2.087e10 * 323.0 * exp(-$Gact1) } ]

   # DEBUG
   #set k2_283 [ expr { 2.087e10 * 283.0 * exp(-$Gact2) } ]
   set k2_283 [ expr { 2.087e10 * 283.0 * exp(-$Gact2/0.562102) } ]
   set k2_298 [ expr { 2.087e10 * 298.0 * exp(-$Gact2) } ]
   set k2_323 [ expr { 2.087e10 * 323.0 * exp(-$Gact2) } ]

   ## FORMAT index1 index2 G1 G2 G= Gact1 Gact2 K1 K2 k1_283 k1_298 k1_323 k2_283 k2_298 k2_323
   set fmtnetw "%5d %5d  %8.3f %8.3f %8.3f %8.3f %8.3f %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e "
   set fmtkmodel "c_%s -> c_%s; k_%s * c_%s - k_%s * c_%s; k_%s = %.3e; k_%s = %.3e;"
   puts $out_netw [ format $fmtnetw $ini $end $ene_ini $ene_end $Gts $Gact1 $Gact2 $Keq1 $Keq2 $k1_283 $k1_298 $k1_323 $k2_283 $k2_298 $k2_323]

   # Printing Tellurium-compatible Kinnetic Model
   puts $out_kmodel [ format $fmtkmodel $ini $end "${ini}_${end}" $ini "${end}_${ini}" $end "${ini}_${end}" $k1_283 "${end}_${ini}" $k2_283] 
   return [llength $min_list_x]


#   puts $out_fh [format $fmt 2 2 $pos2x $pos2y $pos2z 1.0 ]
#   puts $out_fh [format $fmt 3 3 $posMx $posMy $posMz 1.0 ]


}

proc search_pathways { start_indx_list opt_file pathway_file external_indx_list} {
# MAIN PROCEDURE TO SEARCH FOR PATHWAYS BETWEEN MINIMA.
# USE THIS PROCEDURE AS A COMMAND.

   global xOrigen
   global yOrigen
   global zOrigen
   global xdelta
   global ydelta
   global zdelta
   global minlist
   global connect_list

   # FROM INITIAL INDEX LIST OBTAIN THE LIST OF RELEVANT MINIMA (MINLIST).
   set minlist [list ]
   foreach a $start_indx_list {
      foreach b $a {
         if { [ lsearch $minlist $b ] < 0 } {
            lappend minlist $b }
      }
   }

   mol new $opt_file type {pdb} first 0 last -1 step 1 waitfor 1

   set out_fh1 [open "wat_kde_borders.pdb" w]
   set out_fh2 [open "wat_kde_initial_path.pdb" w]
   set out_fh3 [open $pathway_file w]
   set out_netw [open "network.dat" w]
   set out_kmodel [open "kinetic-model.tellurium" w]
   puts $out_netw "index1 index2 G1 G2 G= Gact1 Gact2 K1 K2 k1_283 k1_298 k1_323 k2_283 k2_298 k2_323"

   set cnt 1
   set ntot 0
   set connect_list [list ]
   set cnt 1
   set c_list [list ]

   foreach a $start_indx_list {
      puts "STARTING NEB ANALYSIS No $cnt/[llength $start_indx_list]"
      set ini [lindex $a 0]
      set end [lindex $a 1]
      lappend c_list $ini
      lappend c_list $end
      set sel_ini [ atomselect top "index $ini" ]
      set sel_end [ atomselect top "index $end" ]
      set ene_ini [ $sel_ini get occupancy ]
      set ene_end [ $sel_end get occupancy ]
      set pos1 [ $sel_ini get {x y z} ]
      set pos2 [ $sel_end get {x y z} ]
      set pos1x [lindex $pos1 0 0]
      set pos1y [lindex $pos1 0 1]
      set pos1z [lindex $pos1 0 2]
      set pos2x [lindex $pos2 0 0]
      set pos2y [lindex $pos2 0 1]
      set pos2z [lindex $pos2 0 2]
      set npnts [ grid_NEB $pos1x $pos1y $pos1z  $pos2x $pos2y $pos2z \
                           $xOrigen $yOrigen $zOrigen \
                           $xdelta $ydelta $zdelta \
                           $out_fh3 $out_fh2 $out_fh1 $out_netw $out_kmodel\
                           $ntot $ini $end $ene_ini $ene_end ]
      set ntot [expr { $ntot + $npnts } ]
      set cnt [expr $cnt + 1]
   }
   # Print external terms to kinetic model
   set fmtkmodel "c_%s -> CO + Hb; k_off * c_%s - k_on * CO * Hb;"
   foreach ex $external_indx_list {
	   puts $out_kmodel [format $fmtkmodel $ex $ex]
   }
   puts $out_kmodel ""
   puts $out_kmodel "k_on = 2.5e5; k_off = 1e10"
   puts $out_kmodel ""

   # Print initial concentrations to kinetic model
   set unique_c_list [lsort -unique $c_list]
   puts $unique_c_list
   foreach b $unique_c_list {
	   puts $out_kmodel [format "c_%s = %s" $b "0.0"]
   }

   close $out_fh3
   close $out_fh2
   close $out_fh1
   close $out_netw
   close $out_kmodel
}

proc kde_analysis_preparation { selection batch max_frames top traj output do_build_input } {

   global kde_line_cnt

   set fmt "%11.5f%11.5f%11.5f"
   if { $do_build_input == 1 } {
      set fh [open "tmp1" w]
      set kde_line_cnt 1

      set frm 1
      exec rm -rf $output
      for {set i 0} {$i < $max_frames} {incr i $batch} {
         mol load parm7 $top
         mol addfile $traj step 1 first $i last [expr $i + $batch - 1]  waitfor all type netcdf

         set sel [ atomselect top $selection ]
         set nfr [molinfo top get numframes]

         for {set j 0} {$j < $nfr} {incr j} {
            $sel frame $j
            $sel update
            set xyz_list [ $sel get {x y z} ]
            foreach a $xyz_list {
               set x [ lindex $a 0 ]
               set y [ lindex $a 1 ]
               set z [ lindex $a 2 ]
               puts $fh [ format $fmt $x $y $z ]
               set kde_line_cnt [expr $kde_line_cnt + 1]
            }
            $sel writepdb tmp.pdb
            exec cat tmp.pdb >> $output
            puts "$frm / $max_frames"
            set frm [ expr $frm + 1 ]
         }
         mol delete top
         $sel delete
      }
      close $fh
   }

   exec  sed -i /END/d $output
   exec  sed -i /CRYST/d $output

   mol load pdb $output
   set selmm [ atomselect top "all" ]
   set minmax [ measure minmax $selmm ]

   set fh [open "tmp2" w]
   puts $fh [format $fmt [lindex $minmax 0 0] [lindex $minmax 0 1] [lindex $minmax 0 2] ]
   puts $fh [format $fmt [lindex $minmax 1 0] [lindex $minmax 1 1] [lindex $minmax 1 2] ]
   close $fh

   exec rm -rf "${output}.xyz"
   exec cat tmp2 tmp1 > "${output}.xyz"

}
