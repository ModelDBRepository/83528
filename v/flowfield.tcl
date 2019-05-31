# displays robot actions as flow-field in x,y, with one fixed phi; from file /tmp/obs_flowfield.dat
# belongs to simulator language file: CortexDocking_analyze.enr


set file /tmp/obs_flowfield.dat


    if {[file exists $file]} {
       set fp [open $file r]
    } else {
       puts "no input file $file"
    }

    seek $fp 0

    # read the header
    # discrete dimensions (number of data = xdim*ydim)
    gets $fp line
    set xdim [lindex $line 0]
    set ydim [lindex $line 1]

    # size of the field (xvalues=[0 .. xsize]; yvalues=[-ysize/2 .. ysize/2])
    gets $fp line
    set xsize [lindex $line 0]
    set ysize [lindex $line 1]

    # angle
    gets $fp phi

    for {set x_ct 0} {$x_ct < $xdim} {incr x_ct} {
        for {set y_ct 0} {$y_ct < $ydim} {incr y_ct} {

            gets $fp line
            set xpos($x_ct,$y_ct) [lindex $line 0]
            set ypos($x_ct,$y_ct) [lindex $line 1]
            set forw($x_ct,$y_ct) [lindex $line 2]
            set back($x_ct,$y_ct) [lindex $line 3]
            set left($x_ct,$y_ct) [lindex $line 4]
            set righ($x_ct,$y_ct) [lindex $line 5]
        }
    }

    close $fp

    puts "xdim=$xdim ydim=$ydim"
    puts "xsize=$xsize ysize=$ysize"
    puts "phi=$phi"



set mygreen #000bbb000
set shape "3 7 3"

set fieldscale 20
set arrowscale 7
set fieldheight [expr $fieldscale * $xsize + 50]
set fieldheight [lindex [split $fieldheight "."] 0]
set fieldwidth [expr $fieldscale * $ysize + 50]
set fieldwidth [lindex [split $fieldwidth "."] 0]
canvas .field -height $fieldheight -width $fieldwidth -bg white
pack .field
bind . <q> {destroy .}

# .field create rectangle 25 25 [expr $fieldwidth - 25] [expr $fieldheight - 25] -outline grey


    for {set x_ct 0} {$x_ct < $xdim} {incr x_ct} {
        for {set y_ct 0} {$y_ct < $ydim} {incr y_ct} {

            set xpos_ [expr 25 + $xpos($x_ct,$y_ct) * $fieldscale]
            set ypos_ [expr 25 + ($ysize * 0.5 + $ypos($x_ct,$y_ct)) * $fieldscale]
            set forw_ [expr $forw($x_ct,$y_ct) * $arrowscale]
            set back_ [expr $back($x_ct,$y_ct) * $arrowscale]
            set left_ [expr $left($x_ct,$y_ct) * $arrowscale]
            set righ_ [expr $righ($x_ct,$y_ct) * $arrowscale]

            # the following is for little lines which originate from little circles

            #.field create line $ypos_ $xpos_ [expr $ypos_] [expr $xpos_ - $forw_] -fill black
            #.field create line $ypos_ $xpos_ [expr $ypos_] [expr $xpos_ + $back_] -fill red
            #.field create line $ypos_ $xpos_ [expr $ypos_ + $righ_] [expr $xpos_] -fill blue
            #.field create line $ypos_ $xpos_ [expr $ypos_ - $left_] [expr $xpos_] -fill $mygreen
            #if  {$forw($x_ct,$y_ct) == 1.000000} {
            #    .field create oval [expr $ypos_ - 1] [expr $xpos_ - 1] [expr $ypos_ + 1] [expr $xpos_ + 1] -outline black -fill black
            #}
            #if  {$back($x_ct,$y_ct) == 1.000000} {
            #    .field create oval [expr $ypos_ - 1] [expr $xpos_ - 1] [expr $ypos_ + 1] [expr $xpos_ + 1] -outline red -fill red
            #}
            #if  {$righ($x_ct,$y_ct) == 1.000000} {
            #    .field create oval [expr $ypos_ - 1] [expr $xpos_ - 1] [expr $ypos_ + 1] [expr $xpos_ + 1] -outline blue -fill blue
            #}
            #if  {$left($x_ct,$y_ct) == 1.000000} {
            #    .field create oval [expr $ypos_ - 1] [expr $xpos_ - 1] [expr $ypos_ + 1] [expr $xpos_ + 1] -outline $mygreen -fill $mygreen
            #}

            if  {$forw($x_ct,$y_ct) == 1.000000} {
                .field create line [expr $ypos_ + $forw_ * sin ($phi)] [expr $xpos_ + $forw_ * cos ($phi)] \
                                   [expr $ypos_ - $forw_ * sin ($phi)] [expr $xpos_ - $forw_ * cos ($phi)] -arrow last -arrowshape $shape -fill black
            }
            if  {$back($x_ct,$y_ct) == 1.000000} {
                .field create line [expr $ypos_ - $back_ * sin ($phi)] [expr $xpos_ - $back_ * cos ($phi)] \
                                   [expr $ypos_ + $back_ * sin ($phi)] [expr $xpos_ + $back_ * cos ($phi)] -arrow last -arrowshape $shape -fill red
            }
            if  {$righ($x_ct,$y_ct) == 1.000000} {
                .field create line [expr $ypos_ - $righ_ * cos ($phi)] [expr $xpos_ + $righ_ * sin ($phi)] \
                                   [expr $ypos_ + $righ_ * cos ($phi)] [expr $xpos_ - $righ_ * sin ($phi)] -arrow last -arrowshape $shape -fill blue
            }
            if  {$left($x_ct,$y_ct) == 1.000000} {
                .field create line [expr $ypos_ + $left_ * cos ($phi)] [expr $xpos_ - $left_ * sin ($phi)] \
                                   [expr $ypos_ - $left_ * cos ($phi)] [expr $xpos_ + $left_ * sin ($phi)] -arrow last -arrowshape $shape -fill $mygreen
            }
        }
    }


update
.field postscript -file "flow.eps"
