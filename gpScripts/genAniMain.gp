set term png size 1280,720
set palette grey

do for [ct=0:a]{
        outfile = sprintf('tmpAni/heat%06u.jpg',ct)
        plotfile = sprintf('tmpDat/MatDat%06u.dat',ct)
        titleStr = sprintf('Iteration: %06u',ct);
        set title titleStr
        set output outfile
        plot plotfile matrix using 1:2:(-1*($3>2)) with image
}

