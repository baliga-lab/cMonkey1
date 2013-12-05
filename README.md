 cMonkey biclustering R scripts, version 1.1
========
```
 Thu Dec  8 18:19:55 PST 2005
 Authors: David J. Reiss, ISB    (c) 2005-2006
          Richard Bonneau, ISB
 http://halo.systemsbiology.org/cmonkey
```

 This code is available for non-profit use only, free to be used and modified.
 If results derived using this method are published, please cite:
 Integrated biclustering of heterogeneous genome-wide datasets for the inference
         of global regulatory networks
 by David J Reiss, Nitin S Baliga, Richard Bonneau
 BMC Bioinformatics 2006, 7:280 (2 June 2006)

 Currently supported organisms: halo, hpy, ecoli, yeast.
 All data and code are provided as R-Data binary "object" files.
 The source for the functions included in these files are visible using R
     (http://www.r-project.org) and may be easily updated/overridden.
 Installation instructions below have only been tested on UNIX-based OSes (Linux, OSX).

 INSTALLATION / RUNNING:
 1. Download the relevant organism's RData file and place in a directory called "data/"
 2. Download the "cMonkey_v1.1_source.RData" and "cMonkey.R" files and place them in the 
        current directory
 3. Download (and compile if necessary) MEME and MAST
       (http://metameme.sdsc.edu/mhmm-download.html)
 4. Create links to "meme" and "mast" from within the "progs/" directory.
 5. Start R, and in the R environment:
```
 6. >   install.packages("brlr",repos="http://cran.r-project.org") ## install brlr pkg if not installed
 7. >   load("cMonkey_v1.1_source.RData")  load all cMonkey code
 8. >   run()   enter organism to run on
```

 NOTES:
 A. Ignore all warning messages printed to the screen - only worry if it dies with an error.
 B. All output is placed in the output/ORGANISM/ directory, along with PDF files containing
     bicluster plots and statistics generated during the optimization. Also generated are
     RData dump files, for re-starting a stopped/killed session and for input to run
     the Inferelator (regulatory network inference algorithm).
 C. You may kill/restart a run at any time. To completely re-start a new job, delete the
     output/ORGANISM directory.
 D. At any time you may plot a bicluster (e.g. 2) by typing (in R, after running 
     step 7 above):
```
    >   start.up(); load.latest(); plotCluster.motif( clusterStack[[2]] )
```
 E. Or, to plot them all to a PDF, instead of the "plotCluster.motif()" function, type:
```
    >   plot.clusters( clusterStack )  PDF is placed in output/ORGANISM
```

 ADVANCED USAGE: if you want to adjust the cMonkey parameters, between steps (7) and (8)
                 above, insert:
```
 a. >   start.up(); detach( params )   ## enter the organism to run on
 b. *  change the parameter that you wish to modify in the "params" list; e.g.
       > params$kmax <- 100   ## sets the number of clusters to be learned
 c. >   attach( params )
```

 OTHER ORGANISMS: if you wish to use cMonkey on your own data for another organism,
    you need to inspect the "global.data" data structure, and replace the sub-structures
    seen there. I will be happy to offer help/suggestions/advice.

 PARALELLIZATION: Additional headaches arise if you want run cMonkey parallelized.
     I have not included the necessary ingredients in this distribution, but I will, if
     anyone expresses a desire to do so. The extra steps include:
 1. installing pvm for your machine
 2. setting up pvm-related environment variables
 3. installing the rpvm and snow R libraries
 5. loading the snow-utils code (not included here)
 4. editing the test.opt.all.clusters() function to uncomment the snow-related functions

 CONTACT for assistance: David J. Reiss, ISB -- dreiss@systemsbiology.org

