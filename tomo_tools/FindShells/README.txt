===========
Compilation
===========

  Modify 'make.sh' and run 'sh make.sh'

==============================================================================
Build the input image which is 3D, binarized (0 = void, 1 = solid), compressed
==============================================================================

  It is the role of the application 'tiffSegment' that runs like this:
    ./tiffSegment fileNameTemplate imFirst imLast GreyCut
  
  for exemple:
    ./tiffSegment ~/Desktop/SlicesY/slice%05d.tif 62 1344 30572
    
  The threshold value 'GreyCut' can be known by using imageJ/Fiji (menu Analyze/Histogram)

===============  
Find the Shells
===============

  Take the file 'InputData.txt' as an example. There are comments that may help.
  Run the finding procedure with:
    ./FindShells InputData.txt

===========
Visualizing
===========

  Run:
    ./see InputData.txt
    
  Documentation can be obtained with the touch 'k' 
