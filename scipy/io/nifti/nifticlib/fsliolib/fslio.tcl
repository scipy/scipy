#{{{ copyright and setup 

#   FEAT TCL FSLIO wrappers
#
#   Stephen Smith, FMRIB Image Analysis Group
#
#
#   
#   The fslio.tcl file was originally part of FSL - FMRIB's Software Library
#   http://www.fmrib.ox.ac.uk/fsl
#   fslio.tcl has now been placed in the public domain.
#   
#   Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
#   Imaging of the Brain), Department of Clinical Neurology, Oxford
#   University, Oxford, UK
#   
#   

#}}}

proc imcp { args } {
    global FSLDIR

    regsub -all "\{" $args "" cleanedargs
    regsub -all "\}" $cleanedargs "" cleanedargs

    return [ exec sh -c "${FSLDIR}/bin/imcp $args" ]
}

proc imglob { args } {
    global FSLDIR

    regsub -all "\{" $args "" cleanedargs
    regsub -all "\}" $cleanedargs "" cleanedargs

    return [ exec sh -c "${FSLDIR}/bin/imglob $args" ]
}

proc imln { args } {
    global FSLDIR

    regsub -all "\{" $args "" cleanedargs
    regsub -all "\}" $cleanedargs "" cleanedargs

    return [ exec sh -c "${FSLDIR}/bin/imln $args" ]
}

proc immv { args } {
    global FSLDIR

    regsub -all "\{" $args "" cleanedargs
    regsub -all "\}" $cleanedargs "" cleanedargs

    return [ exec sh -c "${FSLDIR}/bin/immv $args" ]
}

proc imrm { args } {
    global FSLDIR

    regsub -all "\{" $args "" cleanedargs
    regsub -all "\}" $cleanedargs "" cleanedargs

    return [ exec sh -c "${FSLDIR}/bin/imrm $args" ]
}

proc imtest { args } {
    global FSLDIR

    regsub -all "\{" $args "" cleanedargs
    regsub -all "\}" $cleanedargs "" cleanedargs

    return [ exec sh -c "${FSLDIR}/bin/imtest $args" ]
}

proc remove_ext { args } {
    global FSLDIR

    regsub -all "\{" $args "" cleanedargs
    regsub -all "\}" $cleanedargs "" cleanedargs

    return [ exec sh -c "${FSLDIR}/bin/remove_ext $cleanedargs" ]
}

