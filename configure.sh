# This script compiles all Fortran programs in the VASP4CLINT repository
compiler="gfortran"
libraries="-llapack"

cd fortran
#$compiler -Wall $libraries -O1 -o analyze_scalms analyze_scalms.f90
#$compiler -Wall $libraries -O1 -o cut_unitcell cut_unitcell.f90
#$compiler -Wall $libraries -O1 -o eval_bader eval_bader.f90
#$compiler -Wall $libraries -O1 -o mlff_select mlff_select.f90
#$compiler -Wall $libraries -O1 -o modify_xdatcar modify_xdatcar.f90
$compiler  -o partial_dos partial_dos.f90
#$compiler -Wall $libraries -O1 -o rdf_pca lm_good.f rdf_pca.f90
#$compiler -Wall $libraries -O1 -o split_freq split_freq.f90
cd ..
