#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=2


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for circular disc Koiter Steigmann
#----------------------------
cd Validation

echo "Running circular disc Koiter Steigmann validation "
mkdir RESLT
../circular_disc_to_cylinder_ks --validate > OUTPUT_circular_disc
echo "done"
echo " " >> validation.log
echo "Circular disc Koiter Steigmann validation 1" >> validation.log
echo "------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > circular_disc_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
  ../../../bin/fpdiff.py ../validata/circular_disc_ks.dat.gz   \
  circular_disc_results.dat  >> validation.log
fi

echo "Running circular disc Koiter Steigmann validation "
../circular_disc_rotate_ks > OUTPUT_circular_disc_rotate
echo "done"
echo " " >> validation.log
echo "Circular disc Koiter Steigmann validation 2" >> validation.log
echo "------------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > circular_disc_rotate_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
  ../../../bin/fpdiff.py ../validata/circular_disc_rotate_ks.dat.gz   \
  circular_disc_rotate_results.dat  >> validation.log
fi


# Append output to global validation log file
#--------------------------------------------
cat validation.log >> ../../../../validation.log


cd ..



#######################################################################


#Check that we get the correct number of OKs
# validate_ok_count will exit with status
# 0 if all tests has passed.
# 1 if some tests failed.
# 2 if there are more 'OK' than expected.
 . $OOMPH_ROOT_DIR/bin/validate_ok_count

# Never get here
 exit 10
