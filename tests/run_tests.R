# --------------------------------------------------------------
# © 2011 Winged Foot Capital Research, LLC - All rights reserved
# author: Suraj Gupta <suraj@wingedfootcapital.com>
# --------------------------------------------------------------

library( "testthat" )

# convert all warnings to errors
options( warn = 1 )

# run tests
# test_package( "CSSP" )
test_check( "CSSP" )
