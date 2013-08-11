# --------------------------------------------------------------
# © 2011 Winged Foot Capital Research, LLC - All rights reserved
# author: Suraj Gupta <suraj@wingedfootcapital.com>
# --------------------------------------------------------------

library( "testthat" )

# convert all warnings to errors
options( warn = 2 )

# run tests
test_package( "CSSP" )
