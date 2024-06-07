#/bin/bash

# installs header files of the code Lohmann, N. JSON for Modern C++ (Version 3.10.4) [Computer software]. https://github.com/nlohmann

# This version is for the use in the course "advanced programming for scientific computing" by L. Formaggia, MOX, Politecnico di Milano, Italy.

PACS_DIR=../../Examples/

rm -r -f ${PACS_DIR}/include/nlohmann
mkdir ${PACS_DIR}/include/nlohmann
cp -a -v ./include/nlohmann/* ${PACS_DIR}/include/nlohmann
cp -a -v ./include/json.hpp ${PACS_DIR}/include/
