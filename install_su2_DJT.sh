#!/bin/bash

echo "Installing su2_DJT in" ${HOME}
pip3 install .

echo "cleaning up"
rm -r build/ su2_DJT.egg-info/