#!/bin/bash

python3 setup.py clean
python3 setup.py build bdist_wheel
pip3 install --upgrade --no-deps --force-reinstall dist/open_darts_flash*.whl
