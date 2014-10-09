#!/bin/sh
pyside_postinstall.py -install >/dev/null
py.test
