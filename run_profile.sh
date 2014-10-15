#!/bin/sh
python -m cProfile -o Profile.prof ./tests/test_distance.py
echo "You can see the Profile.prof by snakeviz:"
echo "snakeviz Profile.prof"
