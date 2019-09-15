#!/bin/sh
	./force72 $1
	gnuplot stest.gnu
	rm stest/*.dat
