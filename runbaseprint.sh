#!/bin/sh
	./force75 $1
	gnuplot printbasis.gnu
	rm baseprints/*.dat
