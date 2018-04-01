# Additional Information for Developers

## Building a Distribution

When you start from the actual repository source and not a
distribution you can use

        make -f Makefile.build

to run the necessary libtool, autoconf, etc routines and copy a few
additional M4 scripts to the right place to put the directory in a
similar state to that of a distribution directory. Then you can run
the usual

        ./configure; make; sudo make install

If the state of the directory is messed up, then you can run

        make -f Makefile.build clean

to do a brutal cleanup of everything. Yes, there are various clean
commands in Makefile, but this seems more robust and convenient when
developing. After this command `git status` should not list any magic
files that you did not write. Finally, you can use

        make -f Makefile.build dist

to build a distribution `vec-<version>.tar.gz` in a single
command. This merely sets up things and then runs the usual `make
dist`.
