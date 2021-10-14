#!/bin/bash

## Nabbed from Marek Cmero: https://github.com/Oshlack/MINTIE/blob/master/install_linux64.sh
## Marek adapted this from code by Nadia Davidson: https://github.com/Oshlack/JAFFA/blob/master/install_linux64.sh
## This script installs the prerequisite software for the Slinker pipeline
## It will fetch each tool from the web and place it into the tools/ subdirectory.
## Paths to all installed tools can be found in the file workflows/tools.groovy at the
## end of execution of this script. These paths can be changed if a different
## version of software is desired.

cd tools

#a list of which programs need to be installed
commands="bpipe samtools star gffread stringtie"

#installation methods
function bpipe_install {
    wget -O bpipe-0.9.9.5.tar.gz https://github.com/ssadedin/bpipe/releases/download/0.9.9.5/bpipe-0.9.9.5.tar.gz
    tar -zxvf bpipe-0.9.9.5.tar.gz ; rm bpipe-0.9.9.5.tar.gz
    ln -s $PWD/bpipe-0.9.9.5/bin/* $PWD/bin/
}

function samtools_install {
    wget --no-check-certificate http://sourceforge.net/projects/samtools/files/samtools/1.13/samtools-1.13.tar.bz2
    tar -jxvf samtools-1.13.tar.bz2 ; rm samtools-1.13.tar.bz2
    make prefix=$PWD install -C samtools-1.13/
}

function star_install {
    wget --no-check-certificate https://github.com/alexdobin/STAR/archive/refs/tags/2.7.3a.tar.gz
    tar -jxvf 2.7.3a.tar.gz ; rm 2.7.3a.tar.gz
    mv 2.7.3a star_2.7.3a
    make prefix=$PWD -C star-2.7.3a/
}

function gffread_install {
    wget --no-check-certificate https://github.com/gpertea/gffread/archive/refs/tags/v0.12.7.zip
    unzip v0.12.7.zip ; rm v0.12.7.zip
    mv v0.12.7 gffread_0.12.7
    make prefix=$PWD -C gffread_0.12.7/
}

function stringtie_install {
    wget --no-check-certificate https://github.com/gpertea/stringtie/archive/refs/tags/v2.1.7.zip
    unzip v2.1.7.zip
    rm v2.1.7.zip
    mv v2.1.7 stringtie-2.1.7
    make prefix=$PWD -C stringtie-2.1.7/
}

echo "// Path to tools used by the Slinker pipeline" > ../workflows/tools.groovy

for c in $commands ; do
    c_path=`which $PWD/bin/$c 2>/dev/null`
    if [ -z $c_path ] ; then
    echo "$c not found, fetching it"
    ${c}_install
    c_path=`which $PWD/bin/$c 2>/dev/null`
    fi
    echo "$c=\"$c_path\"" >> ../workflows/tools.groovy
done

#loop through commands to check they are all installed
echo "Checking that all required tools were installed:"
Final_message="All commands installed successfully!"
for c in $commands ; do
    c_path=`which $PWD/bin/$c 2>/dev/null`
    if [ -z $c_path ] ; then
    echo -n "WARNING: $c could not be found!!!! "
    echo "You will need to download and install $c manually, then add its path to workflows/tools.groovy"
    Final_message="WARNING: One or more command did not install successfully. See warning messages above. \
                       You will need to correct this before running Slinker."
    else
        echo "$c looks like it has been installed"
    fi
done
echo "**********************************************************"
echo $Final_message