#!/bin/bash
# Script to coordinate the collection of variants from a TVC analysis, naming and concatenation
# of the .txt files containing the variants, and running of cpscChecker.pl to analyze the 
# CPSC sample performance.
#
# Updated for the 1066 R&D server in order to allow for flexible calling.  Won't run the auto
# cpscChecker if 'is_RandD_server' set to 0.  Can feed it a custom sampleKey file if running on 
# 1066 server (seems the index is different than on the ATRF server for some reason.
# 
# Created 3/2/13 - Dave Sims
#
##################################################################################################

#set -x

VERSION="$(basename $0) - v0.8.2"
USAGE="$(cat <<EOT
$VERSION [options]

Program to collect variant call files from an Ion Torrent run and run the cpscChecker utility.  This 
program needs to be run from an experiment results directory, and must be run after the Torrrent Variant
Caller (TVC) plugin has been run.  

    -l    Non-default lookup file for the cpscChecker utility (default: 3)
    -s    Custom sampleKey file rather than the default based on the ion_params_00.json file
    -v    Version information
    -h    This help text

EOT
)"

TVCout="plugin_out/variantCaller_out"
resultsDir="$(pwd)"
colvarsDir="$resultsDir/collectedVariants"
runName="$(echo "$resultsDir" | perl -pe 's/.*user_([P|M]CC-\d+)-\d+.*/$1/')"
cpscSample="IonXpress_001.txt"
cpsc_lookup=3 #Default lookup file for cpscChecker (see cpscChecker for inforamation)
is_RandD_server=0 #Change to 0 for production server with locked pipeline.

while getopts :hvs:l: opt;
do
	case $opt in
		v)
			printf "$VERSION\n"
			exit 0;
			;;
		h)
			printf "$USAGE\n"
			exit 0
			;;
		s)
			customKey="$OPTARG"
			if [[ ! -f $customKey ]]; then
				printf "The custom sampleKey file '$customKey' does not exist\n"
				exit 1
			fi
			;;
		l)
			cpsc_lookup="$OPTARG"
			;;
		:)
			echo "Option -$OPTARG requires an argument!\n"
			printf "$USAGE\n\n"
			exit 1
			;;
		\?)
			printf "Invalid option: -$OPTARG\n"
			printf "$USAGE\n\n"
			exit 1
			;;
	esac
done
shift $((OPTIND - 1))

# Check to make sure TVC has been run and you are in the results dir
if [ ! -d "$TVCout" ]; then
	printf "ERROR: Either you are not running this script from a run results folder, or TVC has not been run on this sample.\n"
	exit 1
fi

# Check to see if varCollector has already been run
if [ -d "$colvarsDir" ]; then
	printf "WARNING: 'collectedVariants' directory found.  Continuing will overwrite results.  Continue [y/n]? "
	read -a overwrite
	case "$overwrite" in
		y|Y) echo "Overwriting old results." && rm -r "$colvarsDir";;
		n|N|"") echo "Exiting..." && exit 1;;
		*) echo "Not a valid choice!" && exit 1;;
	esac
fi

# Make collectedVariants directory and run sampleKeyGen.pl to generate a sample key for the run
mkdir $colvarsDir;
# Use custom key instead of generating one and copy to $colvarsDir
if [[ $customKey ]]; then
	cp $customKey $colvarsDir/
else
	eval "sampleKeyGen -o $colvarsDir/sampleKey.txt /opt/mocha/varCollector/resources/bcIndex.txt ion_params_00.json" 
	if [[ $? -ne 0 ]]; then
		echo "The sampleKeyGen tool had a problem and failed to run correctly!" 
		exit 1
	fi
fi

# Run collectVariants.pl
printf "Running collectVariants.pl on TVC results.  See output in 'collectedVariants'\n"
# TODO: use relative package path.
eval "perl /opt/mocha/varCollector/collectVariants.pl" 
if [[ $? -ne 0 ]]; then
	echo "The collectVariants.pl script had problems and  failed to run!" 
	exit 1
fi

if [ ! -d "$colvarsDir" ]; then 
	printf "ERROR: no 'collectedVariants' directory found.  Something bad happened!\n"
	exit 1
fi

if [[ $is_RandD_server -eq 1 ]]; then
	printf "Running in an R&D environment.  No automatic cpscChecker run.  Please run the cpscChecker tool manually\n"
else
	# Verify conditions are correct for cpscChecker 
	if [ ! -f "$colvarsDir/sampleKey.txt" ]; then
		printf "ERROR: No sampleKey file found.  Can not run cpscChecker!\n"
	elif [ ! -f "$colvarsDir/IonXpress_001.txt" ]; then
		printf "ERROR: no "$cpscSample" file found.  Did you run a CPSC sample?\n";
		exit 1
	fi	
		
	#  If conditions are OK, let's run the script.
	cd $colvarsDir
	printf "Running cpscChecker for '$runName' on sample '$cpscSample'...\n"
	eval "cpscChecker -l $cpsc_lookup $cpscSample"
fi

# Prepare a tarball of results for experiment folder 
cd $colvarsDir
tar vcfz "$runName"_variants.tar.gz * > /dev/null 2>&1

printf "\nScript complete.  Check results in 'collectedVariants' directory\n"
