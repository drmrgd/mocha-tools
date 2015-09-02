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
VERSION="$(basename $0) - v2.3.0"
USAGE="$(cat <<EOT
$VERSION [options]

Program to collect variant call files from an Ion Torrent run and run the cpscChecker utility.  This 
program needs to be run from an experiment results directory, and must be run after the Torrrent Variant
Caller (TVC) plugin has been run.  

    -l    Non-default lookup file for the cpscChecker utility (default: mc)
    -s    Custom sampleKey file rather than the default based on the ion_params_00.json file
    -r    Run plugin without cpscChecker (R&D Server without standard sample IDs)
    -v    Version information
    -h    This help text

EOT
)"

resultsDir="$(pwd)"
colvarsDir="$resultsDir/collectedVariants"
#runName="$(echo "$resultsDir" | perl -pe 's/.*user_((?:[P|M]CC-\d+|MC[12]))-\d+.*/\1/')"
runName="$(echo "$resultsDir" | egrep -o '[PM]C[C123]-[0-9]+')"

cpscSample="IonXpress_001.txt"
cpsc_lookup=mc #Default lookup file for cpscChecker
is_RandD_server=0 #Change to 0 for production server with locked pipeline.

# Get absolute path of scriptname in order to be more flexible later. 
SCRIPT=$(readlink -f $0)
SCRIPTPATH=$(dirname $SCRIPT)

while getopts :hvrs:l: opt;
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
        r)
            is_RandD_server=1
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

# since TSSv4.4 is now allowing to have multiple plugins run and the results residing for them in a numbered dir, 
# grab the dir names. store in a hash, and figure out th elarges key which will be the latest TVC results.
echo "Checking for TVC output dir..."
declare -A tvc_results
for dir in ${resultsDir}/plugin_out/*; do 
    dir=$(basename $dir)
    plugin_results=(${dir//./ })
    if [[ ${plugin_results[0]} != variantCaller_out ]]; then 
        #echo "skipping over '${plugin_results[0]}'..."
        continue
    elif [[ "${plugin_results[1]+isset}" ]]; then
        echo "${plugin_results[1]}  => ${plugin_results[0]}"
        tvc_results[${plugin_results[1]}]=$dir
    else
        tvc_results[0]=$dir
    fi
done

if [[ ${#tvc_results[@]} -eq 0 ]]; then
    echo "ERROR: TVC has not been run on this sample. Please run TVC before running this script!"
    exit 1
else
    largest=0
    for key in ${!tvc_results[@]}; do 
        if [[ $key > $largest ]]; then
            largest=$key
        fi
    done
    tvc_output=${tvc_results[$largest]}
    echo "TVC output directory '$tvc_output' found. Continuing."
    echo
fi

tvc_output="plugin_out/$tvc_output"

# Check to make sure TVC has been run and you are in the results dir
if [[ ! -d $tvc_output ]]; then
	echo "[ ERROR ] Either you are not running this script from a run results folder, or TVC has not been run on this sample."
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
echo "Making a collectedVariants dir..."
mkdir $colvarsDir;
# Use custom key instead of generating one and copy to $colvarsDir
if [[ $customKey ]]; then
	cp $customKey "$colvarsDir/sampleKey.txt"
else
    echo "Generating sampleKey.txt..."
    eval "perl $SCRIPTPATH/sampleKeyGen.pl -o $colvarsDir/sampleKey.txt" 
	if [[ $? -ne 0 ]]; then
		echo "The sampleKeyGen tool had a problem and failed to run correctly!" 
		exit 1
	fi
fi

# Run collectVariants.pl
echo "Running collectVariants.pl on TVC results..." 
eval "perl $SCRIPTPATH/collectVariants.pl" 
if [[ $? -ne 0 ]]; then
	echo "[ ERROR ] The collectVariants.pl script had problems and failed to run!" 
	exit 1
fi

if [ ! -d "$colvarsDir" ]; then 
	echo "[ ERROR ] No 'collectedVariants' directory found. Something bad happened!"
	exit 1
else
    cd $colvarsDir
fi

if [[ $is_RandD_server -eq 1 ]]; then
	echo "[ NOTE ]: Running in an R&D environment.  No automatic cpscChecker run.  Please run the cpscChecker tool manually."
else
	# Verify conditions are correct for cpscChecker 
	if [[ ! -f "$colvarsDir/sampleKey.txt" ]]; then
		printf "[ ERROR ] No sampleKey file found.  Can not run cpscChecker!\n"
    # Going to stick with the tab delimited file here since the VCF file name is not reliable (I rename these based
    # on the sample key and if there is a non-standard sample naming, we'll lose the sample unintentionally.  
    #elif [ ! -f "$cpscSample" ]; then
    elif [ ! -f "$cpscSample" ]; then
		printf "[ ERROR ] No "$cpscSample" file found.  Did you run a CPSC sample?\n";
		exit 1
	fi	
		
	#  If conditions are OK, let's run the script.
	printf "Running cpscChecker for '$runName' on sample '$cpscSample'...\n"
	#eval "perl $SCRIPTPATH/../cpscChecker/cpscChecker.pl -l $cpsc_lookup $cpscSample"
	eval "perl $SCRIPTPATH/../cpscChecker/cpscChecker.pl -l $cpsc_lookup -T $cpscSample -o cpscChecker_output.txt"
fi

# Prepare a zip archive of results for experiment folder 
zip "$runName"_variants.zip * > /dev/null 2>&1

echo -e "\nScript complete. Results written to the 'collectedVariants' directory\n" 
