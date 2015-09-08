#!/usr/bin/python
# Script to coordinate the collection of variants from a TVC analysis, naming and concatenation 
# of the .txt files containing the variants, and running of cpscChecker.pl to analyze the
# CPSC sample performance.
#
# This is a overhaul of the original shell script with Python since there are some issues with
# some of the array implementations.
###############################################################################################
import sys
import re
import argparse
from pprint import pprint

version = '0.0.1_000000'

def get_args():
    parser = argparse.ArgumentParser(
        formatter_class = lambda prog: argparse.HelpFormatter(prog, max_help_position = 100, width=200),
        description='''
        Program to collect variant call files from an Ion Torrent run and run the cpscChecker utility. This 
        program needs to be run from an experiment results directory, and must be run after the Torrent
        Variant Caller (TVC) plugin has been run.
        ''',
        version = '%(prog)s  - ' + version,
        )
    parser.add_argument('-l', '--lookup', default='mc', metavar='<lookup>',
            help='Non-default lookup file for the cpscChecker utility (Default: %(default)s)')
    parser.add_argument('-s', '--samplekey', metavar='<samplekey>',
            help='Custom sample key to use for this run')
    parser.add_argument('-r', '--r_and_d', help='Run plugin without cpscChecker (R&D server without standard sample IDs')
    args = parser.parse_args()
    return args



def main():
    cli_opts = get_args()


if __name__ == '__main__':
    main()
