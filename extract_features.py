import sys
import argparse

def run():
    

def get_parser():
    parser = argparse.ArgumentParser(description =
    'Extract feature vectors from PSMs and MS2 spectrua')
    parser.add_argument('-r', '--result-file', help='supported only msgf+ format')
    parser.add_argument('-s', '--spectrum-file', help='supported only mgf foramt')
    
    return parser

def command_line_runner():
    parser = get_parser()
    args = vars(parser.parse_args()) 
    
    print args
    print "done"

    run()

if __name__ == '__main__':
    command_line_runner()
    
