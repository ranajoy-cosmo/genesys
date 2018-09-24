import sys
from termcolor import colored

class Generic_Class(object):
    pass

def prompt(text, outstream=sys.stdout, color=None):
    if color:
        text = colored(text, color) 
    outstream.write(text)
    outstream.flush()
