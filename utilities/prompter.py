import sys
from termcolor import colored

def prompt(text, outstream=sys.stdout, color=None):
    """
    Prompt any text to the outstream provided immediately. The text is flushed out.
    """
    if color:
        text = colored(text, color) 
    outstream.write(text)
    outstream.flush()
