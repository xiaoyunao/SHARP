"""Public field-design entry: compare baseline and planned-visit-selected alternatives."""
from .compare_campaigns import compare,main

def design(config,start,cache,output):
    return compare(config,start,cache,output)

if __name__=='__main__':main()
