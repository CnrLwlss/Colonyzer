import colonyzer2 as c2
import argparse

def parseArgs(inp=''):
    '''Define console script behaviour, hints and documentation for merging Colonyzer files and adding metadata.'''
    parser=argparse.ArgumentParser(description="Tool for concatanating Colonyer output (.out files) and optionally adding metadata from an experimental description file, a library description file and a file mapping systematic to standard gene names.")

    parser.add_argument("-i","--imdir", type=str, help="Directory containing Colonyzer output (.out files).", default=".")
    parser.add_argument("-m","--mdir", type=str, help="Directory containing metadata files.", default=".")
    parser.add_argument("-e","--edesc", type=str, help="Filename for optional tab-delimited text file describing experiment.  File found in directory specified by -m.", default="ExptDescription.txt")
    parser.add_argument("-l","--ldesc", type=str, help="Filename for optional tab-delimited text file describing library examined during experiment.  File found in directory specified by -m.", default="LibraryDescription.txt")
    parser.add_argument("-g","--g2orf", type=str, help="Filename for optional tab-delimited text file translating systmatic (ORF) to standard gene names.  File found in directory specified by -m.", default="ORF2GENE.txt")
    parser.add_argument("-o","--out", type=str, help="Output filename.", default="ColonyzerOutput.txt")
    parser.add_argument("-f","--fmt", type=str, help="Timestamp string format.", default="%Y-%m-%d_%H-%M-%S")

    if inp=="":
        args = parser.parse_args()
    else:
        args = parser.parse_args(inp.split())
    return(args)

def main():
    inp=parseArgs()
    res=c2.parseAndCombine(inp.imdir,inp.mdir,inp.edesc,inp.ldesc,inp.g2orf,inp.out,inp.fmt)

if __name__ == "__main__":
    main()


