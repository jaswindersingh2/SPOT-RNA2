-------------------------------------------------------------------------

FreeKnot

-------------------------------------------------------------------------

Authors: Jimmy Ka Ho Chiu and Yi-Ping Phoebe Chen

Last updated on 15 Apr 2014



-------------------------------------------------------------------------

Purpose



FreeKnot is a tool for RNA pseudoknot removal. It converts any pseudoknot

into nested substructures in RNA secondary structures. It removes some

crossing stems to eliminate crossings based on certain scoring functions

(details will be provided later in this README file) and reports one or

more optimized pseudoknot-free structures.



-------------------------------------------------------------------------

Platform and pre-requisites



FreeKnot has been tested on various platforms including Linux (Ubuntu),

Mac OS X and Windows. Perl (v5.14 or later) is recommended. Earlier

versions might work but without guarantee. Windows users can download

various Perl distributions for Windows. ViennaRNA package 2.1 is required

for the free energy scoring function.



-------------------------------------------------------------------------

Program/Module Description



BpseqParser.pm, DPParser.pm     - parser to accept bpseq or

                                  dot-parentheses formats as input

BpseqWriter.pm, DPWriter.pm     - writer to output converted results in

                                  bpseq or dot-parentheses formats

ChordModel.pm, CircleGraph.pm   - graphical object for primitive

                                  pseudoknot representation

MIS.pm				- MIS algorithm (for free energy scoring

				  function) 
MWIS.pm                         - MWIS algorithm

ScoringFunctions.pm             - scoring functions

remove_pseudoknot.pl	        - main program for pseudoknot removal

PrimitivePseudoknotExtractor.pm - primitive pseudoknot extraction from

                                  the input secondary structure

BracketPairs.pm                 - processing brackets in input secondary

                                  structure


VertexSubset.pm			- subset objects for storing graph
				  vertices in the MIS algorithm

-------------------------------------------------------------------------

Usage



FreeKnot is executed in console. The command is:



perl remove_pseudoknot.pl -i <secondary structure format of input file>

    -s <scoring function option> <input file path>



Secondary structure format available: dp (dot-parentheses) / bpseq

The secondary structure format for the output file follows that of the

input file. So, if the input file is in bpseq format then the output

file is also in bpseq format. Note that every line of data must end with

a newline character (i.e. \n).



Scoring function options: bp (# of base pairs) / stem (# of base pair

stems) / hb (# of hydrogen bonds) / fe (structure overall free energy)



The results are outputted to the console (stdout) by default. They can be

directed to a file. For example,



perl remove_pseudoknot.pl -i bpseq -s bp input.bpseq > output.bpseq



-------------------------------------------------------------------------
