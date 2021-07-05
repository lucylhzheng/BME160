# BME160
Programming in Life Science


Lab 2
seqCleaner.py reads a DNA string from user input and return a collapsed substring of embedded Ns to: {count}.
Replace block of (N)s with the {count} determined as the output.

fastqParse.py reads the inputted Seqname of FASTQ file and each section of the file is parsed to generate
individual fields of outputs.

coordinateMathsoln.py takes in Three sets of atomic coordinates to generate respective bond lengths and bond angles using programmed mathemtical functions.

converter.py accepts nucleotide sequences and gets converts them to their respective dictionaries of RNA,DNA, and amino acids which outputs the 
the correct conversion based on the type of given data.

Lab 3
proteinParams.py takes in a string of amino acids and calculates the physical-chemical properties related to the given protein sequence.

Lab 4
SequenceAnalysis is a module that contains three classes which are NucParams, FastAreader, and ProteinParams
which can be called from the program titled GenomeAnalyzer. All the calculations of amino acids, nucleotides,
and codons from the data inputted will be done here only if it is called from GenomeAnalyzer and will produce
the output.

Lab 5
SequenceAnalysis-1 is a module that contains four classes which are OrfFinder, NucParams, FastAreader, and ProteinParams
which can be called from the program titled findORFs. The generation of frame shifting, start position, end position, and
gene length is generated within the OrfFinder.

Lab 6
The findUnique.py file takes in fasta file of tRNA sequences from STDIN and generates an output file of the essential sets for each tRNA.
