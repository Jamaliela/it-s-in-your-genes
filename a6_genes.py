######################################################################
# Author: Ela Jamali, Emely Alfaro
# Username: jamalie, alfarozavalae
#
# Assignment: A06: It's in your Genes
#
# Purpose: Determine an amino acid sequence given an input DNA sequence
#Credit: Giorgi Lomia - TA who helped us in this assignment
######################################################################
# Acknowledgements:
#   Original Author: Dr. Jan Pearce
#
#   Idea from: http://www.cs.uni.edu/~schafer/1140/assignments/pa9/index.htm
#
# licensed under a Creative Commons
# Attribution-Noncommercial-Share Alike 3.0 United States License.
####################################################################################


def is_nucleotide(sequence):
    """
    This function checks if the user's input contains the accepted letters or not. if not false is returned
    Returns True if so, and False otherwise

    :param sequence: the user's input that will be evaluated to see if it is part of the list s
    :return: Boolean Values True or False
    """
    s = ["C","A","G","T"]                           # the list of 4 nucleotides
    for i in range(len(sequence)):                  # loop to check if the elements entered are part of the user's input (sequence)
        if sequence[i] not in s:                    # if statement to check if each of the elements are not part of the list
            return False                            # as the loop is running we cant exit it when a condition is evaluated true, therefore we do false first
    return True                                     # returning true when checked that all the string is part of the list


def complement_strand(sequence):
    """
    This function changes the string into their complements and returns the complement.
    If there is a bad input, the function returns "Sequencing Error"

    :param sequence: uses the string entered by the user as input
    :return: complement
    """
    if not is_nucleotide(sequence):
        return "Sequencing Error"

    complement = ""                                 # This can be used to "build" the complement
    for i in range(len(sequence)):                  # loop to go over the length of the sequence input
        if sequence[i] == "T":                      # if function to evaluate if the each of the letters of sequence is T then:
            complement = complement + "A"           # the complement will start building up by adding the new letter A instead of T
        if sequence[i] == "C":                      # if function to evaluate if the each of the letters of sequence is C then:
            complement = complement + "G"           # the complement will continue building up by adding the new letter G instead of C
        if sequence[i] == "A":                      # if function to evaluate if the each of the letters of sequence is A then:
            complement = complement + "T"           # the complement will continue building up by adding the new letter T instead of A
        if sequence[i] == "G":                      # if function to evaluate if the each of the letters of sequence is G then:
            complement = complement + "C"           # the complement will continue building up by adding the new letter C instead of G
    return complement


def mRNA(sequence):
    """
    This function replaces each occurrence of the nucleotide T replaced with the nucleotide U.

    :param sequence: the string that the user input as the dna
    :return: mrna
    """
    mrna = ""
    for i in range(len(sequence)):                  # loop to go over the length of the sequence input
        if sequence[i] == "T":                      # if function to evaluate if each of the letters of sequence is T then:
            mrna = mrna + "U"                       # if a letter T is found then the mrna starts building up to add a U instead of T
        else:                                       # when there is not a letter U found  then
            mrna = mrna + sequence[i]               # the mran will build up by aadding the element evaluated[i]
    return mrna                                     # the final mrna is returned


def chunk_amino_acid(sequence):
    """
    This function takes the mrna and separates it into chunks of three so that those can be returned and the remainder can be discarded.

    :param sequence: this function uses floor division to cut in chunks the string and return only the chunks of three
    :return: list_of_chunks
    """
    list_of_chunks = []                            # starting the variable so the chunks can be build
    n = 0                                          # starting variable at zero so it can be incrementing
    for i in range(len(sequence)//3):              # for loop to go over the different chunks of three and ignore the remainder
        list_of_chunks.append(sequence[n:n+3])     # list to append each of chunks of dna
        n = n+3
    return list_of_chunks


def amino_acid_chunks(threecharseq):
    """
    Expects a three character string as a parameter and returns
    the corresponding single character AminoAcid

    :param threecharseq: a sequence of three characters
    :return: A string representing the amino acid chunk for that sequence
    """

    # ###################################################################### #
    # #  This function is already completed correctly! No changes needed!  # #
    # ###################################################################### #

    # We haven't learned about dictionaries yet, but here is one for the extra curious.
    # You aren't expected to know what this is yet.
    translator = {  "GCA":"A", "GCC":"A", "GCG":"A", "GCU":"A",
                    "AGA":"R", "AGG":"R", "CGA":"R", "CGC":"R", "CGG":"R", "CGU":"R",
                    "GAC":"D", "GAU":"D",
                    "AAC":"N", "AAU":"N",
                    "UGC":"C", "UGU":"C",
                    "GAA":"E", "GAG":"E",
                    "CAA":"Q", "CAG":"Q",
                    "GGA":"G", "GGC":"G", "GGU":"G", "GGG":"G",
                    "CAC":"H", "CAU":"H",
                    "AUA":"I", "AUC":"I", "AUU":"I",
                    "UUA":"L", "UUG":"L", "CUA":"L", "CUC":"L", "CUG":"L", "CUU":"L",
                    "AAA":"K", "AAG":"K",
                    "AUG":"M",
                    "UUC":"F", "UUU":"F",
                    "CCA":"P", "CCC":"P", "CCG":"P", "CCU":"P",
                    "AGC":"S", "AGU":"S", "UCA":"S", "UCC":"S", "UCG":"S", "UCU":"S",
                    "ACA":"T", "ACC":"T", "ACG":"T", "ACU":"T",
                    "UGG":"W",
                    "UAC":"Y", "UAU":"Y",
                    "GUA":"V", "GUC":"V", "GUG":"V", "GUU":"V",
                    "UAA":"*", "UAG":"*", "UGA":"*"
                 }

    if threecharseq in translator.keys():
        return translator[threecharseq]     # Given any 3 letter sequence, this returns the amino acid for that sequence
    else:
        return "?"                          # Returns a question mark if the input is invalid


def sequence_gene(sequence):
    """
    The sequence_gene() function takes a a sequence of nucleotides:
    A, C, G, and T and returns
    the corresponding amino acid sequence.

    :param sequence: a string representing a sequence of nucleotides
    :return: a string representing the amino acid sequence
    """

    # ###################################################################### #
    # #  This function is already completed correctly! No changes needed!  # #
    # ###################################################################### #

    aaseq=""                                                # Amino acid sequence
    if is_nucleotide(sequence):                             # Checks for a valid sequence
        comp_strand = complement_strand(sequence)           # Finds the complement sequence
        mrna = mRNA(comp_strand)                            # Finds the mRNA of the complement
        amino_acid_list = chunk_amino_acid(mrna)            # Chunks the mRNA sequence

        for amino_acid in amino_acid_list:                  # Loops through each chunk
            aaseq = aaseq + amino_acid_chunks(amino_acid)   # Creates the final amino acid sequence
    return aaseq                                            # Returns an empty string for any illegal input


def main():
    """
    The main() function runs the sequence_gene code, which calls all other parts of this code

    :return: None
    """
    sequence = input("Please enter your DNA strand. (DNA strands are formed of C,A,G,T)")
    sequence_gene(sequence)


    # TODO When your code is fixed, the following line will print correctly.
    # TODO You do not need to modify the sequence_gene() function; it is correct already.
    print("The original sequence {0} returns {1}".format("CACGT", sequence_gene("CACGT")))


if __name__ == "__main__":
    main()          # call to main() function
