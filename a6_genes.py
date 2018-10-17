######################################################################
# Author: Ela Jamali, Emely Alfaro
# Username: jamalie, alfarozavalae
#
# Assignment: A06: It's in your Genes
#
# Purpose: Determine an amino acid sequence given an input DNA sequence
#
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
    Checks that the string sequence provided is a valid string
    consisting only of the 4 nucleotides A, C, G, and T
    Returns True if so, and False otherwise

    :param sequence: the user's input that will be evaluated to see if it is part of the list s
    :return: Boolean Values
    """
    s = ["C","A","G","T"]                           # the list of 4 nucleotides
    for i in range(len(sequence)):                  # loop to check if the elements entered are part of the user's input (sequence)
        if sequence[i] not in s:                    # if statement to check if each of the elements are not part of the list
            return False                            # as the loop is running we cant exit it when a condition is evaluated true, therefore we do false first
    return True                                     # returning true when checked that all the string is part of the list

    # FIXME: Finish the docstring above

    # FIXME: Add your code here to complete this function

    # FIXME: Clearly this should not always return True


def complement_strand(sequence):
    """
    Returns the string which will be the second strand of the DNA sequence
    given that Ts complement As, and Cs complement Gs. If given
    a bad input, the function returns "Sequencing Error"

    :param sequence:
    :return:
    """
    # FIXME: Finish the docstring above
    if not is_nucleotide(sequence):
        return "Sequencing Error"

    complement = ""        # This can be used to "build" the complement
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
    # FIXME: Add your code here to compete this function

    # FIXME: Obviously, this should not always return ""

def mRNA(sequence):
    """
    Replaces each occurrence of the nucleotide T replaced with the nucleotide U.

    :param sequence:
    :return:
    """
    # FIXME: Finish the docstring above

    # FIXME: Add your code here to complete this function

    mrna = ""
    for i in range(len(sequence)):                  # loop to go over the length of the sequence input
        if sequence[i] == "T":                      # if function to evaluate if each of the letters of sequence is T then:
            mrna = mrna + "U"                       # if a letter T is found then the mrna starts building up to add a U instead of T
        else:                                       # when there is not a letter U found  then
            mrna = mrna + sequence[i]               # the mran will build up by aadding the element evaluated[i]
    return mrna                                     # the final mrna is returned
    # FIXME: Obviously, this should not always return ""



def chunk_amino_acid(sequence):
    """
    Uses output of mRNA(sequence) and divides it into substrings of length 3,
    ignoring any "extra DNA" at the far end returning the relevant substrings in a list.

    :param sequence:
    :return:
    """
    # FIXME: Finish the docstring above

    # FIXME: Add your code here to complete this function
    list_of_chunks = []
    n = 0
    for i in range(len(sequence)//3):
        list_of_chunks.append(sequence[n:n+3])
        n=n+3

    # FIXME: Obviously, this should not always return an empty string
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
    nucleotide = input("Please enter your DNA strand. (DNA strands are formed of C,A,G,T")
    is_nucleotide(nucleotide)


    # TODO When your code is fixed, the following line will print correctly.
    # TODO You do not need to modify the sequence_gene() function; it is correct already.
    print("The original sequence {0} returns {1}".format("CACGT", sequence_gene("CACGT")))


if __name__ == "__main__":
    main()          # call to main() function
