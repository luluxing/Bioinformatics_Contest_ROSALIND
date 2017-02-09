# Bioinformatics_Contest_ROSALIND

## Introns Detection

The DNA strand consists of exons and introns. The exons are the coding parts of the DNA, while the introns are not. During the process of transcription, DNA sequence of a gene gets copied into an RNA molecule, where the introns are removed during a process called RNA splicing.

Recently, BioLabs Inc. developed an instrument to get some information on what gets spliced out of a gene sequence. Using black magic a complicated biochemical process this instrument produces a number of "reads": substrings of a DNA sequence after the splicing. They hired you to write a program to restore a possible DNA sequence after the splicing operation: that is a sequence that contains as many reads as possible as substrings and that can be obtained from original DNA sequence by a splicing procedure.

In other words, there is a string S, the splicing results in a string T which is a subsequence of S. Your task it to restore T given S and a number of T's substrings.

## Notes

In tests almost all the DNA sequences, exons and introns are generated using real-life data. The reads are emulated using known result of DNA splicing. That means that for each test there exists an answer that satisfies all the reads.

Note, that the tests are approximately sorted by their sizes and not specifically by difficulty. Probably, it will be easier to solve test 6 rather than 4, for example.

## Input Format

The first line of the input contains an original DNA sequence represented as a string of characters A, C, G, T. The second line contains one integer nn --- the number of reads. Each of the next nn lines contains one read each.

## Output Format

The first line of the output should contain a possible DNA sequence after the splicing. Next nn lines should contain one integer each - the position of the corresponding read in the proposed DNA sequence after the splicing. The position are counted starting from one. If your output doesn't contain some of the reads then you can put -1 in the corresponding line.
