# Python_Homework_1 description.

**Python_Homework_1 contains the tool for distance computation between each pair of nucleotide sequences in input FASTA file.** Scoring system attributes values for the letters that are involved and for the type of alteration. The next values are used as penalties:
* <ins>For substitutions:</ins>
  * **-1** for A>G and G>A;
  * **-3** for A>C, C>A, G>T and T>G;
  * **-4** for A>T and T>A;
  * **-5** for G>C and C>G;
  * **0** for C>T and T>C;
* <ins>**-5** for gaps;</ins>
* <ins>For matches:</ins>
  * **10** for A-A;
  * **7** for G-G;
  * **9** for C-C;
  * **8** for T-T.

*As seen above, the most similar sequences will have the maximum score.*
## Help

Arguments:
* -i (--input): the path to input FASTA file, **necessary to specify**
* -t (--threads): number of threads, default value = 1
* -o (--output): the resulting filename, defalt name *'/home/izmaylova/Documents/Matrix.odt'*

The resulting file contains the matrix with pairwise scores for all sequences from input FASTA file.
