# pscan

## Installation

Requires Leiningen to build/run:

http://leiningen.org/#install

(Make sure lein, java and javac are on your $PATH)

    lein repl
    
in the project directory to start the read-eval-print-loop. (This may take some time and produce a long stream of diagnostic messages during the first execution, as Leiningen has to pull in dependencies and compile Java code.)

## Usage

    => (require 'pscan.core)
    
At the REPL to load the project. The var pscan.core/sample-proteins is pre-populated with 50 sample proteins from 5 different UniProt clusters

    => (def a-protein (get pscan.core/sample-proteins 50))
    
Will bind the 50th protien to the var a-protein

    => (:sequence a-protein)
    => (:description a-protein)
    => (:file a-protein)

Sequence, description and source file of a-protein, respectively.

    => (time (def clusters (dbscan pscan.core/sample-proteins pscan.core/protein-metric 250 4)))

Cluster pscan.core/sample-proteins using the DBSCAN algorithm, BLOSUM62 as a similarity score, requiring a score of 250 or greater to cluster, with at least 4 proteins per cluster; bind the result to the var clusters, and print the execution time.

    => (def pscan.core/region-query pscan.core/region-query-s)
    
Region queries (DBSCAN subroutine) in serial instead of parallel

    => (def pscan.core/memoized-dbscan false)
    
DBSCAN without memoization, trading space for time.

    => (doc pscan.core/some-function)
    
See documentation for some-function. (Or see doc/index.html)

## Acknowledgements

This project incorporates code (src/java/match/Match2.java) from Peter Sestoft's implementation of the Needleman-Wunsch algorithm. (http://www.itu.dk/people/sestoft/bsa.html)

The project uses the BioClojure library to read FASTA sequences (http://www.ncbi.nlm.nih.gov/pubmed/24794932)

## License

Copyright © 2014 Joshua Brule

Distributed under the Eclipse Public License either version 1.0 or (at
your option) any later version.
