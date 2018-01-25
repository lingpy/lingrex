# LingRex: Linguistic Reconstruction with LingPy

The basic idea of this repository is to provide initial functionalities for a linguistic reconstruction engine in LingPy. The workflow that we follow to arrive at this ambitious goal involves alignments, and pre-defined cognate sets, and employs correspondence pattern analysis.

## Tentative Workflow

1. identify cognates in the data and align all cognate sets
2. carry out a correspondence pattern analysis in which most frequent correspondence patterns are identified (a network analysis)
3. carry out a simple reconstruction on each of the correspondence patterns, either using parsimony and individual step-matrices for the sound transitions, or majority rules (as a baseline)
4. identify the best candidates for a full reconstruction by determining cognate sets which have enough reflexes in the data (using a user-defined threshold), and add empty slots for each of their position in the alignments
5. compare the positions with the patterns plus their proto-sounds inferred in steps 2 and 3 and add all possible values
6. use decision trees on multi-tiered sequence models to find out how predictive the system is


## Important Note on the Workflow

* make a post-parse of correspondence patterns, after having them clustered the first time. This first step serves only to determine generally possible patterns. The evaluation of each pattern needs to be checked against the data again and all full patterns should be checked on this. The current scores are not very reliable, as they suffer from heavy gaps. What we'll need are post-parses in which we check the real theoretical potential of each pattern, and also get a score on the coverage of each position (e.g., sometimes, we have a long cluster, but only one sound in one slot for a given language. This is of course not what we want).
