MobiDB-lite, long disorder consensus predictor
==============================================
Marco Necci, Damiano Piovesan and Silvio C.E. Tosatto  

Version 3.8.3

Introduction
------------
MobiDB-lite executes up to 9 different disorder predictors, collects the outputs
and calculates a consensus. The consensus is generated by measuring predictors
agreement. At least 62.5% of predictors must agree to assign disorder state to a
residue. Then a mathematical morphology (MM) dilation/erosion processing is
applied. Finally short regions are filtered out.

Requirements
------------
MobiDB-lite works both for 32-bit and 64-bit architectures in Linux system. It
requires Java (version >=7.x) and Python (version 3.x.x). By default it uses
64-bit executables, but 32-bit version can be called by a command line option
(see below).

Usage
-----
The following command prints the help page on screen::

    python3 mobidb_lite.py -h


By default, MobiDB-lite searches the ``binx/`` folder (that includes predictors
executables) in its own directory. So, MobiDB-lite can be executed by just
providing a multi-fasta input::

    python3 mobidb_lite.py /examples/multi.fasta

MobiDB-lite is able to use predictors executables from other paths. See **Configuration**
section for further details.

By default, MobiDB-lite runs on single thread but it can exploits up to 7 threads.
When less than 7 CPUs are available, the "-t" option allows to calibrate the load::

    python3 mobidb_lite.py examples/multi.fasta -t 4

Configuration
-------------

MobiBD-lite is designed to work out of the box. For a basic usage you should not
need to configure anything. Nonetheless, MobiDB-lite offers some configuration
possibilities.

### Bin directories

MobiDB-lite is able to use predictors executables from other paths. Paths of
single predictors can be configured in the ``config.ini`` or custom configuration file.
You should always use **absolute paths** when configuring binx directories.
Notice that when paths ar not set in this file, MobiDB-lite will look for the
executables in the ``binx/`` folder.

### Predictors thresholds
Most predictors' output is a score ranging from 0 to 1. To transform into a binary
prediction a threshold is applied. Thresholds are usually defined in publications
relative to single predictors. Predefined thresholds are obtained from such publications.
However, thresholds of single predictors can be configured in the ``config.ini`` or
custom configuration file.

### Logging
Log file and log level can be configured through command line arguments, ``--log`` and
``--logLevel`` respectively. Default log level is set to ERROR.

Output
------

The output is printed on screen by default. To save it on a file the "-o"
option has to be used::

    python3 mobidb_lite.py examples/multi.fasta -o out_file.txt


The output is provided in 3 different formats. The output format can be chosen
with the ``-f`` option, which accepts one of the following values: `interpro`. 
`fasta`, `vertical`, `extended`, `mobidb4`, `caid`. The default short version 
(option `interpro`) includes the protein name and start/end position of the
disorder regions, one region per line. All elements are separated by a TAB::

    sp|Q15648|MED1_HUMAN	609	706
    sp|Q15648|MED1_HUMAN	652	681	Polar
    sp|Q15648|MED1_HUMAN	792	820
    sp|Q15648|MED1_HUMAN	806	820	Polar
    sp|Q15648|MED1_HUMAN	874	893
    sp|Q15648|MED1_HUMAN	947	1566
    sp|Q15648|MED1_HUMAN	995	1020	Polyampholyte
    sp|Q15648|MED1_HUMAN	1078	1156	Low complexity
    sp|Q15648|MED1_HUMAN	1158	1197	Polar
    sp|Q15648|MED1_HUMAN	1219	1292	Low complexity
    sp|Q15648|MED1_HUMAN	1328	1348	Polar
    sp|Q15648|MED1_HUMAN	1349	1366	Polyampholyte
    sp|Q15648|MED1_HUMAN	1367	1382	Low complexity
    sp|Q15648|MED1_HUMAN	1421	1484	Polar
    sp|Q15648|MED1_HUMAN	1508	1530	Polyampholyte
    sp|Q15648|MED1_HUMAN	1531	1556	Polar


Sequence features can be silenced by using option -sf, inactive by
default. When using this option, sequence features are skipped::

    sp|Q15648|MED1_HUMAN	609	706
    sp|Q15648|MED1_HUMAN	792	820
    sp|Q15648|MED1_HUMAN	874	893
    sp|Q15648|MED1_HUMAN	947	1566


Option `fasta` gives the output in fasta format (new in version 3.8.2).
In this format `D` is disorder, `S` is structure and region subtypes are 
given as figures from 1 to 8::
    
    sp|Q15648|MED1_HUMAN
    SSS...
    SSSSSSSSDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
    D888888888888888888888888888888DDDDDDDDDDDDDDDDDDD
    DDDDDDSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
    SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSDDDDDDDDD
    DDDDD888888888888888SSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
    SSSSSSSSSSSSSSSSSSSSSSSDDDDDDDDDDDDDDDDDDDDSSSSSSS
    SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSDDDD
    DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD111111
    11111111111111111111DDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
    DDDDDDDDDDDDDDDDDDDDDDDDDDD77777777777777777777777
    77777777777777777777777777777777777777777777777777
    777777D8888888888888888888888888888888888888888DDD
    DDD...

Sequence features can be silenced by using option -sf, inactive by
default. When using this option, sequence features are expressed as `D`::

    sp|QDDDDD|MEDD_HUMAN
    SSS...
    SSSSSSSSDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
    DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
    DDDDDDSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
    SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSDDDDDDDDD
    DDDDDDDDDDDDDDDDDDDDSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
    SSSSSSSSSSSSSSSSSSSSSSSDDDDDDDDDDDDDDDDDDDDSSSSSSS
    SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSDDDD
    DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
    DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
    DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
    DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
    DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
    DDD...


Option `vertical` gives the output with one residue per line (new in version 3.8.2).
In this format each line will will have three values: `res  score  state`. Where `res`
is the amino acid one letter code, `score` is the agreement score of the consensus, 
`state` if the predictor decision between order (`S`) and disorder (`D`) ::

Option `extended` provides InterPro extended output (only MobiDB-lite prediction, skip
fully structured proteins). Contains the binary MobiDB-lite prediction and
the IDRs boundaries. Also contains a relaxed consensus (rlx_consensus) calculated
at a threshold of .375 (New in version 3.8.1) ::

    {
      "accession": "sp|Q15648|MED1_HUMAN",
      "consensus": "SSS...SSS",
      "rlx_consensus": "DDD...DDD",
      "regions": [
        [609, 706, "D"],
        [792, 820, "D"],
        [874, 893, "D"],
        [947, 1566, "D"]
      ]
    }


Option `mobidb4` provides an extended output and structured output that is used
to create annotations for MobiDB4 entries::
    {
    
        "acc": "sp|Q15648|MED1_HUMAN",
        "sequence": "MKAQASQAL...",
        "length": 1581,
        "prediction-disorder-mobidb_lite": {
            "regions": [[609, 706], [792, 820], [874, 893], [947, 1566]],
            "scores": [0.75, 0.75, 0.75, 0.875, 0.75, 0.875, ... ],
            "content_count": 767, 
            "content_fraction": 0.485
        },
        "prediction-disorder-th_50": {...},
        "prediction-low_complexity-merge": {...},
        "prediction-disorder-iupl": {...},
        "prediction-disorder-iups": {...},
        "prediction-disorder-espN": {...},
        "prediction-disorder-espD": {...},
        "prediction-disorder-espX": {...},
        "prediction-disorder-glo": {...},
        "prediction-disorder-dis465": {...},
        "prediction-disorder-disHL": {...},
        "prediction-disorder-vsl": {...},
        "prediction-low_complexity-seg": {...},
        "prediction-low_complexity-pfilt": {...},
        "prediction-helix-fess": {...},
        "prediction-sheet-fess": {...},
        "prediction-coil-fess": {...},
        "prediction-rigidity-dynamine": {...},
        "prediction-lip-anchor": {...}
    }


Option `caid` provides a very extended output, containing both scores and regions
for each predictor, but it's restricted to disorder predictors plus DynaMine::

    {
      "accession": "sp|P49137|MAPK2_HUMAN",
      "predictions": [
        {
          "method": "mobidb_lite",
          "scores": [0.875, 1.0, ...
          ],
          "regions": [
            [1, 43, "D"]
          ]
        },
        {...}
      ]
    }

Troubleshooting
---------------
If you are not seeing any output try to use the -fc option to force output
in case of single predictors failures.

MobiDB-lite won't work on OS-X since some of the underlying predictors are
not compiled for it.

If you encounter any bug or strange behaviour, please don't hesitate and
contact us at biocomp@bio.unipd.it.
