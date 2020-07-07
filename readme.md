MobiDB-lite, long disorder consensus predictor
==============================================
Marco Necci, Damiano Piovesan and Silvio C.E. Tosatto
Version 3.8.2

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
`fasta`, `vertical`, `extended`, `mobidb3`, `caid`. The default short version 
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


Option `mobidb3` provides an extended output and structured output that is used
to create annotations for MobiDB3 entries::

    {
      "sequence": "MKA...IGN",
      "mobidb_consensus": {
        "disorder": {
          "predictors": [
            {
              "method": "mobidb_lite",
              "regions": [
                [609, 706, "D_WC"],
                [792, 820, "D_WC"],
                [874, 893, "D_WC"],
                [947, 1566, "D_WC"]
              ],
              "scores": [
                0.75,
                0.75,
                0.75,
                ...
                0.625,
                0.75,
                0.75
              ],
              "dc": 0.48513598987982287
            },
            {
              "method": "mobidb_lite_sub",
              "regions": [
                [652, 681, "D_PO"],
                [806, 820, "D_PO"],
                ...
              ]
            },
            {
              "method": "simple",
              "regions": [
                [1, 16, "D"],
                [89, 89, "D"],
                ...
              ],
              "dc": 0.6299810246679317
            }
          ]
        }
      },
      "mobidb_data": {
        "disorder": {
          "predictors": [
            {
              "method": "iupl",
              "regions": [
                [1, 15, "D"],
                [82, 82, "D"],
                ...
              ]
            },
            {
              "method": "iups",
              "regions": [
                [1, 13, "D"],
                [538, 572, "D"],
                ...
              ]
            },
            {
              "method": "espN",
              "regions": [
                [1, 33, "D"],
                [49, 57, "D"],
                ...
              ]
            },
            {
              "method": "espD",
              "regions": [
                [990, 1581, "D"]
              ]
            },
            {
              "method": "espX",
              "regions": [
                [1, 17, "D"],
                [591, 592, "D"],
                ...
              ]
            },
            {
              "method": "glo",
              "regions": [
                [4, 4, "D"],
                [6, 6, "D"],
                ...
              ]
            },
            {
              "method": "dis465",
              "regions": [
                [1, 16, "D"],
                [51, 55, "D"],
                ...
              ]
            },
            {
              "method": "disHL",
              "regions": [
                [1, 13, "D"],
                [25, 36, "D"],
                ...
              ]
            },
            {
              "method": "vsl",
              "regions": [
                [1, 34, "D"],
                [48, 54, "D"],
                ...
              ]
            },
            {
              "method": "jronn",
              "regions": [
                [1, 21, "D"],
                [82, 98, "D"],
                ...
              ]
            },
            {
              "method": "seg",
              "regions": [
                [8, 23, "D"],
                [551, 574, "D"],
                ...
              ]
            },
            {
              "method": "pfilt",
              "regions": [
                [1094, 1109, "D"],
                [1245, 1256, "D"],
                ...
              ]
            }
          ]
        },
        "ss_populations": {
          "predictors": [
            {
              "method": "fess",
              "type": "helix",
              "scores": [
                0.004,
                0.115,
                0.18,
                ...
                0.233,
                0.083,
                0.003
              ]
            },
            {
              "method": "fess",
              "type": "sheet",
              "scores": [
                0.002,
                0.065,
                ...
                0.418,
                0.206,
                0.002
              ]
            },
            {
              "method": "fess",
              "type": "coil",
              "scores": [
                0.994,
                0.82,
                0.758,
                ...
                0.348,
                0.71,
                0.994
              ]
            },
            {
              "method": "dynamine",
              "type": "coil",
              "scores": [
                0.304,
                0.315,
                0.334,
                ...
                0.216,
                0.25,
                0.257
              ]
            }
          ]
        },
        "lips": {
          "predictors": [
            {
              "method": "anchor",
              "regions": [
                [17, 24, "D"],
                [61, 62, "D"],
                ...
              ]
            }
          ]
        }
      },
      "length": 1581,
      "accession": "sp|Q15648|MED1_HUMAN"
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
