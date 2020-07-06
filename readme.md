MobiDB-lite, long disorder consensus predictor
==============================================
Marco Necci, Damiano Piovesan, Zsuzsanna Dosztányi and Silvio C.E. Tosatto
Version 3.8.1

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
requires Java (version >=7.x) and Python (version 2.7.x). By default it uses
64-bit executables, but 32-bit version can be called by a command line option
(see below).

Usage
-----
The following command prints the help page on screen::

    python mobidb-lite.py -h


By default, MobiDB-lite searches the ``binx/`` folder (that includes predictors
executables) in its own directory. So, MobiDB-lite can be executed by just
providing a multi-fasta input::

    python mobidb-lite.py multifasta_example.fasta

MobiDB-Lite is able to use predictors executables from other paths. See **Configuration**
section for further details.

By default, MobiDB-lite runs on single thread but it can exploits up to 7 threads.
When less than 7 CPUs are available, the "-t" option allows to calibrate the load::

    python mobidb-lite.py multifasta_example.fasta -t 4

Configuration
-------------

MobiBD-Lite is designed to work out of the box. For a basic usage you should not
need to configure anything. Nonetheless, MobiDB-Lite offers some configuration
possibilities.

### Bin directories

MobiDB-Lite is able to use predictors executables from other paths. Paths of
single predictors can be configured in the ``config.ini`` or custom configuration file.
You should always use **absolute paths** when configuring binx directories.
Notice that when paths ar not set in this file, MobiDB-Lite will look for the
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

    python mobidb-lite.py multifasta_example.fasta -o out_file.txt


The output is provided in 3 different formats. The output format can be chosen
with the ``-f`` option, which accepts numbers from 0 to 3. The default short
version (option 0) includes the protein name and start/end position of the
disorder regions, one region per line. All elements are separated by a TAB::

    sp|P49137|MAPK2_HUMAN	1	43
    sp|P49137|MAPK2_HUMAN	10	43	Proline-rich
    sp|P36507|MP2K2_HUMAN	286	310
    sp|Q9H2G9|GO45_HUMAN	1	58
    sp|Q9H2G9|GO45_HUMAN	1	15	Polar
    sp|Q9GZV5|WWTR1_HUMAN	52	117
    sp|Q9GZV5|WWTR1_HUMAN	58	72	Polar
    sp|Q92734|TFG_HUMAN	187	400
    sp|Q92734|TFG_HUMAN	206	220	Polar
    sp|Q92734|TFG_HUMAN	227	347	Polar


Sequence features can be silenced by using option -sf, inactive by
default. When using this option, sequence features are skipped::

	sp|P49137|MAPK2_HUMAN	1	43
    sp|P36507|MP2K2_HUMAN	286	310
    sp|Q9H2G9|GO45_HUMAN	1	58
    sp|Q9GZV5|WWTR1_HUMAN	52	117
    sp|Q92734|TFG_HUMAN	187	400


Option 1 provides InterPro extended output (only MobiDB-Lite prediction, skip
fully structured proteins). Contains the binary MobiDB-Lite prediction and
the IDRs boundaries. Also contains a relaxed consensus (rlx_consensus) calculated
at a threshold of .375 (New in version 3.8.1) ::

    {
      "consensus": "DDD...SSS",
      "rlx_consensus": "DDD...SSS",
      "accession": "sp|P49137|MAPK2_HUMAN",
      "regions": [
        [1, 43, "D"]
      ]
    }


Option 2 provides an extended output and structured output that is used
to create annotations for MobiDB entries::

    {
      "mobidb_data": {
        "ss_populations": {
          "predictors": [
            {
              "type": "helix",
              "method": "fess",
              "regions": [0.004, 0.074, ...]
            },
            {...}
        ]
      },
        "disorder": {
          "predictors": [
            {
              "method": "espX",
              "regions": [
                [1, 60, "D"],
                [...]
              ]
            },
            {...}
          ]
        },
        "lips": {
          "predictors": [
            {
              "method": "anchor",
              "regions": [
                [1, 14, "D"],
                [...]
              ]
            }
          ]
        }
      },
      "length": 400,
      "mobidb_consensus": {
        "disorder": {
          "predictors": [
            {
              "method": "mobidb_lite",
              "dc": 0.1075,
              "scores": [0.875, 1.0, ...],
              "regions": [
                [1, 43, "D_WC"]
              ]
            },
            {
              "method": "simple",
              "regions": [
                [1, 46, "D"],
                [...]
              ]
            }
          ]
        }
      },
      "accession": "sp|P49137|MAPK2_HUMAN"
    }

Option 3 provides a very extended output, containing both scores and regions
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
contact us at biocompup@ngp-net.bio.unipd.it