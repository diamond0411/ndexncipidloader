===========================
NDEx NCI-PID content loader
===========================


.. image:: https://img.shields.io/pypi/v/ndexncipidloader.svg
        :target: https://pypi.python.org/pypi/ndexncipidloader

.. image:: https://img.shields.io/travis/ndexcontent/ndexncipidloader.svg
        :target: https://travis-ci.org/ndexcontent/ndexncipidloader

.. image:: https://coveralls.io/repos/github/ndexcontent/ndexncipidloader/badge.svg?branch=master
        :target: https://coveralls.io/github/ndexcontent/ndexncipidloader?branch=master

.. image:: https://readthedocs.org/projects/ndexncipidloader/badge/?version=latest
        :target: https://ndexncipidloader.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status


Python application that loads NCI-PID data into NDEx_

This tool downloads OWL_ files containing NCI-PID data from: ftp://ftp.ndexbio.org/NCI_PID_BIOPAX_2016-06-08-PC2v8-API/
and performs the following operations:

**1\)** OWL files are converted to extended SIF_ format using Paxtools_ and the SIF_ file is loaded intzo a network

**2\)** A node attribute named **type** is added to each node and is set to one of the following
   by extracting its value from **PARTICIPANT_TYPE** column in SIF_ file:

* **protein** (originally ProteinReference)

* **smallmolecule** (originally SmallMoleculeReference)

* **proteinfamily** (set if node name has **family** and was a **protein**)

* **RnaReference** (original value)

* **ProteinReference;SmallMoleculeReference** (original value)

**3\)** A node attribute named **alias** is added to each node and is loaded from **UNIFICATION_XREF**
column in SIF_ file which is split by `;` into a list. Each element of this list is prefixed with **uniprot:** and t first element is set as the
**represents** value in node and removed from the **alias** attribute. If after
removal, the **alias** attribute value is empty, it is removed.

**4\)** In SIF_ file **INTERACTION_TYPE** defines edge interaction type and **INTERACTION_PUBMED_ID** define
value of **citation** edge attribute. The values in **citation** edge attribute are
prefixed with **pubmed:** Once loaded redundant edges are removed
following these conventions:

* **neighbor-of** edges are removed if they contain no unique citations and an edge of another type exists

* **controls-state-change-map** edges are removed if they contain no unique citations and an edge of type other then **neighbor-of** exists

* **Special case:** After network has been updated following previous two conditions and there exists a **neighbor-of** edge with citations and **one** other edge exists with **no** citations, the citations from **neighbor-of** are added to the other edge and the **neighbor-of** edge is removed

**5\)** An edge attribute named **directed** is set to **True** if edge interaction type is one of the following (otherwise its set to **False**)

.. code-block::

    controls-state-change-of
    controls-transport-of
    controls-phosphorylation-of
    controls-expression-of
    catalysis-precedes
    controls-production-of
    controls-transport-of-chemical
    chemical-affects
    used-to-produce

**6\)** If node name matches **represents** value in node (with **uniprot:** prefix added) then the node name is replaced with gene symbol from `gene_symbol_mapping.json`_

**7\)** If node name starts with **CHEBI** then node name is replaced with value of **PARTICIPANT_NAME** from SIF_ column

**8\)** If node **represents** value starts with **chebi:CHEBI** the **chebi:** is removed

**9\)** If **_HUMAN** in SIF_ file **PARTICIPANT_NAME** column for a given node then this value is replaced by doing a lookup in `gene_symbol_mapping.json`_, unless value in lookup is **-** in which case original name is left

**10\)** Any node with **family** node name is changed as follows if a lookup of node name against **gene_symbol_mapping.json** returns one or more genes

* Node attribute named **member** is added and set to list of genes found in lookup in `gene_symbol_mapping.json`_
* Node attribute named **type** is changed to **proteinfamily**

**11\)** The following network attributes are set

* **name** set to name of OWL_ file with **.owl.gz** suffix removed except for **PathwayCommons.8.NCI_PID.BIOPAX** which is renamed to **NCI PID - Complete Interactions**
* **author** (from **Curated By** column in `networkattributes.tsv`_)
* **labels** (from **PID** column in `networkattributes.tsv`_)
* **organism** is pulled from **organism** attribute of `style.cx`_
* **prov:wasGeneratedBy** is set to html link to this repo with text ndexncipidloader <VERSION> (example: ndexncipidloader 1.2.0)
* **prov:wasDerivedFrom** is set to full path to OWL_ file on ftp site
* **reviewers** (from **Reviewed By** column in `networkattributes.tsv`_)
* **version** is set to Abbreviated month-year (example: MAY-2019)
* **description** is pulled from **description** attribute of `style.cx`_ except for **NCI PID - Complete Interactions** which has a hardcoded description set to `This network includes all interactions of the individual NCI-PID pathways.`
* **type** is set to list of string with single entry **pathway**
* **__normalizationversion** is set to 0.1

Dependencies
------------

* `ndex2 <https://pypi.org/project/ndex2>`_
* `ndexutil <https://pypi.org/project/ndexutil>`_
* `biothings_client <https://pypi.org/project/biothings-client>`_
* `requests <https://pypi.org/project/requests>`_
* `pandas <https://pypi.org/project/pandas>`_


Compatibility
-------------

* Python 3.3+

Installation
------------

.. code-block::

   git clone https://github.com/coleslaw481/ndexncipidloader
   cd ndexncipidloader
   make dist
   pip install dist/ndexncipidloader*whl


Configuration
-------------

The **loadndexncipidloader.py** requires a configuration file in the following format be created.
The default path for this configuration is :code:`~/.ndexutils.conf` but can be overridden with
:code:`--conf` flag.

**Format of configuration file**

.. code-block::

    [<value in --profile (default ndexncipidloader)>]

    user = <NDEx username>
    password = <NDEx password>
    server = <NDEx server(omit http) ie public.ndexbio.org>


**Example configuration file**

.. code-block::

    [ncipid_dev]

    user = joe123
    password = somepassword123
    server = dev.ndexbio.org


Required external tool
-----------------------

Paxtools is needed to convert the OWL files to SIF format.

Please download **paxtools.jar** (http://www.biopax.org/Paxtools/) (requires Java 8+) and
put in current working directory or specify path to **paxtools.jar** with `--paxtools` flag on
**loadnexncipidloader.py**

Usage
-----

For more information invoke :code:`loadndexncipidloader.py -h`

**Example usage**

This example assumes a valid configuration file with paxtools.jar in the working directory.

.. code-block::

   loadndexncipidloader.py sif

**Example usage with sif files already downloaded**

This example assumes a valid configuration file and the SIF files are located in :code:`sif/` directory

.. code-block::

   loadndexncipidloader.py --skipdownload sif

**All flags with descriptions**
positional arguments:
  sifdir                Directory containing .sif files to parse. Under this
                        directory OWL fileswill be downloaded and converted to
                        sif unless --skipdownload flag is set, in which case
                        this script assumes the*.sif files already exist

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  --profile PROFILE     Profile in configuration file to use to load NDEx
                        credentials which meansconfiguration under [XXX] will
                        beused (default ndexncipidloader)
  --logconf LOGCONF     Path to python logging configuration file in format
                        consumable by fileConfig. See
                        https://docs.python.org/3/library/logging.html for
                        more information. Setting this overrides -v|--verbose
                        parameter which uses default logger. (default None)
  --conf CONF           Configuration file to load (default ~/.ndexutils.conf
  --genesymbol GENESYMBOL
                        Use alternate gene symbol mapping file
  --loadplan LOADPLAN   Use alternate load plan file
  --networkattrib NETWORKATTRIB
                        Use alternate Tab delimited file containing PID
                        Pathway Name, reviewed by, curated by and revision
                        data for ncipid networks
  --style STYLE         Path to NDEx CX file to use for stylingnetworks
  --releaseversion RELEASEVERSION
                        Sets version network attribute (default current month
                        and year Example: JUL-2019)
  --singlefile SINGLEFILE
                        Only process file matching name in <sifdir>
  --paxtools PAXTOOLS   Path to paxtools.jar file used to convertowl file to
                        sif file. Ignored if --skipdownload flag is set.
                        Default assumespaxtools.jar is in current working
                        directory
  --skipdownload        If set, skips download of owl filesand conversion. The
                        program assumesthe <sifdir> directory set asthe last
                        argument on the command lineto this program contains
                        sif files
  --skipchecker         If set, skips gene symbol checker thatexamines all
                        nodes of type proteinand verifies they are symbols
  --disablcitededgemerge
                        If set, keeps neighbor-of edges if they contain
                        citations not found in moredescriptive edge
  --getfamilies         If set, code examines owl files and generates mapping
                        of protein families
  --ftphost FTPHOST     FTP host to download owl or sif files from. Ignored if
                        --skipdownload flag set (default ftp.ndexbio.org)
  --ftpdir FTPDIR       FTP directory to download owl or sif files from.
                        Ignored if --skipdownload flag set (default
                        NCI_PID_BIOPAX_2016-06-08-PC2v8-API)
  --verbose, -v         Increases verbosity of logger to standard error for
                        log messages in this module and in
                        ndexutil.tsv.tsv2nicecx2. Messages are output at these
                        python logging levels -v = ERROR, -vv = WARNING, -vvv
                        = INFO, -vvvv = DEBUG, -vvvvv = NOTSET (default is to
                        log CRITICAL)


Via Docker
~~~~~~~~~~~~~~~~~~~~~~

**Example usage**

This example **paxtools.jar** is in current directory, and a configuration
file has been created in current working directory and named :code:`conf`

.. code-block::

   docker run -v `pwd`:`pwd` -w `pwd` coleslawndex/ndexncipidloader:1.0.0 loadndexncipidloader.py --paxtools `pwd`/paxtools.jar --conf conf sif


Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. _NDEx: http://www.ndexbio.org
.. _OWL: https://en.wikipedia.org/wiki/Web_Ontology_Language
.. _Paxtools: https://www.biopax.org/Paxtools
.. _SIF: https://bioconductor.org/packages/release/bioc/vignettes/paxtoolsr/inst/doc/using_paxtoolsr.html#extended-simple-interaction-format-sif-network
.. _uniprot: https://www.uniprot.org/
.. _gene_symbol_mapping.json: https://github.com/ndexcontent/ndexncipidloader/blob/master/ndexncipidloader/gene_symbol_mapping.json
.. _networkattributes.tsv: https://github.com/ndexcontent/ndexncipidloader/blob/master/ndexncipidloader/networkattributes.tsv
.. _style.cx: https://github.com/ndexcontent/ndexncipidloader/blob/master/ndexncipidloader/style.cx
