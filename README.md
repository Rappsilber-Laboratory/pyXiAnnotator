# pyXiAnnotator

Python wrapper around [xiAnnotator](https://github.com/Rappsilber-Laboratory/xiAnnotator) MS spectra annotation for linear and crosslinked peptides with additional data analysis features.

## Requirements
The [xiAnnotator](https://github.com/Rappsilber-Laboratory/xiAnnotator) uses Java. The python wrapper requires following packages. 
```
requests
JPype1
memoized-property
numpy
pytest
python_version=3.7
```
## Build
Clone repository and build. The build method below is deprecated but still works.
```
git clone https://github.com/Rappsilber-Laboratory/pyXiAnnotator.git
cd pyXiAnnotator/
python setup.py sdist bdist_wheel
```
## Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install the wheel.

```bash
pip install dist/pyxiannotator-0.3.5-py3-none-any.whl
```

## Usage
Pseudo code example

```python
import pyxiannotator
import pandas as pd
import os
import pyteomics

# initialise pyxiannotator
annotator = pyxiannotator.XiAnnotator.XiAnnotatorLocal(java_home_dir='/path/to/jdk/', jvmargs="-Xmx800M")

# read file with crosslinked peptide spectrum matches
csm_df = pd.read_csv(csm_file)

# initialise mgf file readers
mgf_readers = {os.path.basename(f): pyteomics.mgf.read(f) for f in glob.glob('path/to/mgf/*.mgf')}


for i, psm in psm_df.iterrows():
    # get the spectrum id
    scan_id = psm[scan_id_column]

    # get the spectrum's peak list file
    mgf_file = psm[peak_list_file_name_column]

    # get the spectrum
    spectrum = mgf_readers[mgf_file].get_by_id(scan_id)

    # get the spectrum's peptide
    pep = pyxiannotator.AnnotatedSpectrum.Peptide(
        pep_seq1=psm['PepSeq1'],
        pep_seq2=psm['PepSeq2'],
        link_pos1=psm['LinkPos1'],
        link_pos2=psm['LinkPos2'],
    )

    peak_list = list(zip(spectrum['m/z array'], spectrum['intensity array']))

    # create the annotation request
    annotation_request = annotator.create_json_annotation_request(
        peak_list=peak_list,
        peptide=pep,
        precursor_charge=int(psm['Charge']),
        precursor_mz=float(psm['ExpMz']),
        fragment_types=fragment_types,
        fragment_tolerance_ppm=frag_ppm,
        fragment_tolerance_abs=frag_da,
        cross_linker=crosslinker,
        custom_settings=custom_settings,
        modifications=modifications
    )

    # annotate
    annotator.request_annotation_json(annotation_request)

    # retrieve the annotated spectrum 
    annotated_spectrum = annotator.get_annotated_spectrum()
```