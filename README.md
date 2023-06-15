# flowplyr
Read and manipulate .fcs, .xml files associated with flow cytometry events for downstream analysis

### test data

* Test data for flowplyr can be accessed `/fg_data/ics_test/flowplyr_test_data`
* For link to test data `flowplyr_test_data.zip` (791.2 MB) contact maintainer.
* The minimal test data includes 27 specimen sample `.fcs` files and a 
`flowplyr_test_data.xml`

```
flowplyr_test_data
├── flowplyr_test_data.xml
└── flowplyr_test_fcs_folder
    ├── Specimen_001_A1_A01_029.fcs
    ├── Specimen_001_A2_A02_030.fcs
    ...
    ├── Specimen_001_A3_A03_031.fcs
    ├── Specimen_001_A4_A04_032.fcs
    └── Specimen_004_C4_C04_028.fcs
```

### preliminary identification of parent_gate and relevant markers

In order to run `extract_flow_events.R`, the user must specify a `parent_gate` and a list of `key_markers`. Only cells (i.e., flow events) within gates for one or more of the `key_markers` are retained in final tabular output, however, the output can include columns for fluorescent intensity for a larger set of `markers` that are not used for subsetting.

To identify an appropriate `parent_gate`, `markers`, and `key_markers`, one can first run the helper script `extract_marker_paths.R`.

#### run commandline `extract_marker_paths.R` on test data

```
TESTPATH=/fh/fast/gilbert_p/fg_data/ics_test/flowplyr_test_data
Rscript extract_marker_paths.R \
  --xml_path $TESTPATH/flowplyr_test_data.xml \
  --fcs_folder_path $TESTPATH/flowplyr_test_fcs_folder/
```

#### abbreviated output

This will load a FlowWorskpace gating_set object from the .xml and fcs files 
and provide a list of marker paths as well as CD4 and CD8 markers lists that 
can be used directly in the next step.

![example_gate_scheme](https://github.com/kmayerb/flowplyr/assets/46639063/77d091f8-a295-4050-8d62-104784bae282)

```
02:29:08 PM Gating Paths:
...
/Time/S/14-/Lv/L/3+
...
02:29:09 PM key CD8 Markers
c("8+/154+", "8+/GzB+", "8+/IFNg+", "8+/IFNg\\IL2", "8+/IL2+",
"8+/IL4+", "8+/IL17a+", "8+/TNFa+")
```


### specification of params.json

Once the `parent_gate` and `key_markers` are identified, to extract flow event data, 
create a params.json file. An example is provided `tests/test_params.json`



#### parameters

specify these parameters in a `param.json` file

* `batch_name`      - (character) name for tracking output 
* `xml_path`        - (character) file path flowJo .xml output
* `fcs_folder_path` - (character) file path to folder with .fcs files
* `parent_gate`     - (character) specifying the parent path for the analysis 
* `xml_keywords`    - (vector) keywords used for extracting relevant information for the .xml file
* `functional_markers`     - (vector) only cells (events) that fall into one of these marker gates will be retains
* `markers`         - (vector) potentially larger set of markers to include in the output
* `output_path`     - (character) path to write output files
* `write_csv`       - (logical) write output as .csv if TRUE
* `write_h5`        - (logical) write output as .h5 if TRUE


#### params.json for the test data

```{json}
{
  "batch_name": "test_flowplyr_batch",
  "xml_path": "/fh/fast/gilbert_p/fg_data/ics_test/flowplyr_test_data/flowplyr_test_data.xml",
  "fcs_folder_path": "/fh/fast/gilbert_p/fg_data/ics_test/flowplyr_test_data/flowplyr_test_fcs_folder/",
  "output_path": "/fh/fast/gilbert_p/fg_data/ics_test/gs_output/flowplyr_test_batch/",
  "parent_gate": "/Time/S/14-/Lv/L/3+",
  "write_csv": false,
  "write_h5": true,
  "xml_keywords": ["$FIL", "Stim", "Sample Order", "EXPERIMENT NAME", "Replicate"],
  "functional_markers" : ["8+/154+", "8+/GzB+", "8+/IFNg+", "8+/IFNg\\IL2", "8+/IL2+",
  "8+/IL4+", "8+/IL17a+", "8+/TNFa+"],
  "markers" : ["8+/154+", "8+/GzB+", "8+/IFNg+", "8+/IFNg\\IL2", "8+/IL2+",
  "8+/IL4+", "8+/IL17a+", "8+/TNFa+"]
}
```

### run commandline `extract_flow_events.R` on test data

```
Rscript extract_flow_events.R --params tests/test_params.json
```

### Outputs

The output if written in .h5 format has three component tables:

*`data/pos`
*`data/fi`
*`data/fcs_index`


```python
import h5py
import numpy as np
import pandas as pd
hdf5_file = '/fh/fast/gilbert_p/fg_data/ics_test/gs_output/flowplyr_test_batch/test_flowplyr_batch.h5'
fcs_ = pd.read_hdf(hdf5_file, key = "data/fcs_index")
hf = h5py.File(hdf5_file , 'r')
pos_ = hf.get('data/pos')
pos_ = np.array(pos_)
fi_ = hf.get('data/fi')
fi_ = np.array(fi_)
cols_pos = hf.get('data/cols_pos')
cols_fi = hf.get('data/cols_fi')
pos = pd.DataFrame(pos_.transpose(), columns = pd.Series(cols_pos).str.decode('UTF-8'))
fi = pd.DataFrame(fi_.transpose(), columns = pd.Series(cols_fi).str.decode('UTF-8'))
assert(fcs_.shape[0] == pos.shape[0])
assert(fcs_.shape[0] == fi.shape[0])
```


### Larger test
```
XML_PATH="/fh/fast/gilbert_p/fg_data/POI_COR/data/Aeras C040-404 BCG Correlates PoI Pilot/new_XMLs_assoc_wFCM_20220908/1814-R-C040 BCG CoR v2-FJ.xml"
FCS_PATH="/fh/fast/gilbert_p/fg_data/POI_COR/data/Aeras C040-404 BCG Correlates PoI Pilot/1814-R-C040-404 BCG"

TESTPATH=/fh/fast/gilbert_p/fg_data/ics_test/flowplyr_test_data
Rscript extract_marker_paths.R \
  --xml_path $TESTPATH/flowplyr_test_data.xml \
  --fcs_folder_path $TESTPATH/flowplyr_test_fcs_folder/

