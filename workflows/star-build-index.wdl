version 1.0

import "../tasks/star-build-index.wdl" as build_index

workflow wf {
    call build_index.GenomeGenerate
}

