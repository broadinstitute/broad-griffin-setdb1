version 1.0

import "../tasks/task-star-build-index.wdl" as build_index

workflow wf {
    call build_index.GenomeGenerate
}

